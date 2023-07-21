u32 _fbo_create(Texture texture) {
    unsigned int fbo;
    {
        glGenFramebuffers(1, &fbo);
        glBindFramebuffer(GL_FRAMEBUFFER, fbo);

        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture._texture_GLuint, 0);

        unsigned int rbo;
        glGenRenderbuffers(1, &rbo);
        glBindRenderbuffer(GL_RENDERBUFFER, rbo);
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, texture.width, texture.height);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, rbo);
        ASSERT(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);   
    }
    return fbo;
}

// TODO: give bunny bones (make this a builtin mesh in the library)
// TODO: IndexedTriangleMesh3D should have an underscored function to grab CPU transformed weights

void app_fbo() {
    char *picking_vertex_shader_source = R""(
            #version 330 core

            uniform mat4 transform;
            uniform float time;
            layout (location = 0) in  vec3 vertex;
            layout (location = 1) in ivec4 boneIndices;
            layout (location = 2) in  vec4 boneWeights;
            uniform mat4 bones[64];

            out vec3 w;

            void main() {
                vec4 tmp_position = vec4(vertex, 1.0);
                tmp_position = boneWeights.x * (bones[boneIndices.x] * tmp_position)
                             + boneWeights.y * (bones[boneIndices.y] * tmp_position)
                             + boneWeights.z * (bones[boneIndices.z] * tmp_position)
                             + boneWeights.w * (bones[boneIndices.w] * tmp_position);
                gl_Position = transform * tmp_position;
                w = vec3(0.0); 
                w[gl_VertexID % 3] = 1.0;
            }
        )"";
    char *picking_fragment_shader_source = R""(
            #version 330 core
            precision highp float;

            uniform bool IndexIfFalse_BarycentricWeightsIfTrue;

            in  vec3 w;
            out vec4 fragColor;
            void main() {
                vec3 color = w;
                if (!IndexIfFalse_BarycentricWeightsIfTrue) {
                    int i = gl_PrimitiveID;
                    for (int d = 0; d < 3; ++d) {
                        color[d] = (i % 256) / 255.0;
                        i /= 256;
                    }
                }
                fragColor = vec4(color, 1.0);
            }
        )"";
    Shader shader = shader_create(picking_vertex_shader_source, picking_fragment_shader_source);

    // TODO: use transform mesh (what you did messes up the normals) -- give better name PositionsAndNormals

    // IndexedTriangleMesh3D mesh = library.meshes.sphere;
    // {
    //     mat4 S = M4_Scaling(0.2, 0.5, 0.2);
    //     mat4 T = M4_Translation(0.0, 0.5, 0.0);
    //     mat4 M = T * S;
    //     for_(i, mesh.num_vertices) mesh.vertex_positions[i] = transformPoint(M, mesh.vertex_positions[i]);

    //     mesh.num_bones = 2;
    //     mesh.bones = (mat4 *) malloc(mesh.num_bones * sizeof(mat4));
    //     for_(i, mesh.num_bones) mesh.bones[i] = globals.Identity;
    //     mesh.bone_indices = (int4 *) malloc(mesh.num_vertices * sizeof(int4));
    //     for_(i, mesh.num_vertices) mesh.bone_indices[i] = { 0, 1, 0, 0 };
    //     mesh.bone_weights = (vec4 *) malloc(mesh.num_vertices * sizeof(vec4));
    //     for_(i, mesh.num_vertices) {
    //         real y = mesh.vertex_positions[i].y;
    //         mesh.bone_weights[i] = { 1.0 - y, y, 0.0, 0.0 };
    //     }
    // }
    IndexedTriangleMesh3D mesh = library.meshes.bean;
    int num_vertices       = mesh.num_vertices;
    vec3 *vertex_positions = mesh.vertex_positions;
    int num_triangles      = mesh.num_triangles;
    int3 *triangle_indices = mesh.triangle_indices;
    int num_bones          = mesh.num_bones;
    mat4 *bones            = mesh.bones;
    int4 *bone_indices     = mesh.bone_indices;
    vec4 *bone_weights     = mesh.bone_weights;
    // mesh._dump_for_library("tmp.txt", "bean");

    while (cow_begin_frame()) {
        static Camera3D camera = { 6.0 };
        camera_move(&camera);
        mat4 P = camera_get_P(&camera);
        mat4 V = camera_get_V(&camera);
        mat4 PV = P * V;
        static real time = 0.0;
        time += 0.0167;

        mesh.bones[0] = M4_RotationAboutXAxis(sin(0.1 * time));


        vec3 ray_origin = { -2.0, 0.0, 0.0 };
        static real a = 0.1;
        static real b = 0.0;
        gui_slider("a", &a, -1.0, 1.0);
        gui_slider("b", &b, -1.0, 1.0);
        vec3 ray_direction = normalized(V3(1.0, a, b));
        mat4 ray_V; {
            vec3 z = -ray_direction;
            vec3 Up = { 0.0, 1.0, 0.0 };
            vec3 y = cross(z, Up);
            vec3 x = cross(y, z);
            mat4 ray_C = M4_xyzo(x, y, z, ray_origin);
            // library.soups.axes.draw(PV * ray_C);
            ray_V = inverse(ray_C);
        }
        mat4 ray_PV = _window_get_P_perspective(RAD(3)) * ray_V;

        int triangleIndex = -1;
        Tri tri;
        vec3 w;
        if (1) {
            static Texture texture = texture_create("fbo");
            static u32 fbo = _fbo_create(texture);

            u8 rgb[3]; { // triangle index
                glBindFramebuffer(GL_FRAMEBUFFER, fbo); {
                    glClearColor(1, 1, 1, 1);
                    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                    {
                        shader_set_uniform(&shader, "bones", num_bones, bones);
                        shader_set_uniform(&shader, "transform", ray_PV);
                        shader_set_uniform(&shader, "time", time);
                        shader_set_uniform(&shader, "IndexIfFalse_BarycentricWeightsIfTrue", false);
                        shader_pass_vertex_attribute(&shader, num_vertices, vertex_positions);
                        shader_pass_vertex_attribute(&shader, num_vertices, bone_indices);
                        shader_pass_vertex_attribute(&shader, num_vertices, bone_weights);
                        shader_draw(&shader, num_triangles, triangle_indices);
                    }

                    {
                        real width, height;
                        _window_get_size(&width, &height);
                        glReadPixels(int(width / 2), int(height / 2), 1, 1, GL_RGB, GL_UNSIGNED_BYTE, rgb);
                    }
                } glBindFramebuffer(GL_FRAMEBUFFER, 0);

                gui_printf("r = %d\n", rgb[0]);
                gui_printf("g = %d\n", rgb[1]);
                gui_printf("b = %d\n", rgb[2]);

                if (1) {
                    glDisable(GL_DEPTH_TEST);
                    library.meshes.square.draw(M4_Translation(0.5, 0.5) * M4_Scaling(0.3), globals.Identity, globals.Identity, {}, texture.name);
                    glEnable(GL_DEPTH_TEST);
                }
            }

            if (rgb[0] + rgb[1] + rgb[2] != 3 * 255) {
                triangleIndex = rgb[0] + 256 * rgb[1] + 256 * 256 * rgb[2];
                tri = mesh.triangle_indices[triangleIndex];

                { // barycentric lookup
                    u8 rgb2[3];

                    glBindFramebuffer(GL_FRAMEBUFFER, fbo); {
                        glClearColor(1, 1, 1, 1);
                        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                        {
                            shader_set_uniform(&shader, "transform", ray_PV);
                            shader_set_uniform(&shader, "time", time);
                            shader_set_uniform(&shader, "IndexIfFalse_BarycentricWeightsIfTrue", true);

                            vec3 _vertex_positions[] = { mesh.vertex_positions[tri[0]], mesh.vertex_positions[tri[1]], mesh.vertex_positions[tri[2]], };
                            int4 _bone_indices[] = { mesh.bone_indices[tri[0]], mesh.bone_indices[tri[1]], mesh.bone_indices[tri[2]], };
                            vec4 _bone_weights[] = { mesh.bone_weights[tri[0]], mesh.bone_weights[tri[1]], mesh.bone_weights[tri[2]], };
                            int3 _triangle_indices[] = { { 0, 1, 2, } };
                            shader_pass_vertex_attribute(&shader, 3, _vertex_positions);
                            shader_pass_vertex_attribute(&shader, 3, _bone_indices);
                            shader_pass_vertex_attribute(&shader, 3, _bone_weights);
                            shader_draw(&shader, 1, _triangle_indices);
                        }

                        {
                            real width, height;
                            _window_get_size(&width, &height);
                            glReadPixels(int(width / 2), int(height / 2), 1, 1, GL_RGB, GL_UNSIGNED_BYTE, rgb2);
                        }

                        for_(d, 3) w[d] = rgb2[d] / 255.0;
                    } glBindFramebuffer(GL_FRAMEBUFFER, 0);

                    if (1) {
                        glDisable(GL_DEPTH_TEST);
                        library.meshes.square.draw(M4_Translation(0.5, -0.5) * M4_Scaling(0.3), globals.Identity, globals.Identity, {}, texture.name);
                        glEnable(GL_DEPTH_TEST);
                    }
                }
            }

            gui_printf("triangleIndex = %d\n", triangleIndex);
            gui_printf("alpha = %lf\n", w[0]);
            gui_printf("beta  = %lf\n", w[1]);
            gui_printf("gamma = %lf\n", w[2]);
        }

        { // drawing mesh
            { // mesh.draw(P, V, globals.Identity);
                shader_set_uniform(&shader, "bones", num_bones, bones);
                shader_set_uniform(&shader, "transform", PV);
                shader_set_uniform(&shader, "time", time);
                shader_set_uniform(&shader, "IndexIfFalse_BarycentricWeightsIfTrue", false);
                shader_pass_vertex_attribute(&shader, num_vertices, vertex_positions);
                shader_pass_vertex_attribute(&shader, num_vertices, bone_indices);
                shader_pass_vertex_attribute(&shader, num_vertices, bone_weights);
                shader_draw(&shader, num_triangles, triangle_indices);
            }
            if (triangleIndex != -1) {
                { // color triangle of interest
                    shader_set_uniform(&shader, "bones", num_bones, bones);
                    shader_set_uniform(&shader, "transform", PV);
                    shader_set_uniform(&shader, "time", time);
                    shader_set_uniform(&shader, "IndexIfFalse_BarycentricWeightsIfTrue", true);

                    vec3 _vertex_positions[] = { mesh.vertex_positions[tri[0]], mesh.vertex_positions[tri[1]], mesh.vertex_positions[tri[2]], };
                    int4 _bone_indices[] = { mesh.bone_indices[tri[0]], mesh.bone_indices[tri[1]], mesh.bone_indices[tri[2]], };
                    vec4 _bone_weights[] = { mesh.bone_weights[tri[0]], mesh.bone_weights[tri[1]], mesh.bone_weights[tri[2]], };
                    int3 _triangle_indices[] = { { 0, 1, 2, } };
                    shader_pass_vertex_attribute(&shader, 3, _vertex_positions);
                    shader_pass_vertex_attribute(&shader, 3, _bone_indices);
                    shader_pass_vertex_attribute(&shader, 3, _bone_weights);
                    shader_draw(&shader, 1, _triangle_indices);
                }

                { // intersection point
                    // TODO: CPU bone stuff
                    vec3 p = mesh._skin(tri, w);
                    draw_ball(P, V, p);
                }
            }
        }

        { // ray
            draw_pipe(P, V, ray_origin, ray_origin + 2.0 * ray_direction, monokai.yellow, 0.5);
        }
    }
}

void main() {
    APPS {
        // APP(eg_shader);
        // APP(eg_kitchen_sink);
        APP(app_fbo);
    }
}
