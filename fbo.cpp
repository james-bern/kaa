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
IntersectionResult GPU_pick(vec3 ray_origin, vec3 ray_direction, IndexedTriangleMesh3D *mesh) {
    int num_vertices       = mesh->num_vertices;
    vec3 *vertex_positions = mesh->vertex_positions;
    int num_triangles      = mesh->num_triangles;
    int3 *triangle_indices = mesh->triangle_indices;
    int num_bones          = mesh->num_bones;
    mat4 *bones            = mesh->bones;
    int4 *bone_indices     = mesh->bone_indices;
    vec4 *bone_weights     = mesh->bone_weights;

    bool DEBUG = false; // draws ray eye view

    static char *picking_vertex_shader_source = R""(
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
    static char *picking_fragment_shader_source = R""(
            #version 330 core
            precision highp float;

            uniform bool IndexIfFalse_BarycentricWeightsIfTrue;

            in  vec3 w;
            out vec4 fragColor;
            void main() {
                vec3 color = w;
                if (!IndexIfFalse_BarycentricWeightsIfTrue) {
                    int i = gl_PrimitiveID;
                    // for (int d = 0; d < 3; ++d) {
                    //     color[d] = (i % 256) / 255.0;
                    //     i /= 256;
                    // }
                    color[0] = (i % 256);
                    color[1] = ((i / 256) % 256);
                    color[2] = ((i / (256 * 256)) % 256);
                    color /= 255.0;
                }
                fragColor = vec4(color, 1.0);
            }
        )"";
    static Shader shader = shader_create(picking_vertex_shader_source, picking_fragment_shader_source);

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

    IntersectionResult result = {};
    {
        static Texture texture = texture_create("fbo");
        static u32 fbo = _fbo_create(texture);
        u8 rgb[3]; { // hit, tri
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

            if (DEBUG) {
                glDisable(GL_DEPTH_TEST);
                library.meshes.square.draw(M4_Translation(0.5, 0.5) * M4_Scaling(0.3), globals.Identity, globals.Identity, {}, texture.name);
                glEnable(GL_DEPTH_TEST);
            }
        }
        if (rgb[0] + rgb[1] + rgb[2] != 3 * 255) { // w, p
            result.hit = true;
            int triangleIndex = rgb[0] + 256 * rgb[1] + 256 * 256 * rgb[2];
            result.tri = triangle_indices[triangleIndex];


            glBindFramebuffer(GL_FRAMEBUFFER, fbo); {
                glClearColor(1, 1, 1, 1);
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                {
                    shader_set_uniform(&shader, "transform", ray_PV);
                    shader_set_uniform(&shader, "time", time);
                    shader_set_uniform(&shader, "IndexIfFalse_BarycentricWeightsIfTrue", true);

                    vec3 _vertex_positions[] = { vertex_positions[result.tri[0]], vertex_positions[result.tri[1]], vertex_positions[result.tri[2]], };
                    int4 _bone_indices[]     = {     bone_indices[result.tri[0]],     bone_indices[result.tri[1]],     bone_indices[result.tri[2]], };
                    vec4 _bone_weights[]     = {     bone_weights[result.tri[0]],     bone_weights[result.tri[1]],     bone_weights[result.tri[2]], };
                    int3 _triangle_indices[] = { { 0, 1, 2, } };
                    shader_pass_vertex_attribute(&shader, 3, _vertex_positions);
                    shader_pass_vertex_attribute(&shader, 3, _bone_indices);
                    shader_pass_vertex_attribute(&shader, 3, _bone_weights);
                    shader_draw(&shader, 1, _triangle_indices);
                }

                u8 rgb2[3]; {
                    real width, height;
                    _window_get_size(&width, &height);
                    glReadPixels(int(width / 2), int(height / 2), 1, 1, GL_RGB, GL_UNSIGNED_BYTE, rgb2);
                }
                for_(d, 3) result.w[d] = rgb2[d] / 255.0;

            } glBindFramebuffer(GL_FRAMEBUFFER, 0);

            result.p = mesh->_skin(result.tri, result.w);

            if (DEBUG) {
                glDisable(GL_DEPTH_TEST);
                library.meshes.square.draw(M4_Translation(0.5, -0.5) * M4_Scaling(0.3), globals.Identity, globals.Identity, {}, texture.name);
                glEnable(GL_DEPTH_TEST);
            }
        }
    }
    return result;
} 

void eg_fbo() {
    IndexedTriangleMesh3D mesh = library.meshes.bean;

    Camera3D camera = { 6.0 };
    real time = 0.0;
    while (cow_begin_frame()) {
        camera_move(&camera);
        mat4 P = camera_get_P(&camera);
        mat4 V = camera_get_V(&camera);
        mat4 PV = P * V;
        time += 0.0167;

        mesh.bones[0] = M4_RotationAboutXAxis(sin(0.1 * time));

        vec3 ray_origin = { -2.0, 0.0, 0.0 };
        vec3 ray_direction; {
            static real a = 0.1;
            static real b = 0.0;
            gui_slider("a", &a, -1.0, 1.0);
            gui_slider("b", &b, -1.0, 1.0);
            ray_direction = normalized(V3(1.0, a, b));
        }

        IntersectionResult result = GPU_pick(ray_origin, ray_direction, &mesh);

        mesh.draw(P, V, globals.Identity);
        draw_pipe(P, V, ray_origin, ray_origin + 2.0 * ray_direction, monokai.yellow, 0.5);
        if (result.hit) draw_ball(P, V, result.p);
    }
}

#if 0
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
        }
    }
}
#endif
