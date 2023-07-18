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

// TODO: give bunny bones (add this 

void app_fbo() {
    char *vertex_shader_source = R""(
            #version 330 core

            uniform mat4 transform;
            uniform float time;
            layout (location = 0) in vec3 vertex_position;
            out vec3 w;

            mat3 rotation3dY(float angle) {
                float s = sin(angle);
                float c = cos(angle);
                return mat3(c, 0.0, -s, 0.0, 1.0, 0.0, s, 0.0, c);
            }

            void main() {
                gl_Position = transform * vec4(rotation3dY(1.0 * sin(0.02 * time) * vertex_position.y) * vertex_position, 1.0);
                w = vec3(0.0); 
                w[gl_VertexID % 3] = 1.0;
            }
        )"";
    char *fragment_shader_source = R""(
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
    Shader shader = shader_create(vertex_shader_source, fragment_shader_source);

    while (cow_begin_frame()) {
        static Camera3D camera = { 6.0 };
        camera_move(&camera);
        mat4 P = camera_get_P(&camera);
        mat4 V = camera_get_V(&camera);
        mat4 PV = P * V;
        static real time = 0.0;
        time += 0.0167;


        IndexedTriangleMesh3D *mesh = &library.meshes.bunny;
        int num_vertices       = mesh->num_vertices;
        vec3 *vertex_positions = mesh->vertex_positions;
        int num_triangles      = mesh->num_triangles;
        int3 *triangle_indices = mesh->triangle_indices;

        vec3 ray_origin = { -2.0, 0.0, 0.0 };
        static real a =  0.4;
        static real b = -0.4;
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
        mat4 ray_PV = _window_get_P_perspective(RAD(5)) * ray_V;

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
                        shader_set_uniform(&shader, "transform", ray_PV);
                        shader_set_uniform(&shader, "time", time);
                        shader_set_uniform(&shader, "IndexIfFalse_BarycentricWeightsIfTrue", false);
                        shader_pass_vertex_attribute(&shader, num_vertices, vertex_positions);
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
                tri = mesh->triangle_indices[triangleIndex];

                { // barycentric lookup
                    u8 rgb2[3];

                    glBindFramebuffer(GL_FRAMEBUFFER, fbo); {
                        glClearColor(1, 1, 1, 1);
                        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                        {
                            shader_set_uniform(&shader, "transform", ray_PV);
                            shader_set_uniform(&shader, "time", time);
                            shader_set_uniform(&shader, "IndexIfFalse_BarycentricWeightsIfTrue", true);

                            vec3 _vertex_positions[] = { mesh->vertex_positions[tri[0]], mesh->vertex_positions[tri[1]], mesh->vertex_positions[tri[2]], };
                            int3 _triangle_indices[] = { { 0, 1, 2, } };
                            shader_pass_vertex_attribute(&shader, 3, _vertex_positions);
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
            { // mesh->draw(P, V, globals.Identity);
                shader_set_uniform(&shader, "transform", PV);
                shader_set_uniform(&shader, "time", time);
                shader_set_uniform(&shader, "IndexIfFalse_BarycentricWeightsIfTrue", false);
                shader_pass_vertex_attribute(&shader, num_vertices, vertex_positions);
                shader_draw(&shader, num_triangles, triangle_indices);
            }
            if (triangleIndex != -1) {
                { // color triangle of interest
                    shader_set_uniform(&shader, "transform", PV);
                    shader_set_uniform(&shader, "time", time);
                    shader_set_uniform(&shader, "IndexIfFalse_BarycentricWeightsIfTrue", true);

                    vec3 _vertex_positions[] = { mesh->vertex_positions[tri[0]], mesh->vertex_positions[tri[1]], mesh->vertex_positions[tri[2]], };
                    int3 _triangle_indices[] = { { 0, 1, 2, } };
                    shader_pass_vertex_attribute(&shader, 3, _vertex_positions);
                    shader_draw(&shader, 1, _triangle_indices);
                }

                { // intersection point
                    // TODO: CPU bone stuff
                    vec3 p = {};
                    for_(d, 3) p += w[d] * vertex_positions[tri[d]];
                    draw_ball(P, V, p, monokai.green);
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
