#define _CRT_SECURE_NO_WARNINGS
typedef double real;
#define GL_REAL GL_DOUBLE
#define JIM_IS_JIM
#define SNAIL_I_SOLEMNLY_SWEAR_I_AM_UP_TO_NO_GOOD
#include "C:\Users\Jim\Documents\GitHub\CS3XX\include.cpp"

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

// TODO: extend to arbitrary ray

void app_fbo() {
    char *vertex_shader_source = R""(
            #version 330 core
            uniform mat4 transform;
            uniform float time;
            layout (location = 0) in vec3 vertex_position;

            mat3 rotation3dY(float angle) {
                float s = sin(angle);
                float c = cos(angle);
                return mat3(c, 0.0, -s, 0.0, 1.0, 0.0, s, 0.0, c);
            }

            void main() {
                gl_Position = transform * vec4(rotation3dY(2.0 * sin(0.05 * time) * vertex_position.y) * vertex_position, 1.0);
            }
        )"";
    char *fragment_shader_source = R""(
            #version 330 core
            precision highp float;
            out vec4 fragColor;
            uniform bool highlight;
            void main() {
                vec3 color = vec3(0.0, 1.0, 1.0);
                if (!highlight) {
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
        static Camera3D camera = { 3.0 };
        camera_move(&camera);
        mat4 P = camera_get_P(&camera);
        mat4 V = camera_get_V(&camera);
        mat4 PV = P * V;
        static real time = 0.0;
        time += 0.0167;

        IndexedTriangleMesh3D *mesh = &library.meshes.teapot;
        int num_vertices       = mesh->num_vertices;
        vec3 *vertex_positions = mesh->vertex_positions;
        int num_triangles      = mesh->num_triangles;
        int3 *triangle_indices = mesh->triangle_indices;

        int triangleIndex = -1; {
            u8 rgb[3]; {
                static Texture texture = texture_create("fbo");
                static u32 fbo = _fbo_create(texture);

                glBindFramebuffer(GL_FRAMEBUFFER, fbo); {
                    glClearColor(1, 1, 1, 1);
                    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                    {
                        shader_set_uniform(&shader, "transform", PV); // TODO
                        shader_set_uniform(&shader, "highlight", false);
                        shader_set_uniform(&shader, "time", time);
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

                glDisable(GL_DEPTH_TEST);
                library.meshes.square.draw(globals.Identity, globals.Identity, globals.Identity, {}, texture.name);
                glEnable(GL_DEPTH_TEST);
            }
            if (rgb[0] + rgb[1] + rgb[2] != 3 * 255) triangleIndex = rgb[0] + 256 * rgb[1] + 256 * 256 * rgb[2];
            gui_printf("triangleIndex = %d\n", triangleIndex);
        }

        if (triangleIndex != -1) { // color triangle of interest
            {
                shader_set_uniform(&shader, "transform", PV); // TODO
                shader_set_uniform(&shader, "highlight", true); // TODO
                shader_set_uniform(&shader, "time", time);
                shader_pass_vertex_attribute(&shader, num_vertices, vertex_positions);
                shader_draw(&shader, 1, &triangle_indices[triangleIndex]);
            }
            eso_begin(globals.Identity, SOUP_POINTS, 0, true);
            eso_color(monokai.yellow);
            eso_vertex(0.0, 0.0);
            eso_end();
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
