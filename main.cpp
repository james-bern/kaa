// TODO: GPU picking
// TODO: make line and spheres show up through the transparent mesh as well
// TODO: see if you can run cow while running an app in VR
// TODO: floor
// TODO: revisit the hessian (kim-style SparseMatrix parallel add?)
// TODO: talking to motors
// TODO: linear blend skinning in a vertex shader
//       get space fish back up and running
// TODO: #define in build.bat for gui stuff
// TODO: split IK line search between frames
// TODO: dFdu sparse matrix (why did this fail last time?)
// TODO: should x be a FixedSizeSelfDestructingArray<vec3>?
// // (waiting on josie) real-world trajectory Figure
// play more with sim params
// play more with ik weights
// // visualization
// pipe and sphere widget
// cable slack visualization
// // fun
// ordering juggling clubs
// // cow
// port MIN, MAX, etc. to be functions

#include "include.cpp"
bool DRAGON = false;

void draw_ball(mat4 P, mat4 V, vec3 s, vec3 color = monokai.white, real scale = 1.0) { library.meshes.sphere.draw(P, V, M4_Translation(s) * M4_Scaling(scale * 0.01), color); };
void draw_pipe(mat4 P, mat4 V, vec3 s, vec3 t, vec3 color = monokai.white, real scale = 1.0) {
    real r = scale * 0.007;
    vec3 e = t - s;
    real mag_e = norm(e);
    vec3 E = { 0.0, 1.0, 0.0 };
    vec3 e_hat = e / mag_e;
    vec3 u = cross(e_hat, E);
    real norm_u = norm(u);
    mat4 R = (norm_u < 1.0e-7) ? globals.Identity : M4_RotationAxisAngle(u / norm_u, ((dot(e_hat, E) < 0) ? 1.0 : -1.0) * asin(norm_u));
    library.meshes.cylinder.draw(P, V, M4_Translation(0.5 * (s + t)) * R * M4_Scaling(r, mag_e, r) * M4_Translation(0.0, -0.5, 0.0), color);
};

#include "widget.cpp"

struct IntersectionResult {
    bool hit;
    vec3 p;
    int3 tri;
    vec3 w;
    real t; // _NOT_ set by GPU version of picking
};
IntersectionResult ray_triangle_intersection(vec3 ray_origin, vec3 ray_direction, Tri tri, const SDVector &x) {
    vec3 a = get(x, tri[0]);
    vec3 b = get(x, tri[1]);
    vec3 c = get(x, tri[2]);
    //                                          p = p          
    // alpha * a + beta * b + gamma * c           = ray_origin + t * ray_direction
    // alpha * a + beta * b + gamma * c - t * ray_direction = ray_origin          
    //                            [ a b c -ray_direction] w = ray_origin          
    //                                        F w = ray_origin          
    vec4 w__t = inverse(M4(a.x, b.x, c.x, -ray_direction.x, a.y, b.y, c.y, -ray_direction.y, a.z, b.z, c.z, -ray_direction.z, 1.0, 1.0, 1.0, 0.0)) * V4(ray_origin.x, ray_origin.y, ray_origin.z, 1);
    bool hit; {
        hit = true;
        for (int k = 0; k < 4; ++k) hit &= w__t[k] > -TINY_VAL;
    }
    IntersectionResult result = {};
    result.tri = tri; // ! (have to copy this over so it makes sense with the ray_mesh_intersection API)
    result.hit = hit;
    if (result.hit) {
        for_(d, 3) result.w[d] = w__t[d];
        result.t = w__t[3];
        result.p = ray_origin + result.t * ray_direction;
    }
    return result;
}
IntersectionResult ray_mesh_intersection(vec3 ray_origin, vec3 ray_direction, const SDVector &x, int num_triangles, int3 *triangle_indices) {
    IntersectionResult result = {}; {
        real min_t = INFINITY;
        for_(triangle_i, num_triangles) {
            Tri tri = triangle_indices[triangle_i];
            IntersectionResult singleRayResult = ray_triangle_intersection(ray_origin, ray_direction, tri, x);
            if (singleRayResult.hit && (singleRayResult.t < min_t)) {
                min_t = singleRayResult.t;
                result = singleRayResult;
            }
        }
    }
    return result;
}


#include "fbo.cpp"





// bones


// C++
// - IK

// C#
// - controllers
// - draw
// -- (GPU)
// --- "skinning" - vertex shader


bool AUTOMATED_SPEED_TEST__QUITS_AFTER_A_COUPLE_SECONDS = false;

const int  MESH_NUMBER_OF_VOLUMETRIC_STACKS_PER_UPPER_SEGMENT = 3;
const int  MESH_NUMBER_OF_VOLUMETRIC_STACKS_PER_LOWER_SEGMENT = 2;
const bool  INCLUDE_DUMMY_SEGMENT                             = false; // FORNOW: bottom segment always 1 stack

int IK_MAX_LINE_SEARCH_STEPS = 8;


const real ROBOT_SEGMENT_LENGTH = 0.1450;
const real ROBOT_SEGMENT_RADIUS = 0.06 / 2;
const int  ROBOT_NUMBER_OF_UPPER_SEGMENTS = 1;
const int  ROBOT_NUMBER_OF_LOWER_SEGMENTS = 4;
const int  ROBOT_NUMBER_OF_SEGMENTS = ROBOT_NUMBER_OF_UPPER_SEGMENTS + ROBOT_NUMBER_OF_LOWER_SEGMENTS + (INCLUDE_DUMMY_SEGMENT ? 1 : 0);
const real ROBOT_LENGTH = ROBOT_NUMBER_OF_SEGMENTS * ROBOT_SEGMENT_LENGTH;
const real MESH_UPPER_STACK_LENGTH = ROBOT_SEGMENT_LENGTH / MESH_NUMBER_OF_VOLUMETRIC_STACKS_PER_UPPER_SEGMENT;
const real MESH_LOWER_STACK_LENGTH = ROBOT_SEGMENT_LENGTH / MESH_NUMBER_OF_VOLUMETRIC_STACKS_PER_LOWER_SEGMENT;
const int  MESH_NUMBER_OF_ANGULAR_SECTIONS = 9;
const int  MESH_NUMBER_OF_NODES_PER_NODE_LAYER = 1 + MESH_NUMBER_OF_ANGULAR_SECTIONS;
const int  _MESH_NUMBER_OF_UPPER_NODE_LAYERS_EXCLUSIVE = ROBOT_NUMBER_OF_UPPER_SEGMENTS * MESH_NUMBER_OF_VOLUMETRIC_STACKS_PER_UPPER_SEGMENT;
const int  _MESH_NUMBER_OF_LOWER_NODE_LAYERS_EXCLUSIVE = ROBOT_NUMBER_OF_LOWER_SEGMENTS * MESH_NUMBER_OF_VOLUMETRIC_STACKS_PER_LOWER_SEGMENT;
const int  MESH_NUMBER_OF_NODE_LAYERS = 1 + _MESH_NUMBER_OF_UPPER_NODE_LAYERS_EXCLUSIVE + _MESH_NUMBER_OF_LOWER_NODE_LAYERS_EXCLUSIVE + (INCLUDE_DUMMY_SEGMENT ? 1 : 0);

////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////

Sim sim;
IndexedTriangleMesh3D dragonHead = _meshutil_indexed_triangle_mesh_load("head.obj", false, true, false);
IndexedTriangleMesh3D dragonBody = _meshutil_indexed_triangle_mesh_load("body.obj", false, true, false);
SDVector u_MAX;
State currentState;
#define MAX_NUM_FEATURE_POINTS 16
Via featurePoints[MAX_NUM_FEATURE_POINTS];

////////////////////////////////////////////////////////////////////////////////

#define DLL_EXPORT extern "C" __declspec(dllexport)
#define delegate DLL_EXPORT
typedef float    UnityVertexAttributeFloat;
typedef uint16_t UnityTriangleIndexInt;
typedef int      UnityGeneralPurposeInt;

////////////////////////////////////////////////////////////////////////////////

#define LEN_U sim.num_cables
#define LEN_X (SOFT_ROBOT_DIM * sim.num_nodes)

delegate UnityGeneralPurposeInt cpp_getNumVertices() { return sim.num_nodes; }
delegate UnityGeneralPurposeInt cpp_getNumTriangles() { return sim.num_triangles; }
bool initialized;


// resets the CPP stuff (featurePoints, simultaion state)
// the caller is responsible for resetting targetPositions and targetEnabled
delegate void cpp_reset() {
    memset(featurePoints, 0, sizeof(featurePoints));
    featurePoints[0] = { { sim.num_nodes - 1, 1.0 } };
    currentState.x = sim.x_rest;
    currentState.u.setZero();
}

delegate void cpp_init() {
    ASSERT(!initialized);
    initialized = true;

    sim = {}; {
        // slice    
        //     2    
        // 3       1
        //     n    
        // 5       0
        //    ...   


        StretchyBuffer<vec3> X = {}; {
            real length = 0.0;
            for_(i, MESH_NUMBER_OF_NODE_LAYERS) {
                for_(a, MESH_NUMBER_OF_ANGULAR_SECTIONS) {
                    real angle = NUM_DEN(a, MESH_NUMBER_OF_ANGULAR_SECTIONS) * TAU;
                    sbuff_push_back(&X, V3(ROBOT_SEGMENT_RADIUS * cos(angle), -length, ROBOT_SEGMENT_RADIUS * sin(angle)));
                }
                sbuff_push_back(&X, V3(0.0, -length, 0.0));
                if (i == MESH_NUMBER_OF_NODE_LAYERS - 1) ASSERT(ARE_EQUAL(length, ROBOT_LENGTH));

                if (i < _MESH_NUMBER_OF_UPPER_NODE_LAYERS_EXCLUSIVE) {
                    length += MESH_UPPER_STACK_LENGTH;
                } else if (i < _MESH_NUMBER_OF_UPPER_NODE_LAYERS_EXCLUSIVE + _MESH_NUMBER_OF_LOWER_NODE_LAYERS_EXCLUSIVE) {
                    length += MESH_LOWER_STACK_LENGTH;
                } else {
                    length += ROBOT_SEGMENT_LENGTH;
                }
            }
        }

        StretchyBuffer<Tet> tets = {}; {
            // https://www.alecjacobson.com/weblog/media/triangular-prism-split-into-three-tetrahedra.pdf
            for_(_c, MESH_NUMBER_OF_NODE_LAYERS - 1) {
                int c = (1 + _c) * MESH_NUMBER_OF_NODES_PER_NODE_LAYER - 1; // center of slice
                int f = c + MESH_NUMBER_OF_NODES_PER_NODE_LAYER;
                int o = _c * MESH_NUMBER_OF_NODES_PER_NODE_LAYER;
                for_line_loop_(_a, _b, MESH_NUMBER_OF_ANGULAR_SECTIONS) {
                    int a = o + _a;
                    int b = o + _b;
                    int d = a + MESH_NUMBER_OF_NODES_PER_NODE_LAYER;
                    int e = b + MESH_NUMBER_OF_NODES_PER_NODE_LAYER;
                    sbuff_push_back(&tets, { a, b, c, d       });
                    sbuff_push_back(&tets, {    b, c, d, e    });
                    sbuff_push_back(&tets, {       c, d, e, f });
                }
            }
        }

        StretchyBuffer<int> pins = {}; {
            for_(i, MESH_NUMBER_OF_NODES_PER_NODE_LAYER) { sbuff_push_back(&pins, i); }
        }

        StretchyBuffer<int> num_vias = {};
        StretchyBuffer<Via> vias = {};
        {
            for_(cable_group, 3) {

                for_(d, 3) {

                    int cable_num_vias = (cable_group == 0) ?
                        (1 + MESH_NUMBER_OF_VOLUMETRIC_STACKS_PER_UPPER_SEGMENT) :
                        (1 + 2 * MESH_NUMBER_OF_VOLUMETRIC_STACKS_PER_LOWER_SEGMENT);

                    sbuff_push_back(&num_vias, cable_num_vias);

                    int i_0; {
                        i_0 = cable_group + 3 * d;
                        if (cable_group >= 1) i_0 +=     (MESH_NUMBER_OF_VOLUMETRIC_STACKS_PER_UPPER_SEGMENT * MESH_NUMBER_OF_NODES_PER_NODE_LAYER);
                        if (cable_group >= 2) i_0 += 2 * (MESH_NUMBER_OF_VOLUMETRIC_STACKS_PER_LOWER_SEGMENT * MESH_NUMBER_OF_NODES_PER_NODE_LAYER);
                    }

                    for_(k, cable_num_vias) {
                        sbuff_push_back(&vias, { { i_0 + (k * MESH_NUMBER_OF_NODES_PER_NODE_LAYER), 1.0 } });
                    }
                }
            }
        }

        SimInput simInput = {}; {
            simInput.num_nodes = X.length;
            simInput.x_rest = (real *) X.data;

            simInput.num_tets = tets.length;
            simInput.tets = tets.data;

            simInput.num_pins = pins.length;
            simInput.pins = pins.data;

            simInput.num_vias = num_vias.data;
            simInput.num_cables = num_vias.length;
            simInput.vias = vias.data;
        }

        sim.allocAndPrecompute(&simInput);
    }
    u_MAX = SDVector(sim.num_cables); {
        for_(j, sim.num_cables) {
            u_MAX[j] = 0.8 * ROBOT_SEGMENT_LENGTH;
            if (j > 3) u_MAX[j] *= 2;

        }
    }
    currentState = State(&sim);
    currentState.enabled__BIT_FIELD = TETS | GRAVITY | PINS | CABLES;

    sim.getNext(&currentState); // FORNOW: global_U_xx

    cpp_reset();
}



void jones() {
    cpp_init();

    Camera3D camera = { 1.7 * ROBOT_LENGTH, RAD(60), 0.0, 0.0, 0.0, -0.5 * ROBOT_LENGTH };
    while (cow_begin_frame()) {
        camera_move(&camera);
        mat4 P = camera_get_P(&camera);
        mat4 V = camera_get_V(&camera);
        mat4 PV = P * V;

        if (gui_button("reset", 'r')) cpp_reset();


        { // parameter sliders
            gui_slider("tetMassDensity", &tetMassDensity, 0.0, 10000.0);
            gui_slider("tetYoungsModulus", &tetYoungsModulus, 2.0, 8.0, false, true);
            gui_slider("tetPoissonsRatio", &tetPoissonsRatio, 0.01, 0.49);
        }
        if (gui_button("disable / enable cables", 'c')) currentState.enabled__BIT_FIELD ^= CABLES;
        if (gui_button("disable cables", 'r')) cpp_reset();


        // TODO: we assume that Sadra "assembles" / sets each motor so its cable is just barely taught
        // - TODO: modify simecableReferenceLengths


        { // control sliders
            for_(j, sim.num_cables) {
                char buffer[] = "u_?";
                buffer[strlen(buffer) - 1] = char('0' + j);
                real a = (j < 3) ? ROBOT_SEGMENT_LENGTH : 2 * ROBOT_SEGMENT_LENGTH;
                gui_slider(buffer, &currentState.u[j], -a, a);
            }
        }
        currentState = sim.getNext(&currentState);

        sim.draw(P * V, &currentState);
        // sim.draw(P, V, globals.Identity, &currentState);
    }
}



// returns whether the ray hit the mesh
// if the ray hit the mesh, functions writes to intersection_position__FLOAT_ARRAY__LENGTH_3
// if the ray hit the mesh and pleaseSetFeaturePoint, updates indexOfFeaturePointToSet-th feature point in C++
//                                                    and writes to feature_point_positions__FLOAT3__ARRAY (unless NULL)
delegate bool cpp_castRay(
        float ray_origin_x,
        float ray_origin_y,
        float ray_origin_z,
        float ray_direction_x,
        float ray_direction_y,
        float ray_direction_z,
        void *intersection_position__FLOAT_ARRAY__LENGTH_3,
        bool pleaseSetFeaturePoint,
        int indexOfFeaturePointToSet,
        void *feature_point_positions__FLOAT3__ARRAY = NULL) {

    vec3 ray_origin = { ray_origin_x, ray_origin_y, ray_origin_z };
    vec3 ray_direction = { ray_direction_x, ray_direction_y, ray_direction_z };
    IntersectionResult result;
    if (DRAGON) {
        result = GPU_pick(ray_origin, ray_direction, &dragonBody);
    } else {
        result = ray_mesh_intersection(ray_origin, ray_direction, currentState.x, sim.num_triangles, sim.triangle_indices);
    }
    if (result.hit) for_(d, 3) (((float *) intersection_position__FLOAT_ARRAY__LENGTH_3)[d]) = float(result.p[d]);
    if (pleaseSetFeaturePoint) {
        ASSERT(indexOfFeaturePointToSet >= 0);
        ASSERT(indexOfFeaturePointToSet < MAX_NUM_FEATURE_POINTS);
        featurePoints[indexOfFeaturePointToSet] = { { { result.tri[0], result.w[0] }, { result.tri[1], result.w[1] }, { result.tri[2], result.w[2] } } };
        if (feature_point_positions__FLOAT3__ARRAY) {
            vec3 tmp = get(currentState.x, featurePoints[indexOfFeaturePointToSet]);
            for_(d, 3) ((float *) feature_point_positions__FLOAT3__ARRAY)[3 * indexOfFeaturePointToSet + d] = float(tmp[d]);
        }
    }

    return result.hit;
}

// solves one step of IK (and physics)
// writes resulting mesh to vertex_positions__FLOAT3_ARRAY, vertex_normals__FLOAT3_ARRAY, triangle_indices__UINT_ARRAY
// also writes feature_point_positions__FLOAT3__ARRAY for use by unity
delegate void cpp_solve(
        int num_feature_points,
        void *_targetEnabled__INT_ARRAY,
        void *_targetPositions__FLOAT3_ARRAY,
        void *vertex_positions__FLOAT3_ARRAY,
        void *vertex_normals__FLOAT3_ARRAY,
        void *triangle_indices__UINT_ARRAY,
        void *feature_point_positions__FLOAT3__ARRAY) {
    ASSERT(num_feature_points <= MAX_NUM_FEATURE_POINTS);
    int *targetEnabled = (int *) _targetEnabled__INT_ARRAY;
    vec3 targetPositions[MAX_NUM_FEATURE_POINTS]; {
        for_(i, num_feature_points) {
            for_(d, 3) targetPositions[i][d] = ((float *) _targetPositions__FLOAT3_ARRAY)[3 * i + d];
        }
    }
    bool relax; {
        relax = true;
        for_(i, MAX_NUM_FEATURE_POINTS) relax &= (!targetEnabled[i]);
    }

    { // step ik
        static real Q_c      = 1.0;
        static real R_c_log  = 0.0036;
        static real R_c_quad = 0.0036;
        static real S_c      = 7.5;
        static real S_eps    = 0.002;
        static real alpha_0  = 0.0045;
        static real _R_c_RELAX  = 32.0;
        static bool project = true;

        // {
        //     gui_slider("Q_c", &Q_c, 0.0, 10.0);
        //     gui_slider("R_c", &R_c, 0.0, 0.1);
        //     gui_slider("S_c", &S_c, 0.0, 10.0);
        //     gui_slider("S_eps", &S_eps, 0.0, 0.01);
        //     gui_slider("alpha_0", &alpha_0, 0.0, 0.01);
        // }
        // gui_checkbox("project", &project, 'b'); // FORNOW (will break dll?)

        auto get_O = [&](State staticallyStableState) -> real {
            SDVector u = staticallyStableState.u;
            SDVector x = staticallyStableState.x;

            if (relax) {
                return _R_c_RELAX * (0.5 * squaredNorm(u));
            } else {
                real Q = 0.0; {
                    for_(i, MAX_NUM_FEATURE_POINTS) if (targetEnabled[i]) Q += Q_c * .5 * squaredNorm(get(x, featurePoints[i]) - targetPositions[i]);
                }

                real R = 0.0; {
                    R += R_c_quad * (0.5 * squaredNorm(u));
                    for_(j, sim.num_cables) {
                        R += R_c_log  * (-log(u_MAX[j] - u[j]));
                        if (u[j] > u_MAX[j]) { R = HUGE_VAL; }
                    }

                }

                real S = 0.0; {
                    SDVector Deltas = sim.computeCableDeltas(x, u);
                    for_(j, LEN_U) {
                        S += ZCQ(S_c, S_eps - Deltas[j]);
                    }

                }

                real result = Q + R + S;
                if (_isnan(result)) { result = INFINITY; }
                return result;
            }
        };

        { // single gradient descent step with backtracking line search
            SDVector dOdu(LEN_U);
            SDVector &u = currentState.u;
            const SDVector &x = currentState.x;
            {
                if (project) { // project
                    SDVector Slacks = sim.computeCableDeltas(x, u);
                    for_(j, LEN_U) Slacks[j] = -MIN(0.0, Slacks[j]);
                    for_(j, LEN_U) {
                        if (Slacks[j] > -.000001) {
                            u[j] += Slacks[j] + .00001;
                        }
                    }
                }
                if (relax) {
                    // TODO scalar multiplication
                    for_(j, sim.num_cables) {
                        dOdu[j] += _R_c_RELAX * u[j];
                    }
                } else {
                    SDVector dQdu; {
                        SDVector dQdx(LEN_X); {
                            for_(i, MAX_NUM_FEATURE_POINTS) if (targetEnabled[i]) add(dQdx, featurePoints[i], Q_c * get(x, featurePoints[i]) - targetPositions[i]);
                        }
                        SDVector L = global_U_xx.solve(dQdx);
                        dQdu = L * global_dFdu;
                    }

                    SDVector dRdu(LEN_U); {
                        for_(j, sim.num_cables) {
                            dRdu[j] += R_c_quad * u[j];
                            dRdu[j] += R_c_log  / (u_MAX[j] - u[j]);
                            if (u[j] > u_MAX[j]) { dRdu[j] = INFINITY; }
                        }
                    }

                    SDVector dSdu(LEN_U); {
                        SDVector Deltas = sim.computeCableDeltas(x, u);
                        for_(j, sim.num_cables) {
                            dSdu[j] += -ZCQp(S_c, S_eps - Deltas[j]);
                        }
                    }

                    for_(j, LEN_U) dOdu[j] = dQdu[j] + dRdu[j] + dSdu[j];
                    for_(j, LEN_U) if (_isnan(dOdu[j])) dOdu[j] = 0; // FORNOW
                }
            }
            {
                real O_curr = get_O(currentState);
                State nextState = currentState;

                int attempt = 0;
                do {
                    for_(j, LEN_U) nextState.u[j] = currentState.u[j] - alpha_0 * pow(.5, attempt) * dOdu[j];
                    nextState = sim.getNext(&nextState);
                    real O_next = get_O(nextState);
                    if (O_next < O_curr) { break; }
                } while (attempt++ < IK_MAX_LINE_SEARCH_STEPS);

                currentState = nextState;
            }
        }
    }

    { // marshall mesh
        SDVector vertex_normals = sim.get_vertex_normals(currentState.x);
        for_(k, LEN_X) {
            ((UnityVertexAttributeFloat *) vertex_positions__FLOAT3_ARRAY)[k] = UnityVertexAttributeFloat(currentState.x[k]);
            ((UnityVertexAttributeFloat *) vertex_normals__FLOAT3_ARRAY)[k] = UnityVertexAttributeFloat(vertex_normals[k]);
        }
        for_(k, 3 * sim.num_triangles) ((UnityTriangleIndexInt *) triangle_indices__UINT_ARRAY)[k] = (UnityTriangleIndexInt) ((int *) sim.triangle_indices)[k];
        for_(i, num_feature_points) {
            vec3 tmp = get(currentState.x, featurePoints[i]);
            for_(d, 3) ((float *) feature_point_positions__FLOAT3__ARRAY)[3 * i + d] = float(tmp[d]);
        }
    }
}

vec3 SPOOF_targetPositions[MAX_NUM_FEATURE_POINTS];
int  SPOOF_targetEnabled[MAX_NUM_FEATURE_POINTS];

void SPOOF_reset() {
    memset(SPOOF_targetPositions, 0, sizeof(SPOOF_targetPositions));
    memset(SPOOF_targetEnabled, 0, sizeof(SPOOF_targetEnabled));
    SPOOF_targetPositions[0] = get(sim.x_rest, featurePoints[0]) + .3 * V3(0, 1, 1);
    SPOOF_targetEnabled[0] = TRUE;
}


void kaa() {
    cpp_init();
    SPOOF_reset();


    UnityVertexAttributeFloat *SPOOF_vertex_positions = (UnityVertexAttributeFloat *) calloc(LEN_X, sizeof(UnityVertexAttributeFloat));
    UnityVertexAttributeFloat *SPOOF_vertex_normals = (UnityVertexAttributeFloat *) calloc(LEN_X, sizeof(UnityVertexAttributeFloat));
    UnityTriangleIndexInt     *SPOOF_triangle_indices = (UnityTriangleIndexInt *) calloc(3 * sim.num_triangles, sizeof(UnityTriangleIndexInt));
    UnityVertexAttributeFloat *SPOOF_feature_point_positions = (UnityVertexAttributeFloat *) calloc(3 * MAX_NUM_FEATURE_POINTS, sizeof(UnityVertexAttributeFloat));

    COW1._gui_hide_and_disable = AUTOMATED_SPEED_TEST__QUITS_AFTER_A_COUPLE_SECONDS;
    Camera3D camera = { 1.7 * ROBOT_LENGTH, RAD(60), 0.0, 0.0, 0.0, -0.5 * ROBOT_LENGTH };
    while (cow_begin_frame()) {
        camera_move(&camera);
        mat4 P = camera_get_P(&camera);
        mat4 V = camera_get_V(&camera);
        mat4 PV = P * V;


        auto draw_sphere = [&](vec3 s, vec3 color = monokai.white) { library.meshes.sphere.draw(P, V, M4_Translation(s) * M4_Scaling(0.01), color); };

        struct CastRayResult {
            bool intersects;
            vec3 intersection_position;
        };
        auto castRay = [&](bool pleaseSetFeaturePoint, int featurePointIndex) -> CastRayResult {
            CastRayResult result = {};
            vec3 ray_origin = camera_get_position(&camera);
            vec3 ray_direction = camera_get_mouse_ray(&camera);
            float intersection_position__FLOAT_ARRAY__LENGTH_3[3]; {
                result.intersects = cpp_castRay(
                        (float) ray_origin.x,
                        (float) ray_origin.y,
                        (float) ray_origin.z,
                        (float) ray_direction.x,
                        (float) ray_direction.y,
                        (float) ray_direction.z,
                        intersection_position__FLOAT_ARRAY__LENGTH_3,
                        pleaseSetFeaturePoint,
                        featurePointIndex,
                        SPOOF_feature_point_positions);
            }
            for_(d, 3) result.intersection_position[d] = intersection_position__FLOAT_ARRAY__LENGTH_3[d];
            return result;
        };

        if (gui_button("reset", 'r')) {
            cpp_reset();
            SPOOF_reset();
        }

        // NOTE: very important to have solved physics at least once so the global hessian is ret to go
        static bool SPOOF_solveIK = true;
        gui_checkbox("SPOOF_solveIK", &SPOOF_solveIK, COW_KEY_SPACE);
        if (SPOOF_solveIK) { // ik
            if ((sim.num_cables >= 0)) {
                float _SPOOF_target_positions__FLOAT_ARRAY[3 * MAX_NUM_FEATURE_POINTS]; {
                    for_(k, _COUNT_OF(_SPOOF_target_positions__FLOAT_ARRAY)) _SPOOF_target_positions__FLOAT_ARRAY[k] = float(((real *) SPOOF_targetPositions)[k]);
                }

                cpp_solve(
                        MAX_NUM_FEATURE_POINTS,
                        SPOOF_targetEnabled,
                        _SPOOF_target_positions__FLOAT_ARRAY,
                        SPOOF_vertex_positions,
                        SPOOF_vertex_normals,
                        SPOOF_triangle_indices,
                        SPOOF_feature_point_positions);
            }
        }




        { // draw scene
            static int tabs = 0;
            if (globals.key_pressed[COW_KEY_TAB]) ++tabs;
            gui_checkbox("DRAGON", &DRAGON, 'd');
            if (!DRAGON) {
                {
                    if (tabs % 3 == 0) {
                        sim.draw(P, V, M4_Identity(), &currentState);
                    } else if (tabs % 3 == 1) { // draw
                        sim.draw(P * V, &currentState);
                    } else { // check stuff being sent to C#
                        eso_begin(PV, SOUP_OUTLINED_TRIANGLES);
                        eso_color(monokai.black);
                        for_(triangle_i, cpp_getNumTriangles()) {
                            for_(d, 3) {
                                eso_vertex(
                                        SPOOF_vertex_positions[3 * SPOOF_triangle_indices[3 * triangle_i + d] + 0],
                                        SPOOF_vertex_positions[3 * SPOOF_triangle_indices[3 * triangle_i + d] + 1],
                                        SPOOF_vertex_positions[3 * SPOOF_triangle_indices[3 * triangle_i + d] + 2]
                                        );
                            }
                        }
                        eso_end();
                    }
                }
            } else { // skinning
                do_once {
                    mat4 RS = M4_RotationAboutXAxis(PI / 2) * M4_Scaling(0.05);
                    dragonHead._applyTransform(RS);
                    dragonBody._applyTransform(M4_Translation(0.0, -0.67, 0.0) * RS);
                };

                const int NUM_BONES = MESH_NUMBER_OF_NODE_LAYERS - 1;

                vec3 boneOrigins          [NUM_BONES + 1];
                vec3 boneOriginsRest      [NUM_BONES + 1];
                vec3 boneXAxisFeaturePoint[NUM_BONES + 1];
                {
                    for_(j, _COUNT_OF(boneOrigins)) {
                        boneOrigins              [j] = get(currentState.x, 9 + j * 10);
                        boneOriginsRest          [j] = get(sim.x_rest    , 9 + j * 10);
                        boneXAxisFeaturePoint    [j] = get(currentState.x, 0 + j * 10);
                    }
                }

                vec3 boneNegativeYAxis[NUM_BONES];
                vec3 bonePositiveXAxis[NUM_BONES];
                {
                    for_(j, _COUNT_OF(boneNegativeYAxis)) {
                        boneNegativeYAxis[j] = normalized(boneOrigins[j + 1] - boneOrigins[j]);
                        bonePositiveXAxis[j] = normalized(boneXAxisFeaturePoint[j] - boneOrigins[j]);
                    }
                }

                { // dragonHead
                    vec3 y = -boneNegativeYAxis[NUM_BONES - 1];
                    vec3 up = { 0.0, 1.0, 0.0 }; 
                    vec3 x = normalized(cross(y, up));
                    vec3 z = cross(x, y);
                    vec3 o = boneOrigins[NUM_BONES];
                    dragonHead.draw(P, V, M4_xyzo(x, y, z, o));
                }

                { // dragonBody
                    do_once {
                        dragonBody.num_bones = NUM_BONES;
                        dragonBody.bones = (mat4 *) malloc(NUM_BONES * sizeof(mat4));
                        dragonBody.bone_indices = (int4 *) malloc(dragonBody.num_vertices * sizeof(int4));
                        dragonBody.bone_weights = (vec4 *) malloc(dragonBody.num_vertices * sizeof(vec4));

                        // assign weights FORNOW hacky nonsense
                        for_(vertex_i, dragonBody.num_vertices) {
                            auto f = [&](int i) {
                                real c = AVG(boneOriginsRest[i].y, boneOriginsRest[i + 1].y);
                                real D = ABS(dragonBody.vertex_positions[vertex_i].y - c);
                                return MAX(0.0, (1.0 / D) - 10.0);
                            };

                            real t = INVERSE_LERP(dragonBody.vertex_positions[vertex_i].y, 0.0, -ROBOT_LENGTH);
                            real b = t * dragonBody.num_bones;

                            int j = MIN(MAX(int(b + 0.5), 0), dragonBody.num_bones - 1);
                            int i = MAX(0, j - 1);
                            int k = MIN(dragonBody.num_bones - 1, j + 1);

                            dragonBody.bone_indices[vertex_i] = { i, j, k };
                            dragonBody.bone_weights[vertex_i] = { f(i), f(j), f(k) };
                            dragonBody.bone_weights[vertex_i] /= sum(dragonBody.bone_weights[vertex_i]);
                        }
                    };
                    for_(bone_i, dragonBody.num_bones) {
                        vec3 y = -boneNegativeYAxis[bone_i];
                        vec3 x = bonePositiveXAxis[bone_i];
                        vec3 z = cross(x, y);
                        mat4 invBind = M4_Translation(-boneOriginsRest[bone_i]);
                        mat4 Bone = M4_xyzo(x, y, z, boneOrigins[bone_i]);
                        dragonBody.bones[bone_i] = Bone * invBind;
                        library.soups.axes.draw(PV * Bone * M4_Scaling(0.08));
                    }
                    dragonBody.draw(P, V, globals.Identity);
                }
            }

            { // ceiling
                real r = 0.3;
                eso_begin(PV, (tabs % 3 == 0) ? SOUP_QUADS : SOUP_OUTLINED_QUADS);
                if (tabs % 3 == 0) {
                    eso_color(1.0, 1.0, 1.0, 0.5);
                } else {
                    eso_color(0.0, 0.0, 0.0, 0.5);
                }
                eso_vertex( r, 0.0,  r);
                eso_vertex( r, 0.0, -r);
                eso_vertex(-r, 0.0, -r);
                eso_vertex(-r, 0.0,  r);
                eso_end();
            }
        }


        { // widget
            bool mouseClickConsumed = false;
            bool mouseHotConsumed = false;
            { // SPOOF_targetPositions
                for_(featurePointIndex, MAX_NUM_FEATURE_POINTS) {
                    if (!SPOOF_targetEnabled[featurePointIndex]) continue;
                    vec3 color = color_kelly(featurePointIndex);
                    vec3 SPOOF_feature_point_position = { SPOOF_feature_point_positions[3 * featurePointIndex + 0], SPOOF_feature_point_positions[3 * featurePointIndex + 1], SPOOF_feature_point_positions[3 * featurePointIndex + 2] };

                    WidgetResult widgetResult = widget(P, V, featurePointIndex, &SPOOF_targetPositions[featurePointIndex], SPOOF_feature_point_position, color);

                    mouseClickConsumed |= widgetResult.mouseClickConsumed; // FORNOW
                    mouseHotConsumed |= widgetResult.mouseHotConsumed; // FORNOW
                    if (widgetResult.pleaseDisableHandle) {
                        SPOOF_targetEnabled[featurePointIndex] = FALSE;
                    }
                    if (widgetResult.recastFeaturePoint) {
                        castRay(true, featurePointIndex);
                    }
                }
            }


            if (!mouseClickConsumed && !mouseHotConsumed) { // SPOOF_intersection_position
                bool pleaseSetFeaturePoint = globals.mouse_left_pressed;
                int featurePointIndex; {
                    for (featurePointIndex = 0; SPOOF_targetEnabled[featurePointIndex]; ++featurePointIndex) {}
                    ASSERT(featurePointIndex < MAX_NUM_FEATURE_POINTS);
                }
                CastRayResult castRayResult = castRay(pleaseSetFeaturePoint, featurePointIndex);
                if (!globals.mouse_left_held && castRayResult.intersects) draw_sphere(castRayResult.intersection_position, color_kelly(featurePointIndex));
                if (castRayResult.intersects && pleaseSetFeaturePoint) {
                    SPOOF_targetEnabled[featurePointIndex] = TRUE;
                    SPOOF_targetPositions[featurePointIndex] = castRayResult.intersection_position;
                }
            }
        }


        { // keyboard input
            if (globals.key_held['a']) SPOOF_targetPositions[0] = transformPoint(M4_RotationAboutYAxis(RAD(1)), SPOOF_targetPositions[0]);
        }

        if (1) { // manual sliders
            for_(j, sim.num_cables) {
                char buffer[] = "u_X";
                buffer[2] = char('0' + j);
                gui_slider(buffer, &currentState.u[j], -(ROBOT_LENGTH / 3), (ROBOT_LENGTH / 3));
            }
        }

        { // FORNOW automated testing
            if (AUTOMATED_SPEED_TEST__QUITS_AFTER_A_COUPLE_SECONDS) {
                static real testTime = 0.0;
                SPOOF_targetPositions[0] += V3(0.003);
                if (testTime > 2.0) {
                    exit(1);
                }
                testTime += .0167;
            }
        }
    }
}

#undef LEN_U
#undef LEN_X

////////////////////////////////////////////////////////////////////////////////

void main() {
    omp_set_num_threads(6);
    APPS {
        // APP(eg_fbo);
        APP(kaa);
        // APP(jones);
    }
}

