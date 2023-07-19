#include <vector>
#include "dense.cpp"
#include "hessian.cpp"
#include <omp.h>

// FORNOW
// FORNOW
// FORNOW
// FORNOW
// FORNOW

////////////////////////////////////////////////////////////////////////////////

// TODO: StretchyHashTable

////////////////////////////////////////////////////////////////////////////////

Hessian global_U_xx;
SDMatrix global_dFdu;

real residualConvergenceThreshold = 1e-5;

real tetMassDensity = 200.0;
real tetYoungsModulus = 1000000;
real tetPoissonsRatio = .47;

real pinSpringConstant = 1.0e7;
real cableSpringConstant = 1.0e9;
real gravitationalConstant = 9.81;


////////////////////////////////////////////////////////////////////////////////

static void *_MEMCPY_ALLOC(void *source, int size) {
    void *result = malloc(size);
    memcpy(result, source, size);
    return result;
}

////////////////////////////////////////////////////////////////////////////////
// ZCQ /////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

const real _ZCQ_DEFAULT_eps__ZERO_CORRESPONDS_TO_UNSMOOTHED = .001;

real ZCQ(real c, real x, real eps = _ZCQ_DEFAULT_eps__ZERO_CORRESPONDS_TO_UNSMOOTHED) {
    return
        (x > eps) ? c * pow(x, 2) - c * eps * x + c * pow(eps, 2) / 3 :
        (x >   0) ? c / (3 * eps) * pow(x, 3) :
        0;
}

real ZCQp(real c, real x, real eps = _ZCQ_DEFAULT_eps__ZERO_CORRESPONDS_TO_UNSMOOTHED) {
    return
        (x > eps) ? 2 * c * x - c * eps :
        (x >   0) ? c / eps * pow(x, 2) :
        0;
}

real ZCQpp(real c, real x, real eps = _ZCQ_DEFAULT_eps__ZERO_CORRESPONDS_TO_UNSMOOTHED) {
    return
        (x > eps) ? 2 * c :
        (x >   0) ? 2 * c / eps * x :
        0;
}

////////////////////////////////////////////////////////////////////////////////

typedef SnailVector<SOFT_ROBOT_DIM> Vec;
typedef SnailMatrix<SOFT_ROBOT_DIM> Mat;
typedef int4 Tet;
typedef int3 Tri;
struct _int_real { int index; real weight; };
struct Via {
    _int_real data[4];
    _int_real &operator [](int index) { return data[index]; }
    const _int_real &operator [](int index) const { return data[index]; }
};

////////////////////////////////////////////////////////////////////////////////

// FORNOW
#define LEN_X (SOFT_ROBOT_DIM * num_nodes)
#define LEN_U num_cables
#define sim_LEN_X (SOFT_ROBOT_DIM * sim->num_nodes)
#define sim_LEN_U sim->num_cables

////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// get and add//////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

Vec get(const SDVector &x, int i) {
    Vec result = {};
    for_(d, SOFT_ROBOT_DIM) {
        result[d] = x[SOFT_ROBOT_DIM * i + d];
    }
    return result;
}

Vec get(const SDVector &x, Via I) {
    Vec result = {};
    for_(ii, SOFT_ROBOT_DIM + 1) {
        for_(d, SOFT_ROBOT_DIM) {
            result[d] += I[ii].weight * x[SOFT_ROBOT_DIM * I[ii].index + d];
        }
    }
    return result;
}

void add(SDVector &L, int i, Vec R) {
    for_(d, SOFT_ROBOT_DIM) {
        L[SOFT_ROBOT_DIM * i + d] += R[d];
    }
}

void add(SDVector &L, Via I, Vec R) {
    for_(ii, SOFT_ROBOT_DIM + 1) { if (I[ii].weight == 0.0) continue;
        for_(d, SOFT_ROBOT_DIM) {
            L[SOFT_ROBOT_DIM * I[ii].index + d] +=  I[ii].weight * R[d];
        }
    }
}

// TODO: Consider manual unroll
void add(Hessian *H, int i, int j, Mat M) {
    if (i < j) return; // fornow
    for_(d_i, SOFT_ROBOT_DIM) {
        int i_H = SOFT_ROBOT_DIM * i + d_i;
        for_(d_j, SOFT_ROBOT_DIM) {
            int j_H = SOFT_ROBOT_DIM * j + d_j;
            if (i_H < j_H) continue; // fornow
            H->add(i_H, j_H, M(d_i, d_j));
        }
    }
}

void add(Hessian *H, Via I, Via J, Mat M) {
    for_(ii, SOFT_ROBOT_DIM + 1) {
        if (I[ii].weight == 0.0) continue;
        int i = I[ii].index;
        for_(jj, SOFT_ROBOT_DIM + 1) {
            if (J[jj].weight == 0.0) continue;
            int j = J[jj].index;
            if (i < j) continue; // fornow
            for_(d_i, SOFT_ROBOT_DIM) {
                int i_H = SOFT_ROBOT_DIM * i + d_i;
                for_(d_j, SOFT_ROBOT_DIM) {
                    int j_H = SOFT_ROBOT_DIM * j + d_j;
                    if (i_H < j_H) continue; // fornow
                    H->add(i_H, j_H, I[ii].weight * J[jj].weight * M(d_i, d_j));
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

enum EnabledFlags {
    TETS    = (1 << 0),
    PINS    = (1 << 1),
    CABLES  = (1 << 2),
    GRAVITY = (1 << 3),
};

struct Sim;
struct State {
    int enabled__BIT_FIELD;
    SDVector x, u;
    SDVector x_pin;
    State() {}
    State(Sim *sim);
};


struct SimInput {
    int num_pins;
    int *pins;

    int num_nodes;
    real *x_rest;

    int num_tets;
    Tet *tets;

    int num_cables;
    int *num_vias;
    Via *vias;
    real *cableSpringConstants;
};

struct Sim {
    int num_pins;
    int num_nodes;
    int num_tets;
    int num_cables;
    int num_cable_vias_total;

    int  *pins;
    Tet  *tets;
    int  *num_vias;
    Via **cables;
    Via  *vias;

    SDVector x_rest;

    FixedSizeSelfDestructingArray<Mat> tet_dXinv;
    SDVector nodeMasses__NOTE_tetDensity_not_yet_multiplied_on;
    SDVector tetRestVolumes;
    SDVector cableReferenceLengths;

    // aesthetics
    int num_triangles;
    Tri *triangle_indices;

    int size_of_x() { return num_nodes * sizeof(Vec); }
    int size_of_u() { return num_cables * sizeof(real); }

    Sim() {}

    void allocAndPrecompute(SimInput *simInput) {
        ASSERT(simInput->num_nodes);
        ASSERT(simInput->num_tets);

        num_pins = simInput->num_pins;
        num_nodes = simInput->num_nodes;
        num_tets = simInput->num_tets;
        num_cables = simInput->num_cables;
        num_cable_vias_total = 0; {
            for_(i, num_cables) num_cable_vias_total += simInput->num_vias[i];
        }

        pins      = (int *) _MEMCPY_ALLOC(simInput->pins, num_pins * sizeof(int));
        tets      = (Tet *) _MEMCPY_ALLOC(simInput->tets, num_tets * sizeof(Tet));
        num_vias  = (int *) _MEMCPY_ALLOC(simInput->num_vias, num_cables * sizeof(int));
        vias      = (Via *) _MEMCPY_ALLOC(simInput->vias, num_cable_vias_total * sizeof(Via));

        cables = (Via **) malloc(num_cables * sizeof(Via *)); {
            Via *cable = vias;
            for_(cable_i, num_cables) {
                cables[cable_i] = cable;
                cable += num_vias[cable_i];
            }
        }

        x_rest = SDVector(LEN_X); {
            memcpy(x_rest.data, simInput->x_rest, size_of_x());
        }

        {
            tet_dXinv = getSimplexEdgeVectorMatrices(x_rest);
            for_(tet_i, num_tets) tet_dXinv[tet_i] = inverse(tet_dXinv[tet_i]);
        }

        nodeMasses__NOTE_tetDensity_not_yet_multiplied_on = SDVector(num_nodes);
        tetRestVolumes = SDVector(num_tets);
        {
            for_(tet_i, num_tets) {
                Tet tet = tets[tet_i];
                #define e(j) (get(x_rest, tet[j]) - get(x_rest, tet[0]))
                tetRestVolumes[tet_i] = ABS(dot(cross(e(1), e(2)), e(3))) / 6;
                #undef e
                for_(node_ii, SOFT_ROBOT_DIM + 1) {
                    int node_i = tet[node_ii];
                    nodeMasses__NOTE_tetDensity_not_yet_multiplied_on[node_i] += tetRestVolumes[tet_i] / (SOFT_ROBOT_DIM + 1);
                }
            }
        }

        cableReferenceLengths = getCableLengths(x_rest);

        { // num_triangles triangle_indices
            StretchyBuffer<Tri> tris = {}; {
                struct { Tri key; int value; } *timesFaceFound = NULL; {
                    hmdefault(timesFaceFound, 0);
                    for_(tet_i, num_tets) {
                        Tet tet = tets[tet_i];
                        for_(d, SOFT_ROBOT_DIM + 1) {
                            int ii[3] = { 0, 1, 2 }; {
                                for (int iii = d; iii < 3; ++iii) ++ii[iii];
                            }
                            int i = tet[ii[0]];
                            int j = tet[ii[1]];
                            int k = tet[ii[2]];
                            Tri key = { MIN(MIN(i, j), k), MID(i, j, k), MAX(MAX(i, j), k) };
                            {
                                int tmp = hmget(timesFaceFound, key);
                                hmput(timesFaceFound, key, tmp + 1);
                            }
                        }
                    }
                }
                for_(i, hmlen(timesFaceFound)) {
                    if (timesFaceFound[i].value == 1) {
                        sbuff_push_back(&tris, timesFaceFound[i].key);
                    }
                }
                hmfree(timesFaceFound);
            }
            num_triangles = tris.length;
            triangle_indices = (Tri *) malloc(num_triangles * sizeof(Tri));
            memcpy(triangle_indices, tris.data, num_triangles * sizeof(Tri));
            ASSERT(num_triangles != 0);

            { // orient boundary
                StretchyBuffer<int> triangleQueue = {};
                bool *trianglePushedToQueue = (bool *) calloc(1, num_triangles * sizeof(bool));
                defer { free(trianglePushedToQueue); };

                auto PUSH_TRIANGLE = [&](int triangle_i) {
                    if (!trianglePushedToQueue[triangle_i]) {
                        sbuff_push_back(&triangleQueue, triangle_i);
                        trianglePushedToQueue[triangle_i] = true;
                    }
                };

                PUSH_TRIANGLE(1);

                StretchyBuffer<int2> positivelyOrientedEdgeQueue = {};

                while ((triangleQueue.length != 0) || (positivelyOrientedEdgeQueue.length != 0)) {
                    if (positivelyOrientedEdgeQueue.length == 0) {
                        Tri tri = triangle_indices[sbuff_pop_back(&triangleQueue)];
                        sbuff_push_back(&positivelyOrientedEdgeQueue, { tri.i, tri.j });
                        sbuff_push_back(&positivelyOrientedEdgeQueue, { tri.j, tri.k });
                        sbuff_push_back(&positivelyOrientedEdgeQueue, { tri.k, tri.i });
                    }

                    int2 edge = sbuff_pop_back(&positivelyOrientedEdgeQueue);

                    bool match = false;
                    for_(triangle_i, num_triangles) {
                        if (trianglePushedToQueue[triangle_i]) continue;

                        Tri *tri = &triangle_indices[triangle_i];
                        for_line_loop_(i, j, 3) {
                            if (((*tri)[i] == edge.i) && ((*tri)[j] == edge.j)) {
                                match = true;
                                { // flip orientation
                                    int tmp = (*tri)[0];
                                    (*tri)[0] = (*tri)[1];
                                    (*tri)[1] = tmp;
                                }
                                break;
                            }
                            if (((*tri)[i] == edge.j) && ((*tri)[j] == edge.i)) {
                                match = true;
                                break;
                            }
                        }
                        if (match) {
                            PUSH_TRIANGLE(triangle_i);
                            break;
                        }
                    }
                }

                for_(triangle_i, num_triangles) ASSERT(trianglePushedToQueue[triangle_i]);
            }
        }
    }

    SDVector getCableLengths(const SDVector &x) {
        SDVector result = SDVector(LEN_U);
        for_(cable_i, num_cables) {
            Via *cable = cables[cable_i];
            for_line_strip_(via_i, via_j, num_vias[cable_i]) {
                result[cable_i] += norm(get(x, cable[via_i]) - get(x, cable[via_j]));
            }
        }
        return result;
    }

    SDVector computeCableDeltas(const SDVector &x, const SDVector &u) {
        SDVector Delta = getCableLengths(x);
        for_(i, num_cables) Delta[i] += (u[i] - cableReferenceLengths[i]);
        return Delta;
    }

    FixedSizeSelfDestructingArray<Mat> getSimplexEdgeVectorMatrices(const SDVector &x) {
        FixedSizeSelfDestructingArray<Mat> result(num_tets);
        for_(tet_i, num_tets) {
            result[tet_i].data[0] = x[3 * tets[tet_i][1] + 0] - x[3 * tets[tet_i][0] + 0];
            result[tet_i].data[1] = x[3 * tets[tet_i][2] + 0] - x[3 * tets[tet_i][0] + 0];
            result[tet_i].data[2] = x[3 * tets[tet_i][3] + 0] - x[3 * tets[tet_i][0] + 0];
            result[tet_i].data[3] = x[3 * tets[tet_i][1] + 1] - x[3 * tets[tet_i][0] + 1];
            result[tet_i].data[4] = x[3 * tets[tet_i][2] + 1] - x[3 * tets[tet_i][0] + 1];
            result[tet_i].data[5] = x[3 * tets[tet_i][3] + 1] - x[3 * tets[tet_i][0] + 1];
            result[tet_i].data[6] = x[3 * tets[tet_i][1] + 2] - x[3 * tets[tet_i][0] + 2];
            result[tet_i].data[7] = x[3 * tets[tet_i][2] + 2] - x[3 * tets[tet_i][0] + 2];
            result[tet_i].data[8] = x[3 * tets[tet_i][3] + 2] - x[3 * tets[tet_i][0] + 2];
        }
        return result;
    }


    void ZERO_AND_COMPUTE(State *state, real *U, SDVector *U_x, Hessian *U_xx, SDMatrix *dFdu = NULL) {
        if (U) *U = 0.0;
        if (U_x) *U_x = SDVector(LEN_X);
        if (U_xx) U_xx->init_or_clear();
        if (dFdu) *dFdu = SDMatrix(LEN_X, LEN_U);

        #define add_I_minus_J(G, I, J, __tmp) do { Vec __segment = __tmp; \
            add(G, I,  __segment); \
            add(G, J, -__segment); \
        } while (0)
        #define add_I_minus_J_times_I_minus_J(H, I, J, __tmp) do { Mat __block = __tmp; \
            add(H, I, I,  __block); \
            add(H, I, J, -__block); \
            add(H, J, J,  __block); \
            add(H, J, I, -__block); \
        } while (0)
        #define add_I_minus_J_times_M_minus_N(H, I, J, M, N, __tmp) do { Mat __block = __tmp; \
            add(H, I, M,  __block); \
            add(H, I, N, -__block); \
            add(H, J, N,  __block); \
            add(H, J, M, -__block); \
        } while (0)

        const SDVector &x     = state->x;
        const SDVector &u     = state->u;
        const SDVector &x_pin = state->x_pin;


        if (state->enabled__BIT_FIELD & PINS) {
            for_(pin_i, num_pins) {
                int pin = pins[pin_i];
                Vec dx = get(x, pin) - get(x_pin, pin);
                if (U) (*U) += pinSpringConstant / 2 * squaredNorm(dx);
                if (U_x) add(*U_x, pin, pinSpringConstant * dx);
                if (U_xx) add(U_xx, pin, pin, pinSpringConstant * identityMatrix<SOFT_ROBOT_DIM>());
            }
        }

        if (state->enabled__BIT_FIELD & GRAVITY) {
            for_(i, num_nodes) {
                Vec minusGravitationalAcceleration = { 0.0, gravitationalConstant, 0.0 };
                if (U) { (*U) += tetMassDensity * nodeMasses__NOTE_tetDensity_not_yet_multiplied_on[i] * dot(minusGravitationalAcceleration, get(x, i)); }
                if (U_x) { add(*U_x, i, tetMassDensity * nodeMasses__NOTE_tetDensity_not_yet_multiplied_on[i] * minusGravitationalAcceleration); }
            }
        }


        if (state->enabled__BIT_FIELD & TETS) { 

            real _lambda = (tetYoungsModulus * tetPoissonsRatio) / ((1 + tetPoissonsRatio) * (1 - 2 * tetPoissonsRatio));
            real _mu = tetYoungsModulus / (2 * (1 + tetPoissonsRatio));

            #define NUM_TETS 2048
            ASSERT(num_tets <= NUM_TETS);
            static real blocks_U[NUM_TETS];
            static Vec blocks_U_x[NUM_TETS][SOFT_ROBOT_DIM + 1];
            static Mat blocks_U_xx[NUM_TETS][SOFT_ROBOT_DIM + 1][SOFT_ROBOT_DIM + 1];
            #undef NUM_TETS
            #pragma omp parallel for default(none) shared(blocks_U_xx) schedule(static, 8)
            for_(tet_i, num_tets) {
                Tet tet = tets[tet_i];
                const Mat &dXinv = tet_dXinv[tet_i];
                Mat F = M3(
                        x[3 * tets[tet_i][1] + 0] - x[3 * tets[tet_i][0] + 0],
                        x[3 * tets[tet_i][2] + 0] - x[3 * tets[tet_i][0] + 0],
                        x[3 * tets[tet_i][3] + 0] - x[3 * tets[tet_i][0] + 0],
                        x[3 * tets[tet_i][1] + 1] - x[3 * tets[tet_i][0] + 1],
                        x[3 * tets[tet_i][2] + 1] - x[3 * tets[tet_i][0] + 1],
                        x[3 * tets[tet_i][3] + 1] - x[3 * tets[tet_i][0] + 1],
                        x[3 * tets[tet_i][1] + 2] - x[3 * tets[tet_i][0] + 2],
                        x[3 * tets[tet_i][2] + 2] - x[3 * tets[tet_i][0] + 2],
                        x[3 * tets[tet_i][3] + 2] - x[3 * tets[tet_i][0] + 2]
                        ) * dXinv;
                Mat Finv = inverse(F);
                Mat FinvT = transpose(Finv);
                real J = determinant(F);
                real logJ = log(J);
                real V = tetRestVolumes[tet_i];
                Mat VdXinvT = V * transpose(dXinv);
                Mat I_dP = (_mu - _lambda * logJ) * FinvT;
                real Ic = squaredNorm(F);

                if (U) {
                    blocks_U[tet_i] = V * (_mu * 0.5 * (Ic - 3.0) - _mu * logJ + _lambda * 0.5 * logJ * logJ);
                }


                if (U_x) {
                    Mat dEdF = _mu * F - I_dP;
                    Vec segment[SOFT_ROBOT_DIM + 1]; {
                        segment[1] = { dEdF(0, 0) * dXinv(0, 0) + dEdF(0, 1) * dXinv(0, 1) + dEdF(0, 2) * dXinv(0, 2), dEdF(1, 0) * dXinv(0, 0) + dEdF(1, 1) * dXinv(0, 1) + dEdF(1, 2) * dXinv(0, 2), dEdF(2, 0) * dXinv(0, 0) + dEdF(2, 1) * dXinv(0, 1) + dEdF(2, 2) * dXinv(0, 2) };
                        segment[2] = { dEdF(0, 0) * dXinv(1, 0) + dEdF(0, 1) * dXinv(1, 1) + dEdF(0, 2) * dXinv(1, 2), dEdF(1, 0) * dXinv(1, 0) + dEdF(1, 1) * dXinv(1, 1) + dEdF(1, 2) * dXinv(1, 2), dEdF(2, 0) * dXinv(1, 0) + dEdF(2, 1) * dXinv(1, 1) + dEdF(2, 2) * dXinv(1, 2) };
                        segment[3] = { dEdF(0, 0) * dXinv(2, 0) + dEdF(0, 1) * dXinv(2, 1) + dEdF(0, 2) * dXinv(2, 2), dEdF(1, 0) * dXinv(2, 0) + dEdF(1, 1) * dXinv(2, 1) + dEdF(1, 2) * dXinv(2, 2), dEdF(2, 0) * dXinv(2, 0) + dEdF(2, 1) * dXinv(2, 1) + dEdF(2, 2) * dXinv(2, 2) };
                        segment[0] = -segment[1] - segment[2] - segment[3];
                    }
                    for_(i, SOFT_ROBOT_DIM + 1) blocks_U_x[tet_i][i] = V * segment[i];
                }

                if (U_xx) {
                    Mat dFdXij[4][3] = {}; {
                        for_(i, 4) {
                            for_(j, 3) {
                                if (i > 0) {
                                    dFdXij[i][j](j, i - 1) = 1;
                                } else {
                                    dFdXij[i][j](j, 0) = dFdXij[i][j](j, 1) = dFdXij[i][j](j, 2) = -1;
                                }
                                dFdXij[i][j] = dFdXij[i][j] * dXinv;
                            }
                        }
                    }
                    Mat dF, dP, tmpM, dH;
                    for_(i, 12) {
                        dF = dFdXij[i / 3][i % 3];
                        dP = _mu * dF;
                        dP = dP + I_dP * transpose(dF) * FinvT;
                        tmpM = Finv * dF;
                        dP = dP + _lambda * (tmpM(0, 0) + tmpM(1, 1) + tmpM(2, 2)) * FinvT;
                        dH = dP * VdXinvT;
                        for_(ii, 3) {
                            for_(jj, 3) {
                                blocks_U_xx[tet_i][ii + 1][i / 3](jj, i % 3) = dH(jj, ii);
                            }
                        }
                        blocks_U_xx[tet_i][0][i / 3](0, i % 3) = -dH(0, 2) - dH(0, 1) - dH(0, 0);
                        blocks_U_xx[tet_i][0][i / 3](1, i % 3) = -dH(1, 2) - dH(1, 1) - dH(1, 0);
                        blocks_U_xx[tet_i][0][i / 3](2, i % 3) = -dH(2, 2) - dH(2, 1) - dH(2, 0);
                    }
                }
            }

            { // serial accumulate
                if (U) {
                    for_(tet_i, num_tets) {
                        *U += blocks_U[tet_i];
                    }
                }

                if (U_x) {
                    for_(tet_i, num_tets) {
                        for_(i, SOFT_ROBOT_DIM + 1) {
                            add(*U_x, tets[tet_i][i], blocks_U_x[tet_i][i]);
                        }
                    }
                }

                if (U_xx) {
                    // #pragma omp parallel for default(none) shared(blocks_U_xx) schedule(static, 8)
                    for_(tet_i, num_tets) {
                        for_(i, SOFT_ROBOT_DIM + 1) {
                            for_(j, SOFT_ROBOT_DIM + 1) {
                                add(U_xx, tets[tet_i][i], tets[tet_i][j], blocks_U_xx[tet_i][i][j]);
                            }
                        }
                    }
                }
            }
        }

        // TODO parallelize cables
        if (state->enabled__BIT_FIELD & CABLES) {

            // TODOLATER move inside
            SDVector cableDeltas = computeCableDeltas(x, u);

            // #define NUM_CABLES 9
            // ASSERT(num_cables <= NUM_CABLES);
            // static real blocks_U[NUM_CABLES];
            // static Vec blocks_U_x[NUM_CABLES][SOFT_ROBOT_DIM + 1];
            // static Mat blocks_U_xx[NUM_CABLES][SOFT_ROBOT_DIM + 1][SOFT_ROBOT_DIM + 1];
            // #undef NUM_CABLES
            // #pragma omp parallel for default(none) shared(blocks_U_xx) schedule(static, 8)
            for_(cable_i, num_cables) {
                Via *cable = cables[cable_i];
                real   Q =   ZCQ(cableSpringConstant, cableDeltas[cable_i]);
                real  Qp =  ZCQp(cableSpringConstant, cableDeltas[cable_i]);
                real Qpp = ZCQpp(cableSpringConstant, cableDeltas[cable_i]);
                if (U) {
                    *U += Q;
                }
                for_line_strip_(via_i, via_j, num_vias[cable_i]) {
                    const Via &I = cable[via_i];
                    const Via &J = cable[via_j];
                    Vec xI_minus_xJ = get(x, I) - get(x, J);
                    Vec firstDerivativeOfNorm__xI_minus_xJ = firstDerivativeOfNorm(xI_minus_xJ);
                    Vec Qpp_firstDerivativeOfNorm__xI_minus_xJ = Qpp * firstDerivativeOfNorm__xI_minus_xJ;
                    if (U_x) { // U_x += Q' dldx
                        add_I_minus_J(*U_x, I, J, Qp * firstDerivativeOfNorm__xI_minus_xJ);
                    }
                    if (U_xx) { // U_xx += Q' d2ldx2
                        add_I_minus_J_times_I_minus_J(U_xx, I, J, Qp * secondDerivativeOfNorm(xI_minus_xJ));
                        for_line_strip_(via_m, via_n, num_vias[cable_i]) { // U_xx += Q'' dldx^T dldx (sparse outer product)
                            const Via &M = cable[via_m];
                            const Via &N = cable[via_n];
                            add_I_minus_J_times_M_minus_N(
                                    U_xx,
                                    I, J,
                                    M, N,
                                    outer(Qpp_firstDerivativeOfNorm__xI_minus_xJ, firstDerivativeOfNorm(get(x, M) - get(x, N))));
                        }
                    }
                    if (dFdu) {
                        // SDVector column(size_of_x()); // FORNOW
                        // add_I_minus_J(column, I, J, -Qpp * firstDerivativeOfNorm__xI_minus_xJ);
                        // for_(r, LEN_X) (*dFdu)(r, cable_i) += column[r];
                        Vec segment = -Qpp_firstDerivativeOfNorm__xI_minus_xJ;
                        for_(kk, SOFT_ROBOT_DIM + 1) {
                            for_(d, SOFT_ROBOT_DIM) {
                                (*dFdu)(SOFT_ROBOT_DIM * I[kk].index + d, cable_i) += I[kk].weight * segment[d];
                                (*dFdu)(SOFT_ROBOT_DIM * J[kk].index + d, cable_i) -= J[kk].weight * segment[d];
                            }
                        }
                    }
                }
            }
        }
    }

    State getNext(State *stateWithCurrXandNextU) {
        // U_x(x + searchDirection) ~ U_x + U_xx searchDirection := 0
        // => getNextX { U_xx searchDirection = -U_x }

        State next = *stateWithCurrXandNextU;

        int iterationOfNewtonWithLineSearch = 0;
        while (1) {

            SDVector searchDirection; {
                SDVector minus_U_x;
                ZERO_AND_COMPUTE(&next, NULL, &minus_U_x, &global_U_xx, &global_dFdu);
                for_(i, LEN_X) minus_U_x[i] *= -1;

                real residual = norm(minus_U_x);
                if (residual < residualConvergenceThreshold) {
                    break;
                }
                if (++iterationOfNewtonWithLineSearch > 50) {
                    printf("phyiscs getNextX failed with residual = %lf\n", residual);
                    break;
                }

                { // get searchDirection (and regularize if needed)
                    int iterationofDynamicRegularization = 0;
                    do {
                        if (iterationofDynamicRegularization == 1) { printf("not a descent direction; regularizing\n"); }
                        searchDirection = global_U_xx.solve(minus_U_x);
                        { // regularize Hessian
                            for_(i, LEN_X) global_U_xx.add(i, i, pow(10, -4 + int(iterationofDynamicRegularization)));
                            ++iterationofDynamicRegularization;
                        }
                        if (iterationofDynamicRegularization == 20) {
                            printf("dynamic regularization failed, just going with whatever we have now\n");
                            break;
                        }
                    } while (dot(searchDirection, minus_U_x) < 0);
                }
            }

            { // line search
                real O0;
                ZERO_AND_COMPUTE(stateWithCurrXandNextU, &O0, NULL, NULL);

                real stepSize = 1.0;
                int iterationOfLineSearch = 0;
                SDVector next_x_0 = next.x;
                while (true) {

                    for_(i, LEN_X) next.x[i] = next_x_0[i] + stepSize * searchDirection[i];

                    real O;
                    ZERO_AND_COMPUTE(&next, &O, NULL, NULL);

                    if (O < O0 + TINY_VAL) {
                        break;
                    }

                    if (++iterationOfLineSearch > 30) {
                        printf("line search failed; taking nominal step.\n");
                        break;
                    }

                    stepSize /= 2;
                }
            }
        }

        return next;
    }

    SDVector get_vertex_normals(const SDVector &x) {
        SDVector result(LEN_X);

        for_(triangle_i, num_triangles) {
            Tri tri = triangle_indices[triangle_i];
            vec3 n = cross(
                    get(x, tri[1]) - get(x, tri[0]),
                    get(x, tri[2]) - get(x, tri[0]));
            for_(d, 3) add(result, tri[d], n);
        }

        for_(node_i, num_nodes) {
            real tmp = norm(get(result, node_i));
            if (!IS_ZERO(tmp)) {
                for_(d, 3) result[3 * node_i + d ] /= tmp;   
            }
        }

        return result;
    }

    void draw(mat4 P, mat4 V, mat4 M, State *state) {
        const SDVector &x = state->x;
        mat4 PV = P * V;

        if (1) {
            if (state->enabled__BIT_FIELD & TETS) {
                SDVector vertex_normals = get_vertex_normals(x);

                static vec3 *vertex_colors;
                if (!vertex_colors) {
                    vertex_colors = (vec3 *) malloc(num_nodes * sizeof(vec3));
                    for_(node_i, num_nodes) {
                        vertex_colors[node_i] = color_rainbow_swirl(-1.3 * x_rest[3 * node_i + 1]);
                    }
                }

                IndexedTriangleMesh3D mesh = {}; {
                    mesh.num_vertices = num_nodes;
                    mesh.num_triangles = num_triangles;
                    mesh.vertex_positions = (vec3 *) x.data;
                    mesh.vertex_normals = (vec3 *) vertex_normals.data;
                    mesh.vertex_colors = vertex_colors;
                    mesh.triangle_indices = triangle_indices;
                }

                mesh.draw(P, V, M);
                // mesh._dump_for_library("snake.txt", "snake"); exit(1);
            }
        }

        if (1) {
            if (state->enabled__BIT_FIELD & CABLES) {

                // static vec3 *vertex_colors;
                // if (!vertex_colors) {
                //     vertex_colors = (vec3 *) malloc(library.meshes.bunny.num_vertices * sizeof(vec3));
                //     for_(node_i, library.meshes.bunny.num_vertices) {
                //         vertex_colors[node_i] = color_rainbow_swirl(-1.3 * library.meshes.bunny.vertex_positions[node_i].y);
                //     }
                //     library.meshes.bunny.vertex_colors = vertex_colors;
                // }
                // for_(i, 256) library.meshes.bunny.draw(P, V, M4_Scaling(0.25));

                // library.meshes.bunny.draw(P, V, M4_Scaling(0.2), monokai.yellow);

                for_(cable_i, num_cables) {
                    Via *cable = cables[cable_i];
                    vec3 color = color_kelly(cable_i);
                    real r = 0.004;
                    for_(via_i, num_vias[cable_i]) library.meshes.sphere.draw(P, V, M4_Translation(get(x, cable[via_i])) * M4_Scaling(r), color);
                    for_line_strip_(via_i, via_j, num_vias[cable_i]) {
                        vec3 x_i = get(x, cable[via_i]);
                        vec3 x_j = get(x, cable[via_j]);
                        vec3 e = x_j - x_i;
                        real mag_e = norm(e);
                        vec3 E = { 0.0, 1.0, 0.0 };
                        mat4 R = M4_RotationFrom(E, e);
                        library.meshes.cylinder.draw(P, V, M4_Translation(0.5 * (x_i + x_j)) * R * M4_Scaling(r, mag_e, r) * M4_Translation(0.0, -0.5, 0.0), color);
                    }
                }
            }
        }
    }

    void draw(mat4 PV, State *state) {
        const SDVector &x = state->x;

        if (state->enabled__BIT_FIELD & TETS) {
            eso_begin(PV, SOUP_OUTLINED_TRIANGLES);
            for_(tet_i, num_tets) {
                // eso_color(color_kelly(tet_i));
                eso_color(monokai.black);
                int ii[] = { 0, 1, 2, 0, 1, 3, 0, 2, 3, 1, 2, 3 };
                for_(iii, 12) {
                    eso_vertex(get(x, tets[tet_i][iii[ii]]));
                }
            }
            eso_end();
        }

        if (state->enabled__BIT_FIELD & CABLES) {
            for_(cable_i, num_cables) {
                Via *cable = cables[cable_i];
                eso_begin(PV, SOUP_LINE_STRIP, 10.0);
                eso_color(color_kelly(cable_i));
                for_(via_i, num_vias[cable_i]) eso_vertex(get(x, cable[via_i]));
                eso_end();
            }
        }

        if (state->enabled__BIT_FIELD & PINS) {
            eso_begin(PV, SOUP_POINTS);
            for_(pin_i, num_pins) {
                eso_color(color_kelly(pin_i));
                eso_vertex(get(x, pins[pin_i]));
            }
            eso_end();
        }
    }
    #undef add_I_minus_J
    #undef add_I_minus_J_times_I_minus_J
    #undef add_I_minus_J_times_M_minus_N
};

State::State(Sim *sim) {
    this->enabled__BIT_FIELD = int(-1);
    this->x = sim->x_rest;
    this->x_pin = sim->x_rest;
    this->u = SDVector(sim_LEN_U);
}

#undef LEN_X
#undef LEN_U
