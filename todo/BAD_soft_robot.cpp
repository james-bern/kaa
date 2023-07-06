#include <vector>
#include "dense.cpp"
#include "hessian.cpp"
#include <omp.h>


// TODO: switch dFdu to SparseMatrix


// FORNOW
// FORNOW
// FORNOW
// FORNOW
// FORNOW

////////////////////////////////////////////////////////////////////////////////

typedef double REAL;

typedef Eigen::Matrix<REAL, 3,  3>  MATRIX3;
typedef Eigen::Matrix<REAL, 9,  9>  MATRIX9;
typedef Eigen::Matrix<REAL, 3,  12> MATRIX3x12;
typedef Eigen::Matrix<REAL, 9,  12> MATRIX9x12;
typedef Eigen::Matrix<REAL, 12, 12> MATRIX12;
typedef Eigen::Matrix<REAL, 2,  1>  VECTOR2;
typedef Eigen::Matrix<REAL, 3,  1>  VECTOR3;
typedef Eigen::Matrix<REAL, 9,  1>  VECTOR9;
typedef Eigen::Matrix<REAL, 12, 1>  VECTOR12;

typedef Eigen::Matrix<int, 2, 1> VECTOR2I;
typedef Eigen::Matrix<int, 3, 1> VECTOR3I;
typedef Eigen::Matrix<int, 4, 1> VECTOR4I;

typedef Eigen::Matrix<REAL, Eigen::Dynamic, 1> VECTOR;
typedef Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> MATRIX;
typedef Eigen::SparseMatrix<REAL> SPARSE_MATRIX;

////////////////////////////////////////////////////////////////////////////////

Hessian global_U_xx;
SDMatrix global_dFdu;

real residualConvergenceThreshold = 1e-5;

real pinSpringConstant = 1.0e7;
real tetMassDensity = 200.0;
real cableSpringConstant = 1.0e9;
real gravitationalConstant = 10.0;

real _tetYoungsModulus = 1000000;
real _tetPoissonsRatio = .47;
real _lambda = (_tetYoungsModulus * _tetPoissonsRatio) / ((1 + _tetPoissonsRatio) * (1 - 2 * _tetPoissonsRatio));
real _mu = _tetYoungsModulus / (2 * (1 + _tetPoissonsRatio));

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
    SDVector nodeMasses;
    SDVector tetRestVolumes;
    SDVector cableReferenceLengths;

    // aesthetics
    int num_triangles;
    int3 *triangle_indices;

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

        nodeMasses = SDVector(num_nodes);
        tetRestVolumes = SDVector(num_tets);
        {
            for_(tet_i, num_tets) {
                Tet tet = tets[tet_i];
                #define e(j) (get(x_rest, tet[j]) - get(x_rest, tet[0]))
                tetRestVolumes[tet_i] = ABS(dot(cross(e(1), e(2)), e(3))) / 6;
                #undef e
                for_(node_ii, SOFT_ROBOT_DIM + 1) {
                    int node_i = tet[node_ii];
                    nodeMasses[node_i] += tetMassDensity * tetRestVolumes[tet_i] / (SOFT_ROBOT_DIM + 1);
                }
            }
        }

        cableReferenceLengths = getCableLengths(x_rest);

        { // num_triangles triangle_indices
            StretchyBuffer<int3> tris = {}; {
                struct { int3 key; int value; } *timesFaceFound = NULL; {
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
                            int3 key = { MIN(MIN(i, j), k), MID(i, j, k), MAX(MAX(i, j), k) };
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
            triangle_indices = (int3 *) malloc(num_triangles * sizeof(int3));
            memcpy(triangle_indices, tris.data, num_triangles * sizeof(int3));
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
                        int3 tri = triangle_indices[sbuff_pop_back(&triangleQueue)];
                        sbuff_push_back(&positivelyOrientedEdgeQueue, { tri.i, tri.j });
                        sbuff_push_back(&positivelyOrientedEdgeQueue, { tri.j, tri.k });
                        sbuff_push_back(&positivelyOrientedEdgeQueue, { tri.k, tri.i });
                    }

                    int2 edge = sbuff_pop_back(&positivelyOrientedEdgeQueue);

                    bool match = false;
                    for_(triangle_i, num_triangles) {
                        if (trianglePushedToQueue[triangle_i]) continue;

                        int3 *tri = &triangle_indices[triangle_i];
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

    static VECTOR9 flatten(const MATRIX3& A) {
        VECTOR9 column;

        unsigned int index = 0;
        for (unsigned int j = 0; j < A.cols(); j++)
            for (unsigned int i = 0; i < A.rows(); i++, index++)
                column[index] = A(i,j);

        return column;
    }

    ///////////////////////////////////////////////////////////////////////
    // convert a VECTOR9 to a MATRIX3 in a consistent way
    ///////////////////////////////////////////////////////////////////////
    static MATRIX3 unflatten(const VECTOR9& v) {
        MATRIX3 A;
        unsigned int index = 0;
        for (unsigned int j = 0; j < A.cols(); j++)
            for (unsigned int i = 0; i < A.rows(); i++, index++)
                A(i,j) = v[index];

        return A;
    }

    MATRIX3 crossProduct(const MATRIX3& F, const int i) {
        return (MATRIX(3,3) <<       0, -F(2,i),  F(1,i),
                F(2,i),       0, -F(0,i),
                -F(1,i),  F(0,i),       0).finished();
    }

    void buildTwistAndFlipEigenvectors(const MATRIX3& U, const MATRIX3& V, MATRIX9& Q)
    {
        // create the twist matrices
        MATRIX3 T0, T1, T2; 
        T0 <<  0, 0, 0,
           0, 0, -1,
           0, 1, 0;   // x-twist
        T1 <<  0, 0, 1,
           0, 0, 0,
           -1, 0, 0;   // y-twist
        T2 <<  0, 1, 0,
           -1, 0, 0,
           0, 0, 0;   // z-twist

        const MATRIX3 Q0 = (1.0 / sqrt(2.0)) * (U * T0 * V.transpose());
        const MATRIX3 Q1 = (1.0 / sqrt(2.0)) * (U * T1 * V.transpose());
        const MATRIX3 Q2 = (1.0 / sqrt(2.0)) * (U * T2 * V.transpose());

        // create the flip matrices
        MATRIX3 L0, L1, L2; 
        L0 <<  0, 0, 0,
           0, 0, 1,
           0, 1, 0;   // x-flip
        L1 <<  0, 0, 1,
           0, 0, 0,
           1, 0, 0;   // y-flip
        L2 <<  0, 1, 0,
           1, 0, 0,
           0, 0, 0;   // z-flip

        const MATRIX3 Q3 = (1.0 / sqrt(2.0)) * (U * L0 * V.transpose());
        const MATRIX3 Q4 = (1.0 / sqrt(2.0)) * (U * L1 * V.transpose());
        const MATRIX3 Q5 = (1.0 / sqrt(2.0)) * (U * L2 * V.transpose());

        Q.col(0) = flatten(Q0);
        Q.col(1) = flatten(Q1);
        Q.col(2) = flatten(Q2);
        Q.col(3) = flatten(Q3);
        Q.col(4) = flatten(Q4);
        Q.col(5) = flatten(Q5);
    }

    void buildScalingEigenvectors(const MATRIX3& U, const MATRIX3& Q,
            const MATRIX3& V, MATRIX9& Q9)
    {
        const VECTOR3 q0 = Q.col(0);
        const VECTOR3 q1 = Q.col(1);
        const VECTOR3 q2 = Q.col(2);

        const MATRIX3 Q0 = U * q0.asDiagonal() * V.transpose();
        const MATRIX3 Q1 = U * q1.asDiagonal() * V.transpose();
        const MATRIX3 Q2 = U * q2.asDiagonal() * V.transpose();

        Q9.col(6) = flatten(Q0);
        Q9.col(7) = flatten(Q1);
        Q9.col(8) = flatten(Q2);
    }

    void ZERO_AND_COMPUTE(State *state, real *E, SDVector *E_x, Hessian *E_xx, SDMatrix *dFdu = NULL) {
        if (E) *E = 0.0;
        if (E_x) *E_x = SDVector(LEN_X);
        if (E_xx) E_xx->init_or_clear(); // TODO: Constructor
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

        if (state->enabled__BIT_FIELD & PINS) {
            const SDVector &x_pin = state->x_pin;
            for_(pin_i, num_pins) {
                int pin = pins[pin_i];
                Vec dx = get(x, pin) - get(x_pin, pin);
                if (E) (*E) += pinSpringConstant / 2 * squaredNorm(dx);
                if (E_x) add(*E_x, pin, pinSpringConstant * dx);
                if (E_xx) add(E_xx, pin, pin, pinSpringConstant * IdentityMatrix<SOFT_ROBOT_DIM>());
            }
        }

        if (state->enabled__BIT_FIELD & GRAVITY) {
            for_(i, num_nodes) {
                Vec minusGravitationalAcceleration = { 0.0, gravitationalConstant, 0.0 };
                if (E) { (*E) += nodeMasses[i] * dot(minusGravitationalAcceleration, get(x, i)); }
                if (E_x) { add(*E_x, i, nodeMasses[i] * minusGravitationalAcceleration); }
            }
        }


        if (state->enabled__BIT_FIELD & TETS) { 
            if (1) {
                // TODO: compute hessian a la hobak
                // TODO: compute clamped hessian a la hobak
                // Does it behave better?--if not revert to safe_soft_robot.cpp
                #define NUM_TETS 2048
                ASSERT(num_tets <= NUM_TETS);
                static real blocks_U[NUM_TETS];
                static Vec blocks_U_x[NUM_TETS][SOFT_ROBOT_DIM + 1];
                static Mat blocks_U_xx[NUM_TETS][SOFT_ROBOT_DIM + 1][SOFT_ROBOT_DIM + 1];
                #undef NUM_TETS

                static std::vector<VECTOR3> _x(num_nodes);
                static std::vector<VECTOR3> _X(num_nodes);
                static std::vector<VECTOR4I> _tets(num_tets);
                static std::vector<MATRIX3> _DmInvs(num_tets);
                static std::vector<MATRIX9x12> _pFpxs(num_tets);
                static bool initialized;
                if (!initialized) {
                    initialized = true;
                    for_(node_i, num_nodes) _X[node_i] = VECTOR3(x_rest[3 * node_i + 0], x_rest[3 * node_i + 1], x_rest[3 * node_i + 2]);
                    for_(tetIndex, num_tets) _tets[tetIndex] = VECTOR4I(tets[tetIndex][0], tets[tetIndex][1], tets[tetIndex][2], tets[tetIndex][3]);
                    for_(tetIndex, num_tets) {
                        { // _DmInvs
                            const VECTOR4I& tet = _tets[tetIndex];
                            MATRIX3 Dm;
                            Dm.col(0) = _X[tet[1]] - _X[tet[0]];
                            Dm.col(1) = _X[tet[2]] - _X[tet[0]];
                            Dm.col(2) = _X[tet[3]] - _X[tet[0]];
                            _DmInvs[tetIndex] = Dm.inverse();
                        }
                        { // _pFpxs
                            const MATRIX3 &DmInv = _DmInvs[tetIndex];
                            const REAL m = DmInv(0, 0);
                            const REAL n = DmInv(0, 1);
                            const REAL o = DmInv(0, 2);
                            const REAL p = DmInv(1, 0);
                            const REAL q = DmInv(1, 1);
                            const REAL r = DmInv(1, 2);
                            const REAL s = DmInv(2, 0);
                            const REAL t = DmInv(2, 1);
                            const REAL u = DmInv(2, 2);

                            const REAL t1 = -m - p - s;
                            const REAL t2 = -n - q - t;
                            const REAL t3 = -o - r - u;

                            MATRIX9x12 PFPu = MATRIX9x12::Zero();
                            PFPu(0, 0)  = t1;
                            PFPu(0, 3)  = m;
                            PFPu(0, 6)  = p;
                            PFPu(0, 9)  = s;
                            PFPu(1, 1)  = t1;
                            PFPu(1, 4)  = m;
                            PFPu(1, 7)  = p;
                            PFPu(1, 10) = s;
                            PFPu(2, 2)  = t1;
                            PFPu(2, 5)  = m;
                            PFPu(2, 8)  = p;
                            PFPu(2, 11) = s;
                            PFPu(3, 0)  = t2;
                            PFPu(3, 3)  = n;
                            PFPu(3, 6)  = q;
                            PFPu(3, 9)  = t;
                            PFPu(4, 1)  = t2;
                            PFPu(4, 4)  = n;
                            PFPu(4, 7)  = q;
                            PFPu(4, 10) = t;
                            PFPu(5, 2)  = t2;
                            PFPu(5, 5)  = n;
                            PFPu(5, 8)  = q;
                            PFPu(5, 11) = t;
                            PFPu(6, 0)  = t3;
                            PFPu(6, 3)  = o;
                            PFPu(6, 6)  = r;
                            PFPu(6, 9)  = u;
                            PFPu(7, 1)  = t3;
                            PFPu(7, 4)  = o;
                            PFPu(7, 7)  = r;
                            PFPu(7, 10) = u;
                            PFPu(8, 2)  = t3;
                            PFPu(8, 5)  = o;
                            PFPu(8, 8)  = r;
                            PFPu(8, 11) = u;

                            _pFpxs[tetIndex] = PFPu;
                        }
                    }
                }

                // FORNOW
                for_(node_i, num_nodes) _x[node_i] = VECTOR3(x[3 * node_i + 0], x[3 * node_i + 1], x[3 * node_i + 2]);

                #pragma omp parallel for default(none) shared(blocks_U_xx) schedule(static, 8)
                for_(tetIndex, num_tets) {
                    const VECTOR4I &tet = _tets[tetIndex];
                    real RestVolume = tetRestVolumes[tetIndex]; // FORNOW
                    MATRIX3 F; {
                        MATRIX3 Ds; {
                            Ds.col(0) = _x[tet[1]] - _x[tet[0]];
                            Ds.col(1) = _x[tet[2]] - _x[tet[0]];
                            Ds.col(2) = _x[tet[3]] - _x[tet[0]];
                        }
                        F = Ds * _DmInvs[tetIndex];
                    }
                    REAL J = F.determinant();
                    REAL logJ = log(J);
                    if (E) blocks_U[tetIndex] = RestVolume * (_mu * 0.5 * (F.squaredNorm() - 3.0) - _mu * logJ + _lambda * 0.5 * logJ * logJ);
                    if (E_x || E_xx) {
                        const MATRIX9x12 &pFpx = _pFpxs[tetIndex];
                        MATRIX3 pJpF; {
                            pJpF.col(0) = F.col(1).cross(F.col(2));
                            pJpF.col(1) = F.col(2).cross(F.col(0));
                            pJpF.col(2) = F.col(0).cross(F.col(1));
                        }
                        if (E_x) {
                            MATRIX3 PK1 = _mu * (F - (1.0 / J) * pJpF) + _lambda * logJ * (1.0 / J) * pJpF;
                            VECTOR12 tetForce = RestVolume * (pFpx.transpose() * flatten(PK1));
                            for_(d, 4) blocks_U_x[tetIndex][d]
                                = { tetForce[3 * d + 0], tetForce[3 * d + 1], tetForce[3 * d + 2] };
                        }
                        if (E_xx) {
                            MATRIX3 U, V;
                            VECTOR3 Sigma;
                            {
                                const Eigen::JacobiSVD<MATRIX3,Eigen::NoQRPreconditioner> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
                                U = svd.matrixU();
                                V = svd.matrixV();
                                Sigma = svd.singularValues();

                                MATRIX3 L = MATRIX3::Identity();
                                L(2,2) = (U * V.transpose()).determinant();

                                const REAL detU = U.determinant();
                                const REAL detV = V.determinant();

                                if (detU < 0.0 && detV > 0)
                                    U = U * L;
                                if (detU > 0.0 && detV < 0.0)
                                    V = V * L;

                                Sigma[2] = Sigma[2] * L(2,2);
                            }
                            VECTOR9 eigenvalues;
                            MATRIX9 eigenvectors;


                            // FORNOW
                            J = Sigma[0] * Sigma[1] * Sigma[2];
                            logJ = log(J);

                            const REAL s0s1inv = 1.0 / (Sigma(0) * Sigma(1));
                            const REAL s0s2inv = 1.0 / (Sigma(0) * Sigma(2));
                            const REAL s1s2inv = 1.0 / (Sigma(1) * Sigma(2));

                            // 0-2 are twist
                            // 3-5 are flip
                            // 6-8 are scaling
                            const REAL front = _lambda * logJ - _mu;
                            eigenvalues[0] = front * s1s2inv + _mu;
                            eigenvalues[1] = front * s0s2inv + _mu;
                            eigenvalues[2] = front * s0s1inv + _mu;
                            eigenvalues[3] = -front * s1s2inv + _mu;
                            eigenvalues[4] = -front * s0s2inv + _mu;
                            eigenvalues[5] = -front * s0s1inv + _mu;

                            // populate matrix for scaling eigenvalues
                            MATRIX3 A;
                            const REAL s0s0 = Sigma(0) * Sigma(0);
                            const REAL s1s1 = Sigma(1) * Sigma(1);
                            const REAL s2s2 = Sigma(2) * Sigma(2);
                            const REAL frontDiag = _lambda * (1.0 - logJ) + _mu;
                            A(0,0) = frontDiag / s0s0 + _mu;
                            A(1,1) = frontDiag / s1s1 + _mu;
                            A(2,2) = frontDiag / s2s2 + _mu;

                            A(0,1) = _lambda * s0s1inv;
                            A(0,2) = _lambda * s0s2inv;
                            A(1,2) = _lambda * s1s2inv;
                            A(1,0) = A(0,1);
                            A(2,0) = A(0,2);
                            A(2,1) = A(1,2);

                            // get the scaling eigenvalues
                            const Eigen::SelfAdjointEigenSolver<MATRIX3> Aeigs(A);
                            for (int d = 0; d < 3; d++)
                                eigenvalues[d + 6] = Aeigs.eigenvalues()[d];

                            // Compute the eigenvectors
                            buildTwistAndFlipEigenvectors(U, V, eigenvectors);
                            buildScalingEigenvectors(U, Aeigs.eigenvectors(), V, eigenvectors);

                            // Clamp the eigenvalues
                            for (int i = 0; i < 9; i++)
                                if (eigenvalues(i) < 0.0)
                                    eigenvalues(i) = 0.0;

                            MATRIX9 hessian = eigenvectors * (RestVolume * eigenvalues).asDiagonal() * eigenvectors.transpose();
                            MATRIX12 H = pFpx.transpose() * (hessian * pFpx);

                            for_(r, 4) for_(c, 4) blocks_U_xx[tetIndex][r][c]
                                = {
                                    H(3 * r + 0, 3 * c + 0), H(3 * r + 0, 3 * c + 1), H(3 * r + 0, 3 * c + 2),
                                    H(3 * r + 1, 3 * c + 0), H(3 * r + 1, 3 * c + 1), H(3 * r + 1, 3 * c + 2),
                                    H(3 * r + 2, 3 * c + 0), H(3 * r + 2, 3 * c + 1), H(3 * r + 2, 3 * c + 2),
                                };
                        }
                    }
                }

                { // serial accumulate
                    if (E) {
                        for_(tet_i, num_tets) {
                            *E += blocks_U[tet_i];
                        }
                    }

                    if (E_x) {
                        for_(tet_i, num_tets) {
                            for_(i, SOFT_ROBOT_DIM + 1) {
                                add(*E_x, tets[tet_i][i], blocks_U_x[tet_i][i]);
                            }
                        }
                    }

                    if (E_xx) {
                        for_(tet_i, num_tets) {
                            for_(i, SOFT_ROBOT_DIM + 1) {
                                for_(j, SOFT_ROBOT_DIM + 1) {
                                    add(E_xx, tets[tet_i][i], tets[tet_i][j], blocks_U_xx[tet_i][i][j]);
                                }
                            }
                        }
                    }
                }
            } else {
                #define NUM_TETS 2048
                ASSERT(num_tets <= NUM_TETS);
                static real blocks_U[NUM_TETS];
                static Vec blocks_U_x[NUM_TETS][SOFT_ROBOT_DIM + 1];
                static Mat blocks_U_xx[NUM_TETS][SOFT_ROBOT_DIM + 1][SOFT_ROBOT_DIM + 1];
                #undef NUM_TETS
                // #pragma omp parallel for default(none) shared(blocks_U_xx) schedule(static, 8)
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
                    real RestVolume = tetRestVolumes[tet_i];
                    Mat VdXinvT = RestVolume * transpose(dXinv);
                    Mat I_dP = (_mu - _lambda * logJ) * FinvT;
                    real Ic = squaredNorm(F);

                    if (E) {
                        blocks_U[tet_i] = RestVolume * (_mu * 0.5 * (Ic - 3.0) - _mu * logJ + _lambda * 0.5 * logJ * logJ);
                    }


                    if (E_x) {
                        Mat dEdF = _mu * F - I_dP;
                        Vec segment[SOFT_ROBOT_DIM + 1]; {
                            segment[1] = { dEdF(0, 0) * dXinv(0, 0) + dEdF(0, 1) * dXinv(0, 1) + dEdF(0, 2) * dXinv(0, 2), dEdF(1, 0) * dXinv(0, 0) + dEdF(1, 1) * dXinv(0, 1) + dEdF(1, 2) * dXinv(0, 2), dEdF(2, 0) * dXinv(0, 0) + dEdF(2, 1) * dXinv(0, 1) + dEdF(2, 2) * dXinv(0, 2) };
                            segment[2] = { dEdF(0, 0) * dXinv(1, 0) + dEdF(0, 1) * dXinv(1, 1) + dEdF(0, 2) * dXinv(1, 2), dEdF(1, 0) * dXinv(1, 0) + dEdF(1, 1) * dXinv(1, 1) + dEdF(1, 2) * dXinv(1, 2), dEdF(2, 0) * dXinv(1, 0) + dEdF(2, 1) * dXinv(1, 1) + dEdF(2, 2) * dXinv(1, 2) };
                            segment[3] = { dEdF(0, 0) * dXinv(2, 0) + dEdF(0, 1) * dXinv(2, 1) + dEdF(0, 2) * dXinv(2, 2), dEdF(1, 0) * dXinv(2, 0) + dEdF(1, 1) * dXinv(2, 1) + dEdF(1, 2) * dXinv(2, 2), dEdF(2, 0) * dXinv(2, 0) + dEdF(2, 1) * dXinv(2, 1) + dEdF(2, 2) * dXinv(2, 2) };
                            segment[0] = -segment[1] - segment[2] - segment[3];
                        }
                        for_(i, SOFT_ROBOT_DIM + 1) blocks_U_x[tet_i][i] = RestVolume * segment[i];
                    }

                    if (E_xx) {
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
                    if (E) {
                        for_(tet_i, num_tets) {
                            *E += blocks_U[tet_i];
                        }
                    }

                    if (E_x) {
                        for_(tet_i, num_tets) {
                            for_(i, SOFT_ROBOT_DIM + 1) {
                                add(*E_x, tets[tet_i][i], blocks_U_x[tet_i][i]);
                            }
                        }
                    }

                    if (E_xx) {
                        for_(tet_i, num_tets) {
                            for_(i, SOFT_ROBOT_DIM + 1) {
                                for_(j, SOFT_ROBOT_DIM + 1) {
                                    add(E_xx, tets[tet_i][i], tets[tet_i][j], blocks_U_xx[tet_i][i][j]);
                                }
                            }
                        }
                    }
                }
            }
        }

        // TODO parallelize cables
        if (state->enabled__BIT_FIELD & CABLES) {
            const SDVector &u = state->u;

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
                if (E) {
                    *E += Q;
                }
                for_line_strip_(via_i, via_j, num_vias[cable_i]) {
                    const Via &I = cable[via_i];
                    const Via &J = cable[via_j];
                    Vec xI_minus_xJ = get(x, I) - get(x, J);
                    Vec firstDerivativeOfNorm__xI_minus_xJ = firstDerivativeOfNorm(xI_minus_xJ);
                    Vec Qpp_firstDerivativeOfNorm__xI_minus_xJ = Qpp * firstDerivativeOfNorm__xI_minus_xJ;
                    if (E_x) { // E_x += Q' dldx
                        add_I_minus_J(*E_x, I, J, Qp * firstDerivativeOfNorm__xI_minus_xJ);
                    }
                    if (E_xx) { // E_xx += Q' d2ldx2
                        add_I_minus_J_times_I_minus_J(E_xx, I, J, Qp * secondDerivativeOfNorm(xI_minus_xJ));
                        for_line_strip_(via_m, via_n, num_vias[cable_i]) { // E_xx += Q'' dldx^T dldx (sparse outer product)
                            const Via &M = cable[via_m];
                            const Via &N = cable[via_n];
                            add_I_minus_J_times_M_minus_N(
                                    E_xx,
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
        // E_x(x + searchDirection) ~ E_x + E_xx searchDirection := 0
        // => getNextX { E_xx searchDirection = -E_x }

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
                        // opt_solve_sparse_linear_system(LEN_X, searchDirection.data, E_xx.numTriplets(), E_xx.dataPointer(), minus_U_x.data);
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
                while (1) {

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
            int3 tri = triangle_indices[triangle_i];
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

    void draw(mat4 P, mat4 RestVolume, mat4 M, State *state) {
        const SDVector &x = state->x;
        mat4 PV = P * RestVolume;

        if (state->enabled__BIT_FIELD & TETS) {
            // TODO should x be a FixedSizeSelfDestructingArray<vec3>?????
            // TODO should x be a FixedSizeSelfDestructingArray<vec3>?????
            // TODO should x be a FixedSizeSelfDestructingArray<vec3>?????
            // TODO should x be a FixedSizeSelfDestructingArray<vec3>?????
            // TODO should x be a FixedSizeSelfDestructingArray<vec3>?????
            SDVector vertex_normals = get_vertex_normals(x);

            static vec3 *vertex_colors;
            if (!vertex_colors) {
                vertex_colors = (vec3 *) malloc(num_nodes * sizeof(vec3));
                for_(node_i, num_nodes) {
                    vertex_colors[node_i] = color_plasma(-1.3 * x_rest[3 * node_i + 1]);
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

            mesh.draw(P, RestVolume, M);
            // mesh._dump_for_library("snake.txt", "snake"); exit(1);
        }

        if (state->enabled__BIT_FIELD & CABLES) {
            for_(cable_i, num_cables) {
                Via *cable = cables[cable_i];
                eso_begin(P * RestVolume * M, SOUP_LINE_STRIP);
                // eso_color(color_kelly(cable_i));
                eso_color(monokai.red);
                for_(via_i, num_vias[cable_i]) eso_vertex(get(x, cable[via_i]));
                eso_end();
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
                eso_begin(PV, SOUP_LINE_STRIP);
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
