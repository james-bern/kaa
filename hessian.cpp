#if 1

struct Hessian {

        Eigen::SparseMatrix<real> eigenSparseMatrix;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<real>, Eigen::Upper> solver;

        void init_or_clear(int = 0, int = 0) {
            FORNOW_NUM_SLOTS = 100000;
            if (array == NULL) {
                array = (OptEntry *) calloc(FORNOW_NUM_SLOTS, sizeof(OptEntry));
                indices = (int *) malloc(FORNOW_NUM_SLOTS * sizeof(int));
            }
            numAdded = 0;
            for_(k, uniqueNumAdded) array[k].val = 0;
            FORNOW_solvedSinceLastClearOrAdd = false;
        }

        void add(int i, int j, real value) {
            ASSERT(numAdded < FORNOW_NUM_SLOTS);
            ASSERT(i >= j);

            std::pair<int, int> ij = std::make_pair(i, j);

            if (numAdded < maxNumAdded) {

                array[indices[numAdded]].val += value;

            } else {

                bool found = false;
                for_(k, uniqueNumAdded) {
                    if (array[k].i == i && array[k].j == j) {
                        found = true;
                        indices[numAdded] = k;
                        array[k].val += value;
                        break;
                    }
                }

                if (!found) {
                    array[uniqueNumAdded].i = i;
                    array[uniqueNumAdded].j = j;
                    array[uniqueNumAdded].val = value;
                    indices[numAdded] = uniqueNumAdded;
                    ++uniqueNumAdded;
                }

            }

            ++numAdded;
            maxNumAdded = MAX(maxNumAdded, numAdded);
            FORNOW_solvedSinceLastClearOrAdd = false;
        }

        SDVector solve(SDVector b) {
            int N = b.N;

            if (!solvedAtLeastOnce) eigenSparseMatrix = Eigen::SparseMatrix<real>(N, N);
            if (!FORNOW_solvedSinceLastClearOrAdd) {
                eigenSparseMatrix.setFromTriplets(array, array + uniqueNumAdded);
                if (!solvedAtLeastOnce) solver.analyzePattern(eigenSparseMatrix);
                solver.factorize(eigenSparseMatrix);
            }

            Eigen::Map<EigenVectorXr> eigenDenseVector(b.data, N);
            // EigenVectorXr eigenDenseVector(N);
            // for_(i, N) eigenDenseVector[i] = b[i];

            SDVector result(N); {
                EigenVectorXr tmp = solver.solve(eigenDenseVector);
                memcpy(result.data, tmp.data(), N * sizeof(real));
            }

            FORNOW_solvedSinceLastClearOrAdd = true;
            solvedAtLeastOnce = true;

            return result;
        }

        bool solvedAtLeastOnce;
        OptEntry *array;
        int *indices;
        int numTriplets() { return uniqueNumAdded; }
        OptEntry *dataPointer() { return array; }
        int maxNumAdded;
        int numAdded;
        int uniqueNumAdded;
        int FORNOW_NUM_SLOTS;
        bool FORNOW_solvedSinceLastClearOrAdd;

};

#else
struct Hessian {
    public:
        void init_or_clear(int = 0, int = 0) {
            FORNOW_NUM_SLOTS = 100000;
            if (array == NULL) {
                array = (OptEntry *) calloc(FORNOW_NUM_SLOTS, sizeof(OptEntry));
                indices = (int *) malloc(FORNOW_NUM_SLOTS * sizeof(int));
            }
            numAdded = 0;
            for_(k, uniqueNumAdded) array[k].val = 0;
            FORNOW_solvedSinceLastClearOrAdd = false;
        }

        void add(int i, int j, real value) {
            ASSERT(numAdded < FORNOW_NUM_SLOTS);
            ASSERT(i >= j);

            if (numAdded < maxNumAdded) {

                array[indices[numAdded]].val += value;

            } else {

                bool found = false;
                for_(k, uniqueNumAdded) {
                    if (array[k].i == i && array[k].j == j) {
                        found = true;
                        indices[numAdded] = k;
                        array[k].val += value;
                        break;
                    }
                }

                if (!found) {
                    array[uniqueNumAdded].i = i;
                    array[uniqueNumAdded].j = j;
                    array[uniqueNumAdded].val = value;
                    indices[numAdded] = uniqueNumAdded;
                    ++uniqueNumAdded;
                }

            }

            ++numAdded;
            maxNumAdded = MAX(maxNumAdded, numAdded);
            FORNOW_solvedSinceLastClearOrAdd = false;
        }

        SDVector solve(SDVector b) {
            int N = b.N;

            if (!solvedAtLeastOnce) eigenSparseMatrix = Eigen::SparseMatrix<real>(N, N);
            if (!FORNOW_solvedSinceLastClearOrAdd) {
                eigenSparseMatrix.setFromTriplets(array, array + uniqueNumAdded);
                if (!solvedAtLeastOnce) solver.analyzePattern(eigenSparseMatrix);
                solver.factorize(eigenSparseMatrix);
            }

            Eigen::Map<EigenVectorXr> eigenDenseVector(b.data, N);
            // EigenVectorXr eigenDenseVector(N);
            // for_(i, N) eigenDenseVector[i] = b[i];

            SDVector result(N); {
                EigenVectorXr tmp = solver.solve(eigenDenseVector);
                memcpy(result.data, tmp.data(), N * sizeof(real));
            }

            FORNOW_solvedSinceLastClearOrAdd = true;
            solvedAtLeastOnce = true;

            return result;
        }

    private:
        bool solvedAtLeastOnce;
        OptEntry *array;
        int *indices;
        int numTriplets() { return uniqueNumAdded; }
        OptEntry *dataPointer() { return array; }
        int maxNumAdded;
        int numAdded;
        int uniqueNumAdded;
        int FORNOW_NUM_SLOTS;
        bool FORNOW_solvedSinceLastClearOrAdd;
        Eigen::SparseMatrix<real> eigenSparseMatrix;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<real>, Eigen::Upper> solver;
};
#endif
