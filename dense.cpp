// TODO: plus, minus, times
// TODO: norm
// TODO: squaredNorm
// TODO: dot
template<typename T> struct FixedSizeSelfDestructingArray {
    T *data;
    int N;

    T &operator [](int index) {
        ASSERT(index >= 0);
        ASSERT(index < N);
        return data[index];
    }

    const T &operator [](int index) const {
        ASSERT(index >= 0);
        ASSERT(index < N);
        return data[index];
    }

    FixedSizeSelfDestructingArray() {
        this->N = 0;
        this->data = NULL;
    }

    FixedSizeSelfDestructingArray(int N) {
        this->N = N;
        this->data = (T *) calloc(N, sizeof(T));
    }

    FixedSizeSelfDestructingArray(const FixedSizeSelfDestructingArray &other) {
        this->N = other.N;
        this->data = (T *) calloc(this->N, sizeof(T));
        memcpy(this->data, other.data, this->N * sizeof(T));
    }

    FixedSizeSelfDestructingArray &operator = (const FixedSizeSelfDestructingArray &other) {
        ASSERT(this != &other);
        if (this->N != other.N) {
            free(this->data);
            this->data = (T *) malloc(other.N * sizeof(T));
        }
        memcpy(this->data, other.data, other.N * sizeof(T));
        this->N = other.N;
        return *this;
    }

    ~FixedSizeSelfDestructingArray() {
        if (N != 0) { free(data); }
    }

    void setZero() {
        memset(data, 0, N * sizeof(T));
    }

};

template <typename T> FixedSizeSelfDestructingArray<T> operator + (const FixedSizeSelfDestructingArray<T> &A, const FixedSizeSelfDestructingArray<T> &B) {
    ASSERT(A.N == B.N);
    FixedSizeSelfDestructingArray<T> result = FixedSizeSelfDestructingArray<T>(A.N);
    for_(k, A.N) {
        result.data[k] = A.data[k] + B.data[k];
    }
    return result;
}
typedef FixedSizeSelfDestructingArray<real> SDVector;

real dot(const SDVector &a, const SDVector &b) {
    ASSERT(a.N == b.N);
    real result = 0;
    for_(i, a.N) {
        result += a[i] * b[i];
    }
    return result;
}
real squaredNorm(const SDVector &a) {
    return dot(a, a);
}
real norm(const SDVector &a) {
    return sqrt(squaredNorm(a));
}


struct SDMatrix {
    int R;
    int C;
    real *data;
    SDMatrix() {}
    SDMatrix(int R, int C) {
        this->R = R;
        this->C = C;
        this->data = (real *) calloc(R * C, sizeof(real));
    }
    real &operator() (int r, int c) {
        ASSERT(data);
        ASSERT(r >= 0);
        ASSERT(c >= 0);
        ASSERT(r < R);
        ASSERT(c < C);
        return data[r * C + c];
    }
    const real &operator() (int r, int c) const {
        ASSERT(data);
        ASSERT(r >= 0);
        ASSERT(c >= 0);
        ASSERT(r < R);
        ASSERT(c < C);
        return data[r * C + c];
    }
};

SDVector operator * (const SDVector &a, const SDMatrix &B) {
    ASSERT(a.N == B.R);
    SDVector result(a.N);
    for_(c, B.C) {
        for_(r, B.R) {
            result[c] += a[r] * B(r, c);
        }
    }
    return result;
}





