

void check_derivatives(SolveInput *solveInput, int enabledBitField = 0, real fd_stepsize = 1e-5) {
    if (!enabledBitField) enabledBitField = solveInput->enabledBitField;

    #define ABSOLUTE_ERROR_THRESHOLD .0001
    #define RELATIVE_ERROR_THRESHOLD .001
    #define TRIGGER_INVALID(a) (isnan(a) || isinf(a))
    #define TRIGGER_ERROR(error) (error > ABSOLUTE_ERROR_THRESHOLD && (2. * error / (fabs(a) + fabs(b))) > RELATIVE_ERROR_THRESHOLD)

    real *U_x    = (real *) calloc(LEN_X, sizeof(real));
    real *U_x_fd = (real *) calloc(LEN_X, sizeof(real));
    SparseMatrixEntry *U_xx_    = 0; arrsetcap(U_xx_, LEN_X * LEN_X); 
    SparseMatrixEntry *U_xx_fd_ = 0; arrsetcap(U_xx_fd_, LEN_X * LEN_X); 
    Matrix U_xx    = { LEN_X, LEN_X, (real *) calloc(LEN_X * LEN_X, sizeof(real)) };
    Matrix U_xx_fd = { LEN_X, LEN_X, (real *) calloc(LEN_X * LEN_X, sizeof(real)) };

    int tmp = solveInput->enabledBitField;
    solveInput->enabledBitField = enabledBitField;
    evaluate_and_add_to(solveInput, 0, U_x, U_xx_);
    solveInput->enabledBitField = tmp;

    finite_difference_and_add_to(U_x_fd, solveInput, enabledBitField, fd_stepsize);
    finite_difference_and_add_to(U_xx_fd_, solveInput, enabledBitField, fd_stepsize);
    for_(k, arrlen(U_xx_)) U_xx(U_xx_[k].row, U_xx_[k].col) += U_xx_[k].val;
    for_(k, arrlen(U_xx_fd_)) U_xx_fd(U_xx_fd_[k].row, U_xx_fd_[k].col) += U_xx_fd_[k].val;

    arrfree(U_xx_);
    arrfree(U_xx_fd_);

    {
        bool PASSES = true;
        for_(k, LEN_X) {
            real a = U_x[k];
            real b = U_x_fd[k];
            real error = fabs(a - b);
            bool INVALID = TRIGGER_INVALID(a) || TRIGGER_INVALID(b);
            bool WRONG = TRIGGER_ERROR(error);
            if (INVALID || WRONG) {
                if (PASSES) { printf(" -- FAIL\n"); }
                PASSES = false;
                if (INVALID) { printf("%3d: nan or inf\n", k); }
                else { printf("%2d: | (%lf) - (%lf) | = %lf\n", k, a, b, error); }
            }
        }
        if (PASSES) { printf(" -- PASS\n"); }
    }
    {
        bool PASSES = true;
        for_(r, LEN_X) for_(c, LEN_X) {
            real a = U_xx(r, c);
            real b = U_xx_fd(r, c);
            real error = fabs(a - b);
            bool INVALID = TRIGGER_INVALID(a) || TRIGGER_INVALID(b);
            bool WRONG = TRIGGER_ERROR(error);
            if (INVALID || WRONG) {
                if (PASSES) { printf(" -- FAIL\n"); }
                PASSES = false;
                if (INVALID) { printf("%3d, %3d: nan or inf\n", r, c); }
                else { printf("%2d, %2d: | (%lf) - (%lf) | = %lf\n", r, c, a, b, error); }
            }
        }
        if (PASSES) { printf(" -- PASS\n"); }
    }



    free(U_x);
    free(U_x_fd);
    free(U_xx.data);
    free(U_xx_fd.data);
}
void finite_difference_and_add_to(real *U_x, SolveInput *solveInput, int enabledBitField, real fd_stepsize = 1e-5) {
    ASSERT(U_x);
    int push = solveInput->enabledBitField;
    solveInput->enabledBitField = enabledBitField;
    {

        real left, right;

        for_(i, LEN_X) {

            real x_i_0 = solveInput->x[i]; {

                left = 0;
                solveInput->x[i] = x_i_0 - fd_stepsize;
                evaluate_and_add_to(solveInput, &left, 0, 0);

                right = 0;
                solveInput->x[i] = x_i_0 + fd_stepsize;
                evaluate_and_add_to(solveInput, &right, 0, 0);

            } solveInput->x[i] = x_i_0;

            U_x[i] += (right - left) / (2 * fd_stepsize);

        }

    }
    solveInput->enabledBitField = push;
}
void finite_difference_and_add_to(SparseMatrixEntry *&U_xx, SolveInput *solveInput, int enabledBitField, real fd_stepsize = 1e-5) {
    ASSERT(U_xx);
    int push = solveInput->enabledBitField;
    solveInput->enabledBitField = enabledBitField;
    {

        real *left = (real *) malloc(LEN_X * sizeof(real));
        real *right = (real *) malloc(LEN_X * sizeof(real));

        for_(c, LEN_X) {

            real x_c_0 = solveInput->x[c]; {

                memset(left, 0, LEN_X * sizeof(real));
                solveInput->x[c] = x_c_0 - fd_stepsize;
                evaluate_and_add_to(solveInput, 0, left, 0);

                memset(right, 0, LEN_X * sizeof(real));
                solveInput->x[c] = x_c_0 + fd_stepsize;
                evaluate_and_add_to(solveInput, 0, right, NO_HESSIAN);

            } solveInput->x[c] = x_c_0;

            for_(r, LEN_X) {
                real value = (right[r] - left[r]) / (2 * fd_stepsize);
                if (!IS_ZERO(value)) arrput(U_xx, { r, c, value } );
            }

        }

        free(left);
        free(right);

    }
    solveInput->enabledBitField = push;
}



else {
    { // fun stuff
        static bool send_sinusoids = false; gui_checkbox("send_sinusoids", &send_sinusoids, '1');
        if (send_sinusoids) {
            for_(j, sim.num_cables) {
                if (j < sim.num_cables / 3) {
                    currentState.u[j] = sin((j * 1.77) + time) * (ROBOT_LENGTH / 6);
                } else {
                    currentState.u[j] = -100.0;
                }
            }
        }

        static bool move_pins = false; gui_checkbox("move_pins", &move_pins, '2');
        if (move_pins) {
            for_(i, sim.num_pins) {
                ((vec3 *) currentState.x_pin.data)[i] += V3(.1 * sin(time), 0, 0);
            }
        }
    }

    currentState = sim.getNext(&currentState);
}

// // OPTIMIZATION
// 65% solving physics
// - 29% evaluate and add to
// 02% computing the hessian again (start here)
// 26% the getNextX for dxdu -- adjoint method
// hide the gui

// // MODEL 
// TODO: experiment with flips of currentState thec:
// TODO: experiment with symmetric tetrahedralization
// TODO: scrub youngs modulus

// TODO: write down all the computation that happens in one frame from a high level
// compute G, H
// getNextX for du
// line search for alpha (or just YOLO it, assuming the GD line search basically always succeeds)
// u <- alpha * du
// getNextX for x

//       could also go from 3 3 3 to 2 4 4
// TODO: hack on gripper (could add long segment)
// TODO: holding right click



#if 0
{ // FORNOW: calling sim.draw
    static bool tabbed;
    gui_checkbox("tabbed", &tabbed, COW_KEY_TAB);

    { // draw
        if (!tabbed) {
            sim.draw(P, V, IdentityMatrix<4>(), &currentState);
        } else {
            sim.draw(P * V, &currentState);
        }

    }
}
#endif
