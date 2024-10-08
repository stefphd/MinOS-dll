/*  File minosc.h 
    C interface to OCPInterface class
    Copyright (C) 2024 Stefano Lovato
*/

#include <stdbool.h> // for bool

#define MINOSC_PREFIX(x) minos_ ## x
#ifndef MINOS_EXPORT_API
    #ifdef _WIN32  // For Windows
        #ifdef MAKE_MINOS
            #define MINOS_EXPORT_API __declspec(dllexport)
        #else
            #define MINOS_EXPORT_API __declspec(dllimport)
        #endif
    #elif defined(__linux__) || defined(__APPLE__)  // For Linux and macOS
        #ifdef MAKE_MINOS
            #define MINOS_EXPORT_API __attribute__((visibility("default")))
        #else
            #define MINOS_EXPORT_API
        #endif
    #else
        #define MINOS_EXPORT_API // Fallback for other platforms
    #endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** @defgroup cinterface C Interface
 *  C Interface to OCPInterface class.
 *  @{
 */

/** 
 * \brief OCPInterface C type.
 * 
 * C type for OCPInterface instance.
 */
typedef void* OCP_t;

/**
 * \brief Option keywords.
 * 
 * Keywords to set option for the solver. Use OCPInterface::set_option to set an option using OCPIneterface::Opt as a option key.
 * - `MAX_ITER` (value): maximum number of iterations (default 3000).
 * - `MU_INIT` (value): initial NLP barrier parameter (default depends on the NLP solver). This applies only when usng interior-point and may be used for a warm start of the NLP solver.
 * - `FLAG_HESSIAN` (value): true to use approximated Hessian instead of exact one (default false). This option has no effect if the optimal-control problem is built skipping the Hessian calculation.
 * - `DISPLAY` (value): true to display the solver output messages (default true).
 * - `PRINT_ITERSOL` (value): print output file every `PRINT_ITERSOL` iterations, 0 to disable (default 0).
 * - `LOGFILE` (string): log file name for the NLP solver (default is `<name>.log`). Use an empty string for default and "none" for no log file.
 * - `NLPSOLVER` (string): NLP solver to use (default is `ipopt`)
 */
enum Opt {
    MAX_ITER,
    MU_INIT,
    FLAG_HESSIAN,
    DISPLAY,
    PRINT_ITERSOL,
    LOGFILE,
    NLPSOLVER
};

/**
 * @brief Instantiates a new OCP.
 * 
 * Instantiates a new OCP with `N` mesh points, initial time `ti`, 
 * and final time `tf`. The default mesh consists in an equally-spaced mesh grid.
 * The OCP functions are loaded in runtime from a dynamic library `<name>`; returns
 * a non-zero value if fails.
 * 
 * \param ocp_ptr Pointer to OCP instance.
 * \param name Problem name (i.e. dynamic library name implementing OCP functions, without extension).
 * \param N Number of mesh points.
 * \param ti Initial time.
 * \param tf Final time.
 * 
 * \return Non-zero value in the case of failure.
 */
MINOS_EXPORT_API int MINOSC_PREFIX(new)(
    OCP_t* ocp_ptr, 
    const char* name,
    int N,
    double ti,
    double tf
);

/**
 * @brief Instantiates a new OCP.
 * 
 * Instantiates a new OCP with `N` mesh points, initial time `ti`, 
 * and final time `tf`. The default mesh consists in an equally-spaced mesh grid.
 * The OCP functions are loaded in runtime from a dynamic library `<name>`; returns
 * the OCP instance.
 * 
 * \param name Problem name (i.e. dynamic library name implementing OCP functions, without extension).
 * \param N Number of mesh points.
 * \param ti Initial time.
 * \param tf Final time.
 * 
 * \return The OCP instance (NULL if failed).
 */
MINOS_EXPORT_API OCP_t MINOSC_PREFIX(new2)(
    const char* name,
    int N,
    double ti,
    double tf
);

/**
 * @brief Deletes OCP.
 * Deletes OCP and free allocated memory.
 * 
 * \param ocp_ptr Pointer to OCP instance.
 */
MINOS_EXPORT_API void MINOSC_PREFIX(free)(
    OCP_t* ocp_ptr
);

/** 
 * \brief Gets the problem dimensions.
 * 
 * Get the optimal-control problem dimensions, such as number of states, controls, 
 * parameters, path constraints, boundary conditions, etc.
 * Use NULL to not retrieve a result.
 * 
 * \param ocp OCP instance.
 * \param nx Number of states.
 * \param nu Number of controls.
 * \param np Number of static parameters.
 * \param nc Number of path constraints.
 * \param nb Number of boundary conditions.
 * \param nq Number of integral constraints.
 * \param nz Number of NLP variables.
 * \param ng Number of NLP constraints.
 * \param nnzj Number of non-zeros elements in NLP Jacobvian.
 * \param nnzh Number of non-zeros elements in NLP Lagragian Hessian.
 * \param na Number of auxdata.
 */
MINOS_EXPORT_API void MINOSC_PREFIX(get_dims)(
    const OCP_t ocp,
    int *nx,
    int *nu,
    int *np,
    int *nc,
    int *nb,
    int *nq,
    int *nz,
    int *ng,
    int *nnzj,
    int *nnzh,
    int *na
);

/**
 * \brief Gets the number of mesh points.
 * 
 * Returns the number of mesh points `N`.
 * 
 * \param ocp OCP instance.
 * 
 * \return The number of mesh points.
 */
MINOS_EXPORT_API int MINOSC_PREFIX(get_N)(
    const OCP_t ocp
);


/**
 * \brief Sets the guess.
 * 
 * Sets the guess for the optimal-control problem. Guess for the Lagrange multipliers may be optionally provided. 
 * Note that all Lagrange multipliers should be specified for a warm start of the NLP solving.
 * Use NULL to not specified a guess.
 * 
 * \param ocp OCP instance.
 * \param x0 State guess, specified as an array of length `N*nx`, with `x0[nx*k+i]` the i-th state at the k-th mesh point.
 * \param u0 Control guess, specified as an array of length `N*nu`, with `u0[nu*k+i]` the i-th control at the k-th mesh point.
 * \param p0 Static parameter guess, specified as an array of length `np`.
 * \param lam_x0 State multiplier guess, specified as an array of length `N*nx`, with `lam_x0[nx*k+i]` the i-th state multiplier at the k-th mesh point.
 * \param lam_u0 Control multiplier guess, specified as an array of length `N*nu`, with `lam_u0[nu*k+i]` the i-th control multiplier at the k-th mesh point.
 * \param lam_p0 Static parameter multiplier guess, specified as an array of length `np`.
 * \param lam_f0 Dynamic constraint multiplier guess, specified as an array of length `(N-1)*nx`, with `lam_f0[nx*k+i]` the i-th multiplier at the k-th mesh point.
 * \param lam_c0 Path constraint multiplier guess, specified as an array of length `N*nc`, with `lam_c0[nc*k+i]` the i-th multiplier at the k-th mesh point.
 * \param lam_b0 Boundary condition multiplier guess, specified as an array of length `nb`.
 * \param lam_q0 Integral constraint multiplier guess, specified as an array of length `nq`.
 */
MINOS_EXPORT_API void MINOSC_PREFIX(set_guess)(
    const OCP_t ocp,
    double *x0,
    double *u0,
    double *p0,
    double *lam_x0,
    double *lam_u0,
    double *lam_p0,
    double *lam_f0,
    double *lam_c0,
    double *lam_b0,
    double *lam_q0
);

/**
 * \brief Sets the bounds.
 * 
 * Sets the bounds for the optimal-control problem.
 * Use NULL to not specified a bounds.
 * 
 * \param ocp OCP instance.
 * \param lbx State lower bounds, specified as an array of length `nx`.
 * \param ubx State upper bounds, specified as an array of length `nx`.
 * \param lbu Control lower bounds, specified as an array of length `nu`.
 * \param ubu Control upper bounds, specified as an array of length `nu`.
 * \param lbp Static parameter lower bounds, specified as an array of length `np`.
 * \param ubp Static parameter upper bounds, specified as an array of length `np`.
 * \param lbc Path constraint lower bounds, specified as an array of length `nc`.
 * \param ubc Path constraint upper bounds, specified as an array of length `nc`.
 * \param lbb Boundary condition lower bounds, specified as an array of length `nb`.
 * \param ubb Boundary condition upper bounds, specified as an array of length `nb`.
 * \param lbq Integral constraint lower bounds, specified as an array of length `nq`.
 * \param ubq Integral constraint upper bounds, specified as an array of length `nq`.
 */
MINOS_EXPORT_API void MINOSC_PREFIX(set_bounds)(
    const OCP_t ocp,
    double *lbx,
    double *ubx,
    double *lbu,
    double *ubu,
    double *lbp,
    double *ubp,
    double *lbc,
    double *ubc,
    double *lbb,
    double *ubb,
    double *lbq,
    double *ubq
);

/**
 * \brief Sets a value option.
 * 
 * Sets a valued option for the solver. See Opt for a list of the available options.
 * 
 * \param ocp OCP instance.
 * \param optkey Keyword for the option to set.
 * \param val Value of the option.
 * 
 * \return Non-zero value in the case of errors.
 */
MINOS_EXPORT_API int MINOSC_PREFIX(set_option_val)(
    const OCP_t ocp,
    int optkey,
    double val
);

/** \brief Sets a string option.
 * 
 * Sets a string option for the solver. See Opt for a list of the available options.
 * 
 * \param ocp OCP instance.
 * \param optkey Keyword for the option to set.
 * \param str String of the option.
 * 
 * \return Non-zero value in the case of errors.
 */
MINOS_EXPORT_API int MINOSC_PREFIX(set_option)(
    const OCP_t ocp,
    int optkey,
    const char* str
);

/** 
 * \brief Sets the auxdata.
 * 
 * Sets the (optional) auxdata for the optimal-control problem.
 * 
 * \param ocp OCP instance.
 * \param auxdata Auxdata, specified as an array of length `na`.
 */
MINOS_EXPORT_API void MINOSC_PREFIX(set_auxdata)(
    const OCP_t ocp,
    double *auxdata
);

/**
 * \brief Solves the problem.
 * 
 * Solves the optimal-control problem.
 * 
 * \param ocp OCP instance.
 * 
 * \return Exit code of the employed NLP solver.
 * 
 */
MINOS_EXPORT_API int MINOSC_PREFIX(solve)(
    const OCP_t ocp
);


/**
 * \brief Gets the solution.
 * 
 * Gets the optimal solution of the optimal-control problem. 
 * Allocation of the required memory should be performed by the user according to the length of the arrays.
 * Use NULL to not retrieve a result.
 * 
 * \param ocp OCP instance.
 * \param J_opt Final cost.
 * \param t Time, specified as an array of length `N`.
 * \param x_opt State, specified as an array of length `N*nx`, with `x_opt[nx*k+i]` the i-th state at the k-th mesh point.
 * \param u_opt Control, specified as an array of length `N*nu`, with `u_opt[nu*k+i]` the i-th control at the k-th mesh point.
 * \param p_opt Static parameter, specified as an array of length `np`.
 * \param lamx_opt State multiplier, specified as an array of length `N*nx`, with `lamx_opt[nx*k+i]` the i-th state multiplier at the k-th mesh point.
 * \param lamu_opt Control multiplier, specified as an array of length `N*nu`, with `lamu_opt[nu*k+i]` the i-th control multiplier at the k-th mesh point.
 * \param lamp_opt Static parameter multiplier, specified as an array of length `np`.
 * \param lamf_opt Dynamic constraint multiplier, specified as an array of length `(N-1)*nx`, with `lamf_opt[nx*k+i]` the i-th multiplier at the k-th mesh point.
 * \param lamc_opt Path constraint multiplier, specified as an array of length `N*nc`, with `lamc_opt[nc*k+i]` the i-th multiplier at the k-th mesh point.
 * \param lamb_opt Boundary condition multiplier, specified as an array of length `nb`.
 * \param lamq_opt Integral constraint multiplier, specified as an array of length `nq`.
 * \param f_opt Dynamic constraint, specified as an array of length `(N-1)*nx`, with `f_opt[nx*k+i]` the i-th dynamic constraint at the k-th mesh point.
 * \param c_opt Path constraint, specified as an array of length `N*nc`, with `c_opt[nc*k+i]` the i-th path constraint at the k-th mesh point.
 * \param b_opt Boundary condition, specified as an array of length `nb`.
 * \param q_opt Integral constraint, specified as an array of length `nq`.
 * \param l_opt Lagrange cost, specified as an array of length `(N-1)`.
 * \param m_opt Mayer cost.
 * \param gradx_opt NLP cost gradient w.r.t. the state, specified as an array of length `N*nx`, with `gradx_opt[nx*k+i]` the gradient w.r.t. the i-th state at the k-th mesh point.
 * \param gradu_opt NLP cost gradient w.r.t. the control, specified as an array of length `N*nu`, with `gradu_opt[nu*k+i]` the gradient w.r.t. the i-th control at the k-th mesh point.
 * \param gradp_opt NLP cost gradient w.r.t. the static parameter, specified as an array of length `np`.
 * \param irj Row indexes of non-zero elements of NLP Jacobian, specified as an array of length `nnzj`. CSC format is used.
 * \param jcj Column indexes of non-zero elements of NLP Jacobian, specified as an array of length `nnzj`. CSC format is used.
 * \param jac Values of non-zero elements of NLP Jacobian, specified as an array of length `nnzj`.
 * \param irh Row indexes of non-zero elements of NLP Lagragian Hessian, specified as an array of length `nnzh`. CSC format is used.
 * \param jch Row indexes of non-zero elements of NLP Lagragian Hessian, specified as an array of length `nnzh`. CSC format is used.
 * \param hess Values of non-zero elements of NLP Lagragian Hessian, specified as an array of length `nnzh`.
 */
MINOS_EXPORT_API void MINOSC_PREFIX(get_sol)(
    const OCP_t ocp,
    double *J_opt, double *t,
    double *x_opt, double *u_opt, double *p_opt,
    double *lamx_opt, double *lamu_opt, double *lamp_opt,
    double *lamf_opt, double *lamc_opt, double *lamb_opt, double *lamq_opt,
    double *f_opt, double *c_opt, double *b_opt, double *q_opt,
    double *l_opt, double *m_opt,
    double *gradx_opt, double *gradu_opt, double *gradp_opt,
    size_t *irj, size_t *jcj, double *jac,
    size_t *irh, size_t *jch, double *hess
);

/**
 * \brief Outputs the solution.
 * 
 * Outputs the solution to a C string. Allocation is performed internally.
 * You can use this function to print the C string to a stream.
 * 
 * \param ocp OCP instance.
 * 
 * \return A C string.
 */
MINOS_EXPORT_API const char* MINOSC_PREFIX(print_sol)(
    const OCP_t ocp
);

/**
 * \brief Sets the mesh.
 * 
 * Sets the mesh fractions for the optimal-control problem. The default mesh consists in an equally-spaced mesh grid.
 * This option is usefull when mesh refinements is required.
 * 
 * \param ocp OCP instance.
 * \param mesh Mesh fractions, specify by as an array of length `N-1` that sums to unity.
 */
MINOS_EXPORT_API void MINOSC_PREFIX(set_mesh)(
    const OCP_t ocp,
    double *mesh
);

/**
 * \brief Gets the CPU times.
 * 
 * Gets the CPU times required to solve the optimal-control problem.
 * 
 * \param ocp OCP instance.
 * \param tcpu_tot Total CPU time in seconds.
 * \param tcpu_alg CPU time required for algorithm execution in seconds.
 * \param tcpu_eval CPU time required for NLP function evaluation in seconds.
 */
MINOS_EXPORT_API void MINOSC_PREFIX(get_cpu_time)(
    const OCP_t ocp,
    double* tcpu_tot,
    double* tcpu_alg,
    double* tcpu_eval
);

/**
 * \brief Gets the NLP barrier parameter.
 * 
 * Gets the current barrier parameter `mu` used by the NLP solver. This may be usefull for a warm start of the NLP solver by setting the option MU_INIT.
 * This applies only when using interior-point algorithms.
 * 
 * \param ocp OCP instance.
 * 
 * \return Current barrier parameter.
 */
MINOS_EXPORT_API double MINOSC_PREFIX(get_mu_curr)(
    const OCP_t ocp
);

/** 
 * \brief Gets the number of iterations.
 * 
 * Gets the number of iterations required to solve the optimal-control problem.
 * 
 * \param ocp OCP instance.
 * 
 * \return Number of iterations.
 */
MINOS_EXPORT_API int MINOSC_PREFIX(get_num_iter)(
    const OCP_t ocp
);

/** 
 * \brief Gets the convergence history.
 * 
 * Gets the convergence history of the objective, infinite norm of constraint violation, and infinite norm of optimality.
 * Allocation of the required memory should be performed by the user according to the length of the arrays.
 * 
 * \param ocp OCP instance.
 * \param obj Convergence history for objective, specified as an array of length returned by `get_num_iter()+1`.
 * \param inf_pr Convergence history for constraint violation, specified as an array of length returned by `get_num_iter()+1`.
 * \param inf_du Convergence history for optimality, specified as an array of length returned by `get_num_iter()+1`.
 */
MINOS_EXPORT_API void MINOSC_PREFIX(get_history)(
    const OCP_t ocp,
    double* obj,
    double* inf_pr,
    double* inf_du
);

/**
 * \brief Gets MinOS version.
 * 
 * Gets the current version of MinOS.
 * 
 * \return MinOS version in the format X.Y.Z (X=major version, Y=minor version, Z=release number).
 */
MINOS_EXPORT_API const char* MINOSC_PREFIX(get_version)(
);

/**
 * \brief Sets the printing function.
 * 
 * Sets the printing function to be used by the solver to print the output messages. 
 * The default one is the standard `printf`. 
 * Other functions may be specified by the user for a custom printing (e.g. to a file etc.).
 * 
 * \param ocp OCP instance.
 * \param ext_print_funptr Pointer to the printing function, with prototype `int (const char* fmt, ...)`
 */
MINOS_EXPORT_API void MINOSC_PREFIX(set_printfun)(
    const OCP_t ocp,
    int (*ext_print_funptr)(const char *fmt, ...)
);

/** 
 * \brief Sets the interrupt handling function.
 * 
 * Sets the interrupt handling function which is checked by the solver every iterations.
 * The default one is a function that catches the CTRL+C interrupt and stops the solving.
 * Other functions may be specified to change the default behaviour.
 * 
 * \param ocp OCP instance.
 * \param ext_int_funptr Pointer to the interrupt handling function, with prototype `bool (void)`. This function returns true when required to stop the solving.
 */
MINOS_EXPORT_API void MINOSC_PREFIX(set_interruptfun)(
    const OCP_t ocp,
    bool (*ext_int_funptr)(void)
);

#ifdef __cplusplus
} /* extern "C" */
#endif

/** @} */