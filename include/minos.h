/*  File minos.h 
    Declarations of OCPInterface class
    Copyright (C) 2024 Stefano Lovato
*/

#ifndef _MINOS_H
#define _MINOS_H

#include <iostream> // for std::cout
#include <ostream> // for std::ostream
#include <sstream> // for std::ostringstream 
#include <vector> // for std::vector
#include <string> // for std::string
#include <cstring> // for memcpy
#include <cmath> // for NAN
#include <fstream> // for std::ifstream 
#include <sstream> // std::ostringstream
#include <cstdio> // for standard printf
#include <atomic> // for std:atomic_bool

#ifndef MINOS_EXPORT_API
    #ifdef _WIN32  // For Windows
        #ifdef BUILD_MINOS
            #define MINOS_EXPORT_API __declspec(dllexport)
        #else
            #define MINOS_EXPORT_API __declspec(dllimport)
        #endif
    #elif defined(__linux__) || defined(__APPLE__)  // For Linux and macOS
        #ifdef BUILD_MINOS
            #define MINOS_EXPORT_API __attribute__((visibility("default")))
        #else
            #define MINOS_EXPORT_API
        #endif
    #else
        #define MINOS_EXPORT_API // Fallback for other platforms
    #endif
#endif

/** 
 * \brief Sets and solves the optimal-control problem.
 * 
 * This class provides an interfaces to the user for an easy setting and solving 
 * of the optimal-control problem.
 */
class MINOS_EXPORT_API OCPInterface {

private:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    /* Friend class for IPOPT interface */
    friend class IPOPTSolver;

    /* Friend class for KNITRO interface */
    friend class KNITROSolver;

    /* Friend class for WORHP interface */
    friend class WORHPSolver;

    /* Friend class for SNOPT interface */
    friend class SNOPTSolver;

    /* Friend */

    /* Typedefs for OCP functions */
    typedef int         (*EvalFunc)     (const double**, double**, int*, double*, int);
    typedef int         (*AllocFunc)    (void);
    typedef void        (*FreeFunc)     (int);
    typedef int         (*NInFunc)      (void);
    typedef const int*  (*SpInFunc)     (int);
    typedef const int*  (*SpOutFunc)    (int);
    typedef int         (*WorkFunc)     (int*, int*, int*, int*);

    /* Structure to hold the OCP function data */
    struct ocpfun_t
    {
        EvalFunc eval = NULL;
        AllocFunc alloc = NULL;
        FreeFunc free = NULL;
        NInFunc nin = NULL;
        SpInFunc spin = NULL;
        SpOutFunc spout = NULL;
        WorkFunc work = NULL;
        int mem = 0;
        const double **arg = NULL;
        double **res = NULL;
        int *iw = NULL;
        double *w = NULL;
    };

    /* Typedef for NLP solve function */
    typedef int         (*SolveFunc)    (OCPInterface*);

    /* Friend function for handling CTRL+C */
    friend void signal_callback_handler(
        int signum
    );
#endif

    /* Name of the problem */
    std::string name;

    /* Dimensions */
    int N, nx, nu, np, nc, nb, nq, nz, ng, na; 

    /* Initial time, final time, time step, current time, time step, and mesh fractions */
    double ti, tf, t, h, *mesh;

    /* Cost gradient pattern */
    int nnzrcg, *irrcg = NULL, *krcg = NULL;
    int nnzbcg, *irbcg = NULL, *kbcg = NULL;

    /* Constraint jacobian pattern */
    int nnzdj, *irdj = NULL, *jcdj = NULL; // dyn jac
    int nnzpj, *irpj = NULL, *jcpj = NULL; // path jac
    int nnzbj, *irbj = NULL, *jcbj = NULL; // bcs jac
    int nnzqj, *irqj = NULL, *jcqj = NULL, *kjq = NULL, *kj = NULL; // int jac
    int nnzj; // nnz of constraint jac

    /* Lagragian hessian pattern */
    int nnzhb, nnzhi, nnzh; // hessian nnz
    int *irhb = NULL, *jchb = NULL, *irhi = NULL, *jchi = NULL; // hessian patterns
    int *kh = NULL, *khb = NULL, *khi = NULL; // hessian indexes
    
    /* Variable bounds */
    double *lbx = NULL, *ubx = NULL; 
    double *lbu = NULL, *ubu = NULL; 
    double *lbp = NULL, *ubp = NULL; 
    double *lbc = NULL, *ubc = NULL;
    double *lbb = NULL, *ubb = NULL;
    double *lbq = NULL, *ubq = NULL;

    /* Initial guess */
    double *x0 = NULL, *u0 = NULL, *p0 = NULL;
    double *lam_x0 = NULL, *lam_u0 = NULL, *lam_p0 = NULL;
    double *lam_f0 = NULL, *lam_c0 = NULL, *lam_b0 = NULL, *lam_q0 = NULL;

    /* auxdata */
    double *auxdata = NULL;
    
    /* Solution */
    double *z_opt = NULL, *lamz_opt = NULL, *lamg_opt = NULL;
    double J_opt = 0, *grad_opt = NULL, *g_opt = NULL, *jac_opt = NULL, *hess_opt = NULL;
    double mu_curr = -1, *obj_history = NULL, *infpr_history = NULL, *infdu_history = NULL;
    int num_iter = -1;

    /* Timing stats */
    double tcpu_tot = NAN, tcpu_eval = NAN;

    /* Options */
    int max_iter = -1; // max iterations, -1 for default value
    double mu_init = -1; // initil barrier parameter, -1 for default value
    bool flag_hessian = false; // flag hessian (true for approx Hessian)
    bool display = true; // flag display
    int print_itersol = 0; // flag print sol to file every printsol iters (<=0 for disable)
    std::string logfile; // NLP logfile name
    std::string nlpsolver = "ipopt";  // NLP solver name

    /* Flags for lambda guess */
    bool flag_lamx, flag_lamu, flag_lamp, flag_lamf, flag_lamc, flag_lamb, flag_lamq;

    /* Constant strings */
    static constexpr char print_headerFormat[] = "%-10s %-18s %-18s %-18s\n"; // using printf
    static constexpr char print_dataFormat[] = "%-10d %-+18.6e %-+18.6e %-+18.6e\n"; // using printf
    static constexpr char print_strFormat[] = "%-25s %s"; // using printf
    static constexpr char print_intFormat[] = "%-25s %d"; // using printf
    static constexpr char print_doubleFormat[] = "%-25s %-+12.5e"; // using printf
    static constexpr char print_indexFormat[] = "%3d"; // using printf
    static constexpr char print_rangeFormat[] = "%-25s [%-+12.5e, %-+12.5e]"; // using printf
    static constexpr char print_varnameFormat[] = "%-12s"; // using printf
    static constexpr char print_vardataFormat[] = "%-+12.5e"; // using printf

    /* Pointer to printing function */
    int (*print_funptr)(const char *fmt, ...) = &printf; // use standard printf by default

    /* Pointer to interrupt handling function */
    bool (*int_funptr)() = &OCPInterface::check_interrupt; // use internal check_interrupt by default

    /* Pointer to OCP libraries */
    void* ocp_lib = NULL;

    /* Pointers to OCP functions */
    ocpfun_t ocp_dyn, ocp_path, ocp_bcs, ocp_int, ocp_runcost, ocp_bcscost;
    ocpfun_t ocp_dyn_jac, ocp_path_jac, ocp_bcs_jac, ocp_int_jac, ocp_runcost_grad, ocp_bcscost_grad;
    ocpfun_t ocp_hessb, ocp_hessi;

public:

    /** 
     * \brief Creates a new OCPInterface object.
     * 
     * Creates a new OCPInterface object with `N` mesh points, initial time `ti`, 
     * and final time `tf`. The default mesh consists in an equally-spaced mesh grid.
     * The OCP functions are loaded in runtime from a dynamic library `<name>`; an exception
     * is thrown if load fails. You may catch the exception with try-catch.
     * 
     * \param name Problem name (i.e. dynamic library name implementing OCP functions, without extension).
     * \param N Number of mesh points.
     * \param ti Initial time.
     * \param tf Final time.
     */
    OCPInterface(
        std::string name,
        int N, 
        double ti, 
        double tf
    );

    /**
     * \brief Destroys the OCPInterface object.
     * 
     * Destroys the OCPInterface object and frees the allocated memory.
     */
    ~OCPInterface();

    /** 
     * \brief Gets the problem dimensions.
     * 
     * Get the optimal-control problem dimensions, such as number of states, controls, 
     * parameters, path constraints, boundary conditions, etc.
     * Use NULL to not retrieve a result.
     * 
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
    void get_dims(
        int *nx,
        int *nu,
        int *np = NULL,
        int *nc = NULL,
        int *nb = NULL,
        int *nq = NULL,
        int *nz = NULL,
        int *ng = NULL,
        int *nnzj = NULL,
        int *nnzh = NULL,
        int *na = NULL
    );

    /**
     * \brief Gets the number of mesh points.
     * 
     * Returns the number of mesh points `N`.
     * 
     * \return The number of mesh points.
     */
    int get_N();
 
    /**
     * \brief Sets the guess.
     * 
     * Sets the guess for the optimal-control problem. Guess for the Lagrange multipliers may be optionally provided. 
     * Note that all Lagrange multipliers should be specified for a warm start of the NLP solving.
     * Use NULL to not specified a guess.
     * 
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
    void set_guess(
        double *x0,
        double *u0,
        double *p0 = NULL,
        double *lam_x0 = NULL,
        double *lam_u0 = NULL,
        double *lam_p0 = NULL,
        double *lam_f0 = NULL,
        double *lam_c0 = NULL,
        double *lam_b0 = NULL,
        double *lam_q0 = NULL
    );
    
    /**
     * \brief Sets the bounds.
     * 
     * Sets the bounds for the optimal-control problem.
     * Use NULL to not specified a bounds.
     * 
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
    void set_bounds(
        double *lbx,
        double *ubx,
        double *lbu,
        double *ubu,
        double *lbp = NULL,
        double *ubp = NULL,
        double *lbc = NULL,
        double *ubc = NULL,
        double *lbb = NULL,
        double *ubb = NULL,
        double *lbq = NULL,
        double *ubq = NULL
    );

    /** 
     * \brief Sets the auxdata.
     * 
     * Sets the (optional) auxdata for the optimal-control problem.
     * 
     * \param auxdata Auxdata, specified as an array of length `na`.
     */
    void set_auxdata(
        double *auxdata = NULL
    );

    /**
     * \brief Solves the problem.
     * 
     * Solves the optimal-control problem.
     * 
     * \return Exit code of the employed NLP solver.
     * 
     */
    int solve();
    
    /**
     * \brief Gets the solution.
     * 
     * Gets the optimal solution of the optimal-control problem. 
     * Allocation of the required memory should be performed by the user according to the length of the arrays.
     * Use NULL to not retrieve a result.
     * 
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
     * \param irj Row indexes of non-zero elements of NLP Jacobian, specified as an array of length `nnzj`. COO format is used.
     * \param jcj Column indexes of non-zero elements of NLP Jacobian, specified as an array of length `nnzj`. COO format is used.
     * \param jac Values of non-zero elements of NLP Jacobian, specified as an array of length `nnzj`.
     * \param irh Row indexes of non-zero elements of NLP Lagragian Hessian, specified as an array of length `nnzh`. COO format is used.
     * \param jch Row indexes of non-zero elements of NLP Lagragian Hessian, specified as an array of length `nnzh`. COO format is used.
     * \param hess Values of non-zero elements of NLP Lagragian Hessian, specified as an array of length `nnzh`.
     */
    void get_sol(
        double *J_opt, double *t,
        double *x_opt, double *u_opt, double *p_opt = NULL,
        double *lamx_opt = NULL, double *lamu_opt = NULL, double *lamp_opt = NULL,
        double *lamf_opt = NULL, double *lamc_opt = NULL, double *lamb_opt = NULL, double *lamq_opt = NULL,
        double *f_opt = NULL, double *c_opt = NULL, double *b_opt = NULL, double *q_opt = NULL,
        double *l_opt = NULL, double *m_opt = NULL,
        double *gradx_opt = NULL, double *gradu_opt = NULL, double *gradp_opt = NULL,
        int *irj = NULL, int *jcj = NULL, double *jac = NULL,
        int *irh = NULL, int *jch = NULL, double *hess = NULL
    );

    /**
     * \brief Prints the solution.
     * 
     * Prints the optimal solution of the optimal-control problem to an output stream.
     * 
     * \param os Output stream.
     * \param ocpInterface OCPInterface object to print.
     * 
     * \return Output stream.
     */
    friend MINOS_EXPORT_API std::ostream& operator<<(
        std::ostream& os, 
        OCPInterface& ocpInterface
    );

    /**
     * \brief Prints the solution.
     * 
     * Prints the optimal solution of the optimal-control problem to an output stream.
     * Equivalent to the operator <<.
     * 
     * \param os Output stream.
     */
    void print_sol(
        std::ostream &os = std::cout
    );

    /**
     * \brief Converts the solution to a string.
     * 
     * Converts the optimal solution of the optimal-control problem to a string.
     * 
     * \return String representation of the optimal solution.
     */
    std::string toString();

    /**
     * \brief Sets the mesh.
     * 
     * Sets the mesh fractions for the optimal-control problem. The default mesh consists in an equally-spaced mesh grid.
     * This option is usefull when mesh refinements is required.
     * 
     * \param mesh Mesh fractions, specify by as an array of length `N-1` that sums to unity.
     */
    void set_mesh(
        double *mesh = NULL
    );

    /**
     * \brief Gets the name.
     * 
     * Gets the name of the optimal-control problem.
     * 
     */
    std::string get_name();

    /**
     * \brief Gets the initial and final times.
     * 
     * Gets the initial and final times of the optimal-control problem. Use NULL to not retrieve a result.
     * 
     * \param ti Initial time.
     * \param tf Final time.
     * 
     */
    void get_t(
        double* ti = NULL,
        double* tf = NULL
    );

    /**
     * \brief Gets the auxdata.
     * 
     * Gets the auxdata of the optimal-control problem. 
     * Allocation of the required memory should be performed by the user according to the length of the array.
     * 
     * \param auxdata Auxdata.
     * 
     */
    void get_auxdata(
        double* auxdata = NULL
    );

    /**
     * \brief Gets the bounds.
     * 
     * Gets the bounds of the optimal-control problem. Use NULL to not retrieve a result.
     * Allocation of the required memory should be performed by the user according to the length of the array.
     * 
     * \param lbx State lower bounds.
     * \param ubx State upper bounds
     * \param lbu Control lower bounds.
     * \param ubu Control upper bounds.
     * \param lbp Static parameter lower bounds.
     * \param ubp Static parameter upper bounds.
     * \param lbc Path constraint lower bounds.
     * \param ubc Path constraint upper bounds.
     * \param lbb Boundary condition lower bounds.
     * \param ubb Boundary condition upper bounds.
     * \param lbq Integral constraint lower bounds.
     * \param ubq Integral constraint upper bounds.
     */
    void get_bounds(
        double *lbx = NULL,
        double *ubx = NULL,
        double *lbu = NULL,
        double *ubu = NULL,
        double *lbp = NULL,
        double *ubp = NULL,
        double *lbc = NULL,
        double *ubc = NULL,
        double *lbb = NULL,
        double *ubb = NULL,
        double *lbq = NULL,
        double *ubq = NULL
    );

    /**
     * \brief Gets the mesh.
     * 
     * Gets the mesh fractions of the optimal-control problem.
     * Allocation of the required memory should be performed by the user according to the length of the array.
     * 
     * \param mesh Mesh fractions, specify by as an array of length `N-1`.
     */
    void get_mesh(
        double *mesh
    );

    /**
     * \brief Gets the CPU times.
     * 
     * Gets the CPU times required to solve the optimal-control problem.
     * 
     * \param tcpu_tot Total CPU time in seconds.
     * \param tcpu_alg CPU time required for algorithm execution in seconds.
     * \param tcpu_eval CPU time required for NLP function evaluation in seconds.
     */
    void get_cpu_time(
      double& tcpu_tot,
      double& tcpu_alg,
      double& tcpu_eval
    );

    /**
     * \brief Gets the CPU times.
     * 
     * Gets the CPU times required to solve the optimal-control problem.
     * 
     * \param tcpu_tot Total CPU time in seconds.
     * \param tcpu_alg CPU time required for algorithm execution in seconds.
     * \param tcpu_eval CPU time required for NLP function evaluation in seconds.
     */
    void get_cpu_time(
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
     * \return Current barrier parameter.
     */
    double get_mu_curr();

    /** 
     * \brief Gets the number of iterations.
     * 
     * Gets the number of iterations required to solve the optimal-control problem.
     * 
     * \return Number of iterations.
     */
    int get_num_iter();

    /** 
     * \brief Gets the convergence history.
     * 
     * Gets the convergence history of the objective, infinite norm of constraint violation, and infinite norm of optimality.
     * Allocation of the required memory should be performed by the user according to the length of the arrays.
     * 
     * \param obj Convergence history for objective, specified as an array of length returned by `get_num_iter()+1`.
     * \param inf_pr Convergence history for constraint violation, specified as an array of length returned by `get_num_iter()+1`.
     * \param inf_du Convergence history for optimality, specified as an array of length returned by `get_num_iter()+1`.
     */
    void get_history(
        double* obj,
        double* inf_pr,
        double* inf_du
    );

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
     * \brief Sets a value option.
     * 
     * Sets a valued option for the solver. See OCPInterface::Opt for a list of the available options.
     * 
     * \param optkey Keyword for the option to set, see also OCPInterface::Opt.
     * \param val Value of the option.
     * 
     * \return True for success.
     */
    bool set_option(
        int optkey,
        double val
    );

    /** \brief Sets a string option.
     * 
     * Sets a string option for the solver. See OCPInterface::Opt for a list of the available options.
     * 
     * \param optkey Keyword for the option to set.
     * \param str String of the option.
     * 
     * \return True for success.
     */
    bool set_option(
        int optkey,
        std::string str
    );

    /**
     * \brief Gets a value option.
     * 
     * Gets a valued option for the solver. See OCPInterface::Opt for a list of the available options.
     * 
     * \param optkey Keyword for the option to set, see also OCPInterface::Opt.
     * \param val Value of the option.
     */
    void get_option(
        int optkey,
        double* val
    );

    /** \brief Gets a string option.
     * 
     * Gets a string option for the solver. See OCPInterface::Opt for a list of the available options.
     * 
     * \param optkey Keyword for the option to set.
     * \param str String of the option.
     */
    void get_option(
        int optkey,
        std::string& str
    );

    /**
     * \brief Sets the printing function.
     * 
     * Sets the printing function to be used by the solver to print the output messages. 
     * The default one is the standard `printf`. 
     * Other functions may be specified by the user for a custom printing (e.g. to a file etc.).
     * 
     * \param ext_print_funptr Pointer to the printing function, with prototype `int (const char* fmt, ...)`
     */
    void set_printfun(
        int (*ext_print_funptr)(const char *fmt, ...)
    );

    /** 
     * \brief Sets the interrupt handling function.
     * 
     * Sets the interrupt handling function which is checked by the solver every iterations.
     * The default one is a function that catches the CTRL+C interrupt and stops the solving.
     * Other functions may be specified to change the default behaviour.
     * 
     * \param ext_int_funptr Pointer to the interrupt handling function, with prototype `bool (void)`. This function returns true when required to stop the solving.
     */
    void set_interruptfun(
        bool (*ext_int_funptr)(void)
    );

    /**
     * \brief Converts matrix from COO to CSC formats
     * 
     * Converts a sparse matrix from coordinate format (COO) to compressed sparse column (CSC) format.
     * 
     * \param n Numver of rows.
     * \param m Number of columns.
     * \param nnz Number of non-zero values.
     * \param ir1 Row indexes (COO).
     * \param jc1 Column indexes (COO).
     * \param v1 Non-zero values (COO).
     * \param ir2 Row indexes (CSC).
     * \param jc2 Column indexes (CSC).
     * \param v2 Non-zero values (CSC).
     */
    static void coo2csc(
        int n, int m, int nnz, 
        int *ir1, int *jc1, double *v1, // COO format
        size_t *ir2, size_t *jc2, double *v2  // CSC format
    );

private:

    /** Set optimal solution as guess */
    void set_optsol_as_guess(
    );

    /** Init cost gradient */
    void init_cost_gradient(
    );

    /** Init constraint jacobian */
    void init_constr_jac(
    );

    /** Init lagragian hessian */
    void init_lag_hessian(
    );

    /** Get pattern of constrain jacobian */
    int get_pattern_jac(
        int* irj,
        int* jcj
    );
    
    /** Get pattern of constrain jacobian - 64bit int version */
    long int get_pattern_jac(
        long int* irj,
        long int* jcj
    );
    
    /** Get pattern of lagragian hessian */
    int get_pattern_hess(
        int* irh,
        int* jch
    );
    
    /** Get pattern of lagragian hessian - 64bit int version */
    long int get_pattern_hess(
        long int* irh,
        long int* jch
    );

    /** Allocate memory */
    void alloc_mem();

    /** Init constant arguments */
    void init_args();

    /** Init memory to store data */
    void init_data_mem();

    /** Init memory for history */
    void init_hist_mem(
        int len
    );

    /** Check if guess for lambda is given */
    bool check_lambda_guess();

    /** Get NLP variable and constraint boundaries */
    void get_nlp_bounds(
        double* lbz,
        double* ubz,
        double* lbg,
        double* ubg
    );

    /** Get NLP initial point */
    void get_nlp_init(
        bool init_z,
        double* z,
        bool init_lamz,
        double* lamz,
        bool init_lamg,
        double* lamg
    );

    /** Eval objective */
    int eval_obj(
        const double* z, 
        double* obj
    );

    /** Eval objective gradient */
    int eval_obj_grad(
        const double* z, 
        double* grad
    );

    /** Eval constraint */
    int eval_constr(
        const double* z, 
        double* g
    );

    /** Eval constraint jacobian */
    int eval_constr_jac(
        const double* z, 
        double* jac
    );

    /** Eval lagragian hessian */
    int eval_hessian(
        const double* z, 
        const double sigma, 
        const double* lamg, 
        double* hess
    );

    /** Eval OCP functions */
    void eval_ocpfuncs(
    );

    /** Print iteration */
    bool print_iter(
        const int iter,
        const double obj,
        const double inf_pr,
        const double inf_du,
        const bool force
    );

    /** print statistics */
    void print_stats(
        const int status,
        const int succed,
        const int acceptable,
        const int userstop,
        const int maxiter,
        const int infeasible
    );

    /** Check NLP option file */
    bool check_nlpopt_file(
        const char* filename
    );

    /** Load OCP library */
    void* load_ocplib(
        std::string name
    );

    /** Load NLP library */
    void* load_nlplib(
        std::string name
    );

    /** Import NLP solve function */
    void* import_nlpsolve(
        void* libhandle
    );

    /** Free the specified library */
    void free_library(
        void* libhandle
    );

    /**@name Static methods for utilities */
    //@{

    /** Get n, m, and nnz of function input/output */
    static void get_sizes(
        const int* (*sparsity_fun) (int), 
        const int id,
        int *n,
        int *m,
        int *nnz
    );

    /** Get pattern indexes of function input/output */
    static void get_pattern(
        const int* (*sparsity_fun) (int), 
        int id,
        int *ir,
        int *jc
    );

    /** Allocate CASADI memory */
    static void alloc_casadi_mem(
        ocpfun_t *ocpfun
    );

    /** Deallocate CASADI memory */
    static void dealloc_casadi_mem(
        ocpfun_t *ocpfun
    );

    /** Search value in array */
    static int search_value(
        int arr[],
        int size,
        int value
    );

    /** Check for file existing */
    static bool exist_file(
        const char *filename
    );

    /** Implementation of format strings */
    static std::string format_str(
        const char *fmt
        , ...
    );

    /** Get header for MinOS */
    static std::string get_header();

    /** Check interrutpt during NLP solution.
     * Return true if interrupt detected, false otherwise. 
     * Currently, this function catches CTRL+C events to stop the NLP solving.
     * The function may be overridded using set_interruptfun.
     */
    static bool check_interrupt();

    /** Get mapping from unsorted to sorted indexes */
    static void get_sorted_mapping(
        int *arr, 
        int *id, 
        int n
    );

    /** Global flag that is true when CTRL+C detected. 
     * This is initialized to false every time OCPInterface::solve is called.
     */
    static std::atomic_bool is_interrupt_requested;

    //@}

};

#endif
