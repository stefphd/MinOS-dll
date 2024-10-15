/*  File ipoptsolver.cpp 
    NLP solver interface for IPOPT
    Copyright (C) 2024 Stefano Lovato
*/

#include "ipoptsolver.h" // base class with all method implementations
#include "macros.h" // for nlp macros

/**
* IPOPT solver class implementation
*/

/** Constructor */
IPOPTSolver::IPOPTSolver(
    OCPInterface* ocpInterface
) {
    this->ocpInterface = ocpInterface; // store the ocpInterface pointer
    lamzu = new double[ocpInterface->nz+1];
    lamzl = new double[ocpInterface->nz+1];
}

/** Destructor */
IPOPTSolver::~IPOPTSolver(){   
}

/** Free memory allocated */
void IPOPTSolver::free_mem() {
    // free mem
    delete[] lamzu;
    delete[] lamzl;
    this->ocpInterface = NULL;
}

bool IPOPTSolver::get_nlp_info(
   Index&          n,
   Index&          m,
   Index&          nnz_jac_g,
   Index&          nnz_h_lag,
   IndexStyleEnum& index_style
) {
   // Num of NLP vars and constraints
   n = ocpInterface->nz; 
   m = ocpInterface->ng;

   // NNZ of jac g
   nnz_jac_g = ocpInterface->nnzj;

   // NNZ of hess l
   nnz_h_lag = ocpInterface->nnzh;

   // Use the standard c index style for row/col entries
   index_style = C_STYLE;

   return true;
}

bool IPOPTSolver::get_bounds_info(
   Index   n,
   Number* zl,
   Number* zu,
   Index   m,
   Number* gl,
   Number* gu
) {
    ocpInterface->get_nlp_bounds((double*) zl, (double*) zu, (double*) gl, (double*) gu);
    return true;
}

bool IPOPTSolver::get_starting_point(
   Index   n,
   bool    init_z,
   Number* z,
   bool    init_lamz,
   Number* lamzl,
   Number* lamzu,
   Index   m,
   bool    init_lamg,
   Number* lamg
) { 
    double *lamz = NULL;
    if (init_lamz) lamz =  new double[n+1];
    ocpInterface->get_nlp_init(init_z, (double*) z,
                            init_lamz, (double*) lamz,
                            init_lamg, (double*) lamg);
    // Compute lamzl, lamzu using POS and NEG parts of lamz
    if (init_lamz && lamz) {
        for (int k = 0; k < n; ++k) {
            Number lamzk = lamz[k]; // temporary k-th z lambda (either > or < 0)
            lamzl[k] = NEGPART(lamzk); // negative part is lamzl
            lamzu[k] = POSPART(lamzk); // positive part is lamzk
        }
    }
    if (init_lamz) delete[] lamz;
    return true;
}

bool IPOPTSolver::eval_f(
   Index         n,
   const Number* z,
   bool          new_z,
   Number&       obj
) {
    return (ocpInterface->eval_obj((const double*) z, (double*) &obj) == 0);
}

bool IPOPTSolver::eval_grad_f(
   Index         n,
   const Number* z,
   bool          new_z,
   Number*       grad
) {
    return (ocpInterface->eval_obj_grad((const double*) z, (double*) grad) == 0);
}

bool IPOPTSolver::eval_g(
   Index         n,
   const Number* z,
   bool          new_z,
   Index         m,
   Number*       g
) {
    return (ocpInterface->eval_constr((const double*) z, (double*) g) == 0);
}

bool IPOPTSolver::eval_jac_g(
   Index         n,
   const Number* z,
   bool          new_z,
   Index         m,
   Index         nnz,
   Index*        ir,
   Index*        jc,
   Number*       jac
) {
    // calc jacobian indexes
    if(ir && jc) {
        return (ocpInterface->get_pattern_jac(ir, jc) == nnz); // check also the nnz returned by get_pattern_jac
    }

    // calc jacobian values
    if (jac) {
        return (ocpInterface->eval_constr_jac((const double*) z, (double*) jac) == 0);
    }

    return true;
}

bool IPOPTSolver::eval_h(
    Index         n,
    const Number* z,
    bool          new_z,
    Number        sigma,
    Index         m,
    const Number* lamg,
    bool          new_lamg,
    Index         nnz,
    Index*        ir,
    Index*        jc,
    Number*       hess
) {
    // calc hessian indexes
    if(ir && jc) {
        return (ocpInterface->get_pattern_hess(ir, jc) == nnz); // check also the nnz returned by get_pattern_hess
    }

    // calc hessian values
    if (hess) {
        return (ocpInterface->eval_hessian((const double*) z, (const double) sigma, (const double*) lamg, (double*) hess) == 0);
    }

    return true;
}

void IPOPTSolver::finalize_solution(
   SolverReturn               status,
   Index                      n,
   const Number*              z,
   const Number*              lamzl,
   const Number*              lamzu,
   Index                      m,
   const Number*              g,
   const Number*              lamg,
   Number                     obj_value,
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq
) {
    // copy solution
    ocpInterface->J_opt = obj_value;
    memcpy(ocpInterface->z_opt, z, ocpInterface->nz*sizeof(double));
    for (int i = 0; i < ocpInterface->nz; ++i) ocpInterface->lamz_opt[i] = lamzu[i] - lamzl[i];
    memcpy(ocpInterface->lamg_opt, lamg, ocpInterface->ng*sizeof(double));
    // eval nlp functions at solution
    ocpInterface->eval_ocpfuncs();
    // solution data
    ocpInterface->mu_curr = ip_data->curr_mu();
    ocpInterface->num_iter = ip_data->iter_count();
    // set guess solution to optimal if success
    if ((status == SUCCESS) || (status == STOP_AT_ACCEPTABLE_POINT))
        ocpInterface->set_optsol_as_guess();

    // save cpu time
    ocpInterface->tcpu_tot = ip_data->TimingStats().OverallAlgorithm().TotalCpuTime();
    ocpInterface->tcpu_eval = ip_data->TimingStats().TotalFunctionEvaluationCpuTime();
}

bool IPOPTSolver::intermediate_callback(
   AlgorithmMode              mode,
   Index                      iter,
   Number                     obj_value,
   Number                     inf_pr,
   Number                     inf_du,
   Number                     mu,
   Number                     d_norm,
   Number                     regularization_size,
   Number                     alpha_du,
   Number                     alpha_pr,
   Index                      ls_trials,
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq
) {
    // avoid repeating iter
    if (iter<=ocpInterface->num_iter) { return true; }
    // eval current solution data
    if ((ocpInterface->print_itersol>0) && ((iter % ocpInterface->print_itersol) == 0)) {
        Ipopt::TNLP::get_curr_iterate(ip_data, ip_cq, false, ocpInterface->nz, ocpInterface->z_opt, lamzl, lamzu, ocpInterface->ng, ocpInterface->g_opt, ocpInterface->lamg_opt);
        for (int i = 0; i < ocpInterface->nz; ++i) ocpInterface->lamz_opt[i] = lamzu[i] - lamzl[i];
    }    
    // print iter
    return ocpInterface->print_iter(iter, obj_value, inf_pr, inf_du, false); // false to not force print
}

/** IPOPTSolver::callSolve */
int IPOPTSolver::callSolve(
    OCPInterface *ocpInterface
){
    /* Create an instance of the IpoptApplication
       We are using the factory, since this allows us to compile this
       example with an Ipopt Windows DLL */
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    
    /* IPOPT options */
    SET_IPOPT_STR_OPTION(app, linear_solver, "mumps"); // use mumps linear solver
    SET_IPOPT_STR_OPTION(app, print_timing_statistics, "yes");
    SET_IPOPT_INT_OPTION(app, print_level, 0);
    SET_IPOPT_STR_OPTION(app, sb, "yes"); // suppress ipopt banner
    SET_IPOPT_INT_OPTION(app, file_print_level, 5);
    //if (ocpInterface->logfile.empty()) ocpInterface->logfile = "ipopt.log";
    if (ocpInterface->logfile != "none")
        SET_IPOPT_STR_OPTION(app, output_file, ocpInterface->logfile);
    else 
        SET_IPOPT_STR_OPTION(app, output_file, "");

    if (ocpInterface->check_lambda_guess())
        SET_IPOPT_STR_OPTION(app, warm_start_init_point, "yes");
    if (ocpInterface->max_iter<0) 
        ocpInterface->max_iter = MAX_DEFAULT_ITER; // use default if <0
    SET_IPOPT_INT_OPTION(app, max_iter, ocpInterface->max_iter);
    if (ocpInterface->flag_hessian)
        SET_IPOPT_STR_OPTION(app, hessian_approximation, "limited-memory");
    else   
        SET_IPOPT_STR_OPTION(app, hessian_approximation, "exact"); 
    if (ocpInterface->mu_init>0) // set mu_init if >0
        SET_IPOPT_NUM_OPTION(app, mu_init, ocpInterface->mu_init);
    /* Load IPOPT option file if any */
    // IPOPT read "ipopt.opt" by default, however check still necessary to notify user
    const char optfile[] = "ipopt.opt";
    if (ocpInterface->check_nlpopt_file(optfile))
        SET_IPOPT_STR_OPTION(app, option_file_name, optfile);

    /* Init history */
    ocpInterface->init_hist_mem(ocpInterface->max_iter);

    /* Initialize the IpoptApplication and process the options */
    ApplicationReturnStatus status;
    status = app->Initialize();
    if( status != Solve_Succeeded ) { return (int) status; } // sth wrong

    /* Create IPOPTSolver object */
    Ipopt::SmartPtr<IPOPTSolver> ipoptSolver = new IPOPTSolver(ocpInterface);

    /* Call IPOPT solve */
    status = app->OptimizeTNLP(ipoptSolver);

    /* Print last iter */
    Ipopt::SmartPtr<Ipopt::IpoptCalculatedQuantities> ip_cq = app->IpoptCQObject();
    double inf_pr = ip_cq->curr_primal_infeasibility(Ipopt::ENormType::NORM_MAX);
    double inf_du = ip_cq->curr_dual_infeasibility(Ipopt::ENormType::NORM_MAX);
    if (ocpInterface->display)
        ocpInterface->print_iter(ocpInterface->num_iter, ocpInterface->J_opt, inf_pr, inf_du, true);

    /* Print exit status */
    if (ocpInterface->display)
        ocpInterface->print_stats(status, Solve_Succeeded, Solved_To_Acceptable_Level, User_Requested_Stop, Maximum_Iterations_Exceeded, Infeasible_Problem_Detected);

    /* Explicitly free memory */
    ipoptSolver->free_mem();

    /* Return status */
    return (int) status;
}

/** A wrapper around IPOPTSolver::callSolve to load the static function from the dynamic library in runtime */
#ifdef __cplusplus
extern "C" {
#endif
MINOS_EXPORT_API int callSolve(
    OCPInterface* ocpInterface
) {
    return IPOPTSolver::callSolve(ocpInterface);
}
#ifdef __cplusplus
} /* extern "C" */
#endif
