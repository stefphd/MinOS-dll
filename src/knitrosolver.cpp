/*  File knitrosolver.cpp 
    NLP solver interface for KNITRO
    Copyright (C) 2024 Stefano Lovato
*/

#define MAKE_MINOS

#include "knitrosolver.h" // base class with all method implementations
#include "macros.h" // for nlp macros

/** KNITROSolver implementation */
int KNITROSolver::callbackEvalFC (
    KN_context_ptr             kc,
    CB_context_ptr             cb,
    KN_eval_request_ptr const  evalRequest,
    KN_eval_result_ptr  const  evalResult,
    void              * const  ocpInterface_p
) {
    OCPInterface* ocpInterface = CASTOCPINTERFACE(ocpInterface_p);
    if (evalRequest->type != KN_RC_EVALFC) { return -1; }
    int exit;
    // Eval objective
    exit = ocpInterface->eval_obj(evalRequest->x, evalResult->obj);
    // Eval constraint
    exit = ocpInterface->eval_constr(evalRequest->x, evalResult->c);
    // Return
    return exit;
}

int KNITROSolver::callbackEvalGA (
    KN_context_ptr             kc,
    CB_context_ptr             cb,
    KN_eval_request_ptr const  evalRequest,
    KN_eval_result_ptr  const  evalResult,
    void              * const  ocpInterface_p
) {
    OCPInterface* ocpInterface = CASTOCPINTERFACE(ocpInterface_p);
    if (evalRequest->type != KN_RC_EVALGA) { return -1; }
    int exit;
    // Eval gradient (dense)
    exit = ocpInterface->eval_obj_grad(evalRequest->x, evalResult->objGrad);
    // Eval jacobian (sparse)
    exit = ocpInterface->eval_constr_jac(evalRequest->x, evalResult->jac);
    // Return
    return exit;
}

int KNITROSolver::callbackEvalH(
    KN_context_ptr             kc,
    CB_context_ptr             cb,
    KN_eval_request_ptr const  evalRequest,
    KN_eval_result_ptr  const  evalResult,
    void              * const  ocpInterface_p
) {
    OCPInterface* ocpInterface = CASTOCPINTERFACE(ocpInterface_p);
    if (evalRequest->type != KN_RC_EVALH
        && evalRequest->type != KN_RC_EVALH_NO_F) { return -1; }
    int exit;
    // Eval lagragian hessian
    exit = ocpInterface->eval_hessian(evalRequest->x,  *(evalRequest->sigma), evalRequest->lambda, evalResult->hess);
    // Return
    return exit;
}

int KNITROSolver::callbackNewPoint(
    KN_context_ptr        kc,
    const double * const  x,
    const double * const  lambda,
    void              * const  ocpInterface_p
) {
    OCPInterface* ocpInterface = CASTOCPINTERFACE(ocpInterface_p);
    // Get iteration data
    int iter;
    double obj_value, inf_pr, inf_du;
    KN_get_number_iters(kc, &iter); // get number iters
    KN_get_obj_value(kc, &obj_value); // get objective value
    KN_get_abs_feas_error(kc, &inf_pr); // get feasibility error
    KN_get_abs_opt_error(kc, &inf_du); // get optimality error

    // fix iter
    iter--;

    // avoid repeating iter
    if (iter<=ocpInterface->num_iter) { return 0; } 
    // eval current solution data
    if ((ocpInterface->print_itersol>0) && ((iter % ocpInterface->print_itersol) == 0)) {
        KN_get_var_primal_values_all(kc, ocpInterface->z_opt);
        KN_get_var_dual_values_all(kc, ocpInterface->lamz_opt);
        KN_get_con_dual_values_all(kc, ocpInterface->lamg_opt);
        ocpInterface->eval_constr(ocpInterface->z_opt, ocpInterface->g_opt);
    }
    // print iter
    bool exit = ocpInterface->print_iter(iter, obj_value, inf_pr, inf_du, false); // false to not force print

    // Return
    return (exit) ? 0 : KN_RC_USER_TERMINATION;
}

void KNITROSolver::finalize_solution(
    KN_context_ptr        kc,
    int                   status,
    OCPInterface*         ocpInterface
) {
    // get solution
    KN_get_obj_value(kc, &(ocpInterface->J_opt));
    KN_get_var_primal_values_all(kc, ocpInterface->z_opt);
    KN_get_var_dual_values_all(kc, ocpInterface->lamz_opt);
    KN_get_con_dual_values_all(kc, ocpInterface->lamg_opt);
    // eval nlp functions at solution
    ocpInterface->eval_ocpfuncs();
    // solution data
    KN_get_number_iters(kc, &ocpInterface->num_iter);
    ocpInterface->obj_history[ocpInterface->num_iter] = ocpInterface->J_opt; // add last for consistency with IPOPT
    //ocpInterface->mu_curr = -1; // get function not available for KNITRO
    // set guess solution to optimal if success
    if ( (status == KN_RC_OPTIMAL_OR_SATISFACTORY) || 
        ((status>=KN_RC_NEAR_OPT) && (status<KN_RC_INFEASIBLE)) )
        ocpInterface->set_optsol_as_guess();
    // save cpu time
    KN_get_solve_time_real(kc, &ocpInterface->tcpu_tot);
    ocpInterface->tcpu_eval = 0; // get function not available for KNITRO
}


/** KNITROSolver::callSolve */
int KNITROSolver::callSolve(
    OCPInterface *ocpInterface
){
    /* Declare varibales */
    KN_context* kc; // KNITRO context
    CB_context* cb; // Callback context
    int status; // KNITRO status output

    /* Init KNITRO context */
    status = KN_new(&kc);
    if (status || !kc) {
        /* Print exit status */
        if (ocpInterface->display)
            ocpInterface->print_stats(status, KN_RC_OPTIMAL_OR_SATISFACTORY, KN_RC_NEAR_OPT, KN_RC_USER_TERMINATION, KN_RC_ITER_LIMIT_FEAS, KN_RC_INFEASIBLE);
        return status;
    }

    /* Set KNITRO option for log */
    KN_set_int_param(kc, KN_PARAM_OUTLEV, KN_OUTLEV_ITER);
    KN_set_int_param(kc, KN_PARAM_OUTMODE, KN_OUTMODE_FILE);
    //if (ocpInterface->logfile.empty()) ocpInterface->logfile = "knitro.log";
    if (ocpInterface->logfile == "none") ocpInterface->logfile = ""; // no log file
    KN_set_char_param(kc, KN_PARAM_OUTNAME, ocpInterface->logfile.c_str());
    
    /* Get dims */
    int nz = ocpInterface->nz;
    int ng = ocpInterface->ng;
    int nnzj = ocpInterface->nnzj;
    int nnzh = ocpInterface->nnzh;

    /* Get NLP bounds and init */
    double *z0 = NULL, *lamz0 = NULL, *lamg0 = NULL;
    double *lbz = NULL, *ubz = NULL, *lbg = NULL, *ubg = NULL;
    z0 = new double[1+nz];
    lbz = new double[1+nz];
    ubz = new double[1+nz];
    lbg = new double[1+ng];
    ubg = new double[1+ng];
    if (ocpInterface->check_lambda_guess()) {
        lamz0 = new double[1+nz];
        lamg0 = new double[1+ng];
    }
    ocpInterface->get_nlp_init(z0 != NULL, z0, 
                       lamz0 != NULL, lamz0,
                       lamg0 != NULL, lamg0); 
    ocpInterface->get_nlp_bounds(lbz, ubz, lbg, ubg);

    /* Get sparse patterns */
    int *irj = new int[1+nnzj];
    int *jcj = new int[1+nnzj];
    int *irh = new int[1+nnzh];
    int *jch = new int[1+nnzh];
    ocpInterface->get_pattern_jac(irj, jcj);
    ocpInterface->get_pattern_hess(irh, jch);

    /* KNITRO problem definition */
    KN_add_vars(kc, nz, NULL);
    KN_add_cons(kc, ng, NULL);
    if (z0) KN_set_var_primal_init_values_all(kc, z0);
    if (lamz0) KN_set_var_dual_init_values_all(kc, lamz0);
    if (lamg0) KN_set_con_dual_init_values_all(kc, lamg0);
    KN_set_var_lobnds_all(kc, lbz);
    KN_set_var_upbnds_all(kc, ubz);
    KN_set_con_lobnds_all(kc, lbg);
    KN_set_con_upbnds_all(kc, ubg);

    /* KNITRO callbacks */
    KN_add_eval_callback_all(kc, KNITROSolver::callbackEvalFC, &cb);
    KN_set_cb_grad(kc, cb, KN_DENSE, NULL, nnzj, irj, jcj, KNITROSolver::callbackEvalGA);
    if (!ocpInterface->flag_hessian) 
        KN_set_cb_hess(kc, cb, nnzh, irh, jch, KNITROSolver::callbackEvalH);
    KN_set_newpt_callback(kc, KNITROSolver::callbackNewPoint, (void*) ocpInterface);
    KN_set_cb_user_params(kc, cb, (void*) ocpInterface);
    
    /* KNITRO options */
    KN_set_int_param(kc, KN_PARAM_ALG, KN_ALG_BAR_DIRECT); // interior/direct algorithm
    KN_set_obj_goal(kc, KN_OBJGOAL_MINIMIZE); // minimize objective
    if (ocpInterface->flag_hessian) {
        KN_set_int_param(kc, KN_PARAM_HESSOPT, KN_HESSOPT_LBFGS); // limited-memory BFGS Hessian
        KN_set_int_param(kc, KN_PARAM_LMSIZE, 5); // default 10, decrease a bit for faster hessian approx
    }
    else    
        KN_set_int_param(kc, KN_PARAM_HESSOPT, KN_HESSOPT_EXACT); // exact hessian
    if (ocpInterface->check_lambda_guess())
        KN_set_int_param(kc, KN_PARAM_STRAT_WARM_START, KN_STRAT_WARM_START_YES); // warm start
    if (ocpInterface->max_iter<0) 
        ocpInterface->max_iter = MAX_DEFAULT_ITER; // use default if <0
    KN_set_int_param(kc, KN_PARAM_MAXIT, ocpInterface->max_iter); // max iter
    if (ocpInterface->mu_init>0)
        KN_set_double_param(kc, KN_PARAM_BAR_INITMU, ocpInterface->mu_init); // initial barrier parameter

    /* Load KNITRO option file if any */
    const char optfile[] = "knitro.opt";
    if (ocpInterface->check_nlpopt_file(optfile))
        KN_load_param_file(kc, optfile); //load opt file
    
    /* Init history */
    ocpInterface->init_hist_mem(ocpInterface->max_iter);

    /* KNITRO solver */
    status = KN_solve(kc);

    /* Finalize solution */
    KNITROSolver::finalize_solution(kc, status, ocpInterface);

    /* Print last iter */
    int iter;
    double obj_value, pr, du;
    KN_get_number_iters(kc, &iter); // get number iters
    KN_get_obj_value(kc, &obj_value); // get objective value
    KN_get_abs_feas_error(kc, &pr); // get feasibility error
    KN_get_abs_opt_error(kc, &du); // get optimality error    
    iter--; // fix iter
    if (ocpInterface->display)
        ocpInterface->print_iter(ocpInterface->num_iter, ocpInterface->J_opt, pr, du, true);

    /* Print exit status */
    if (ocpInterface->display)
        ocpInterface->print_stats(status, KN_RC_OPTIMAL_OR_SATISFACTORY, KN_RC_NEAR_OPT, KN_RC_USER_TERMINATION, KN_RC_ITER_LIMIT_FEAS, KN_RC_INFEASIBLE);

    /* Free KNITRO mem */
    KN_free(&kc); 

    /* Free mem */
    delete[] z0;  delete[] lbz; delete[] ubz; 
    delete[] lbg; delete[] ubg; 
    delete[] irj; delete[] jcj;
    delete[] irh; delete[] jch;

    /* Return */
    return status;
}

/** A wrapper around KNITROSolver::callSolve to load the static function from the dynamic library in runtime */
#ifndef MINOS_DECLSPEC
#define MINOS_DECLSPEC __declspec(dllexport)
#endif
#ifdef __cplusplus
extern "C" {
#endif
MINOS_DECLSPEC int callSolve(
    OCPInterface* ocpInterface
) {
    return KNITROSolver::callSolve(ocpInterface);
}
#ifdef __cplusplus
} /* extern "C" */
#endif
