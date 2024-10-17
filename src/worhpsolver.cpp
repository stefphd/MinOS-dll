/*  File worhpsolver.cpp 
    NLP solver interface for WORHP
    Copyright (C) 2024 Stefano Lovato
*/

#include "worhpsolver.h" // base class with all method implementations
#include "macros.h" // for nlp macros

/* Forward declarations */
std::unique_ptr<std::ofstream> WORHPSolver::logfile;

/**
* WORHP solver class implementation
*/

/** Constructor */
WORHPSolver::WORHPSolver(
    OCPInterface* ocpInterface
) {
    /* store the ocpInterface pointer */
    this->ocpInterface = ocpInterface;

    /* Redirect WORHP print to a log file */
    WORHPSolver::logfile = std::unique_ptr<std::ofstream>(new std::ofstream());
    if (ocpInterface->logfile != "none" )
        WORHPSolver::logfile->open(ocpInterface->logfile);
    SetWorhpPrint(WORHPSolver::print_logfile);

    /*
    * Properly zeros everything, or else the following routines
    * could get confused.
    */
    WorhpPreInit(&opt, &wsp, &par, &cnt);

    /* Init param */
    InitParams(&cnt.status, &par);

    /* Initialize dimensions */
    init_dim();

    /* Allocate memory */
    alloc_mem();

    /* Data structure initialisation. */
    WorhpInit(&opt, &wsp, &par, &cnt);
    if (cnt.status != FirstCall) {
        status = cnt.status;
        return;
    }

    /* All functions are assumed nonlinear by default */

    /* Get NLP bounds and init */
    ocpInterface->get_nlp_bounds(opt.XL, opt.XU, opt.GL, opt.GU); // bounds
    ocpInterface->get_nlp_init(true, opt.X, false, NULL, false, NULL); // init x
    if (ocpInterface->check_lambda_guess()) // init lambda and mu
        ocpInterface->get_nlp_init(false, NULL, true, opt.Lambda, true, opt.Mu);

    /* Get jacobian and hessian sparsity */
    set_sparsity();
}

/** Destructor */
WORHPSolver::~WORHPSolver(){   
    WorhpFree(&opt, &wsp, &par, &cnt);
    delete[] mapj; delete[] jac; 
    delete[] maph; delete[] hess;
}

/** Initialize dimensions */
void WORHPSolver::init_dim(){  
    /* Specify number of variables and constraints. */
    opt.n = ocpInterface->nz;
    opt.m = ocpInterface->ng;
    /* Specify nonzeros of derivative matrices. */
    wsp.DF.nnz = ocpInterface->nz;        /* dense, structure init by solver */
    wsp.DG.nnz = ocpInterface->nnzj;      /* sparse, structure init by solver */
    
    /* sparse, structure init by user 
    WORHP wants Hessian with sparse strickly lower triangular part 
    and full diagonal part: we need to increase nnzh by zero diagonal entries.
    */
    wsp.HM.nnz = 0; // default
    if (!ocpInterface->flag_hessian) {
        int c = ocpInterface->nz; // number of all diagonal terms
        // get pattern to count number of strickly lower triangular entries
        int *irh = new int[1+ocpInterface->nnzh];
        int *jch = new int[1+ocpInterface->nnzh];
        ocpInterface->get_pattern_hess(irh, jch);
        // count number of off diagonal term - hessian is already lower triangular here
        for (int i = 0; i < ocpInterface->nnzh; ++i)
            if (jch[i] != irh[i]) ++c; // increment for off-diagonal term
        wsp.HM.nnz = c; // set nnz of hessian (>=nnzh)
        // free mem
        delete[] irh; delete[] jch;
    }
}

/** Initialize dimensions */
void WORHPSolver::alloc_mem(){  
    mapj = new int[1+wsp.DG.nnz] { -1 }; 
    jac = new double[1+ocpInterface->nnzj];
    maph = new int[1+wsp.HM.nnz] { -1 };
    hess = new double[1+ocpInterface->nnzh];
}

/* Set sparsity */
void WORHPSolver::set_sparsity() {
    /* WORHP by default employs coordinate storage format (CS) with major column: 
    we need to sort arrays with increasing column indexes and save the mapping
    from unsorted to sorted indexes.
    */
    // jacobian
    if (wsp.DG.NeedStructure) {
        // get jacobian pattern
        int *irj = new int[1+ocpInterface->nnzj];
        int *jcj = new int[1+ocpInterface->nnzj];
        ocpInterface->get_pattern_jac(irj, jcj);
        // find mapping from unsourted to sorted
        OCPInterface::get_sorted_mapping(jcj, mapj, ocpInterface->nnzj);
        // use FORTRAN indexing, CS format
        int ir, jc; // tmp row and column indexes
        for (int i = 0; i < ocpInterface->nnzj; ++i) {
            ir = irj[mapj[i]]; jc = jcj[mapj[i]];
            wsp.DG.row[i] = ir + 1;
            wsp.DG.col[i] = jc + 1;
        }
        delete[] irj; delete[] jcj;
    }
    // hessian
    if (wsp.HM.nnz>0 && wsp.HM.NeedStructure) {
        // get hessian pattern
        int *irh = new int[1+ocpInterface->nnzh];
        int *jch = new int[1+ocpInterface->nnzh];
        ocpInterface->get_pattern_hess(irh, jch);
        // find mapping from unsourted to sorted
        int *tmpmap = new int[1+ocpInterface->nnzh]; // store temp mapping from unsorted to sorted hessian
        OCPInterface::get_sorted_mapping(jch, tmpmap, ocpInterface->nnzh);
        // order first strickly lower triangular part, and then all diagonal terms
        int c = 0; // counter of strickly lower triangular term
        int ir, jc; // temp indexes
        // use FORTRAN indexing, CS format
        for (int i = 0; i < ocpInterface->nnzh; ++i) {
            ir = irh[tmpmap[i]]; jc = jch[tmpmap[i]];
            if (ir > jc) { //strickly lower part
                wsp.HM.row[c] = ir+1; 
                wsp.HM.col[c] = jc+1; 
                maph[c] = tmpmap[i]; // mapping from unsort to sort
                ++c; // increment counter
            } else if (ir == jc) // diagonal term
                maph[ir+wsp.HM.nnz-ocpInterface->nz] = tmpmap[i];
        }
        // add diagonal indexes, use FORTRAN indexing
        for (int i = 0; i < ocpInterface->nz; ++i) {
            wsp.HM.row[i+wsp.HM.nnz-ocpInterface->nz] = i+1;
            wsp.HM.col[i+wsp.HM.nnz-ocpInterface->nz] = i+1;
        }
        // free
        delete[] irh; delete[] jch;
        delete[] tmpmap;
    }
}

/* Intermediate callback */
bool WORHPSolver::int_callback() {
    // avoid repeating iter
    if (wsp.MajorIter<=ocpInterface->num_iter) { return true; } 
    // eval current solution data
    if ((ocpInterface->print_itersol>0) && ((wsp.MajorIter % ocpInterface->print_itersol) == 0)) {
        memcpy(ocpInterface->z_opt, opt.X, ocpInterface->nz*sizeof(double));
        memcpy(ocpInterface->lamz_opt, opt.Lambda, ocpInterface->nz*sizeof(double));
        memcpy(ocpInterface->lamg_opt, opt.Mu, ocpInterface->ng*sizeof(double));
        ocpInterface->eval_constr(ocpInterface->z_opt, ocpInterface->g_opt);
    }
    // get inf pr and inf du
    switch (par.Algorithm) {
        case 1: // WORHP-SQP
            //ocpInterface->mu_curr not used for SQP
            inf_pr = wsp.NormMax_CV;
            inf_du = wsp.OptiMax;
            break;
        case 2: // WORHP-IP
            inf_pr = wsp.IP_FeasMax;
            inf_du = wsp.IP_ComplMax;
            break;
    }
    // print iter
    return ocpInterface->print_iter(wsp.MajorIter, opt.F/wsp.ScaleObj, inf_pr, inf_du, false);
}

/* Solve */
int WORHPSolver::solve() {
    bool userRequestStop = false; // flag for stop requested by user
    // Reset timers
    ResetTimer(&totTimer);
    ResetTimer(&evalTimer);
    // solving loop
    StartTimer(&totTimer);
    while(!userRequestStop && cnt.status < TerminateSuccess && cnt.status > TerminateError) {
        /* Solve action */
        if (GetUserAction(&cnt, callWorhp))
            Worhp(&opt, &wsp, &par, &cnt);
            /* No DoneUserAction! */
        /* Evaluation action */
        StartTimer(&evalTimer);
        if (GetUserAction(&cnt, iterOutput)) {
            IterationOutput(&opt, &wsp, &par, &cnt);
            if (!int_callback()) userRequestStop = true; // user requested stop
            DoneUserAction(&cnt, iterOutput);
        }
        if (GetUserAction(&cnt, evalF)) {
            ocpInterface->eval_obj(opt.X, &opt.F);
            opt.F *= wsp.ScaleObj;
            DoneUserAction(&cnt, evalF);
        }
        if (GetUserAction(&cnt, evalG)) {
            ocpInterface->eval_constr(opt.X, opt.G);
            DoneUserAction(&cnt, evalG);
        }
        if (GetUserAction(&cnt, evalDF)) {
            ocpInterface->eval_obj_grad(opt.X, wsp.DF.val);
            for (int i = 0; i < wsp.DF.nnz; ++i)
                wsp.DF.val[i] *= wsp.ScaleObj;
            DoneUserAction(&cnt, evalDF);
        }
        if (GetUserAction(&cnt, evalDG)) {
            ocpInterface->eval_constr_jac(opt.X, jac);
            for (int i = 0; i < wsp.DG.nnz; ++i)
                wsp.DG.val[i] = jac[mapj[i]];
            DoneUserAction(&cnt, evalDG);
        }
        if (GetUserAction(&cnt, evalHM)) {
            ocpInterface->eval_hessian(opt.X, wsp.ScaleObj, opt.Mu, hess);
            for (int i = 0; i < wsp.HM.nnz; ++i) {
                int idx = maph[i];
                if (idx<0) wsp.HM.val[i] = 0.0; // diagonal zero entry
                else wsp.HM.val[i] = hess[maph[i]]; 
            }
            DoneUserAction(&cnt, evalHM);
        }
        StopTimer(&evalTimer);
        /* WorhpFidif action */
        if (GetUserAction(&cnt, fidif))
            WorhpFidif(&opt, &wsp, &par, &cnt);
            /* No DoneUserAction! */
        
    }
    StopTimer(&totTimer);
    StopTimer(&evalTimer);

    /* Final task */
    StatusMsg(&opt, &wsp, &par, &cnt);
    WORHPSolver::logfile->close();

    /* Return */
    if (userRequestStop) { return -1; }
    return cnt.status;
}

/* Finalize the solution */
void WORHPSolver::finalize_solution() {
    // copy solution
    memcpy(ocpInterface->z_opt, opt.X, ocpInterface->nz*sizeof(double));
    memcpy(ocpInterface->lamz_opt, opt.Lambda, ocpInterface->nz*sizeof(double));
    memcpy(ocpInterface->lamg_opt, opt.Mu, ocpInterface->ng*sizeof(double));
    // eval nlp functions at solution
    ocpInterface->eval_ocpfuncs();
    // solution data
    ocpInterface->num_iter = wsp.MajorIter;
    if (par.Algorithm==2) ocpInterface->mu_curr = wsp.IP_Barrier; // only for WORHP-IP
    // set guess solution to optimal if success
    if (cnt.status >= TerminateSuccess)
        ocpInterface->set_optsol_as_guess();
    // save cpu time
    ocpInterface->tcpu_tot = GetTimer(&totTimer);
    ocpInterface->tcpu_eval = GetTimer(&evalTimer);
}

void WORHPSolver::print_logfile(
    int mode, 
    const char msg[]
) {
    if (logfile->is_open())
        (*logfile) << msg << std::endl;
}

/** WORHPSolver::callSolve */
int WORHPSolver::callSolve(
    OCPInterface *ocpInterface
){
    /* Create WORHP solver object */
    WORHPSolver solver(ocpInterface);
    int exitstatus = solver.cnt.status;
    if (exitstatus<0) { // error init
        /* Print exit status */
        if (ocpInterface->display)
            ocpInterface->print_stats(exitstatus, OptimalSolution, AcceptableSolution, -1, MaxIter, LocalInfeas);
        return exitstatus;
    }

    /* Options */
    solver.par.Algorithm = 1;
    if (ocpInterface->max_iter<0) 
        ocpInterface->max_iter = MAX_DEFAULT_ITER; // use default if <0
    solver.par.MaxIter = ocpInterface->max_iter;
    if (ocpInterface->flag_hessian) {   
        solver.par.UserHM = false;
        solver.par.BFGSmethod = 2;
    } else 
        solver.par.UserHM = true;
    if ((ocpInterface->mu_init>0) && solver.par.Algorithm==2) // set mu_init if >0 and for WORHP-IP
        solver.par.IP_BarrierInit = ocpInterface->mu_init;

    /* Load WORHP option file if any */
    const char optfile[] = "worhp.xml";
    if (ocpInterface->check_nlpopt_file(optfile))
        ReadParamsNoInit(&solver.status, "worhp.xml", &solver.par);

    /* Init history */
    ocpInterface->init_hist_mem(ocpInterface->max_iter);

    /* Solve */
    exitstatus = solver.solve();

    /* Finalize solution */
    solver.finalize_solution();

    /* Print last iter */
    ocpInterface->print_iter(ocpInterface->num_iter,  solver.opt.F/solver.wsp.ScaleObj, solver.inf_pr, solver.inf_du, true);

    /* Print exit status */
    if (ocpInterface->display)
        ocpInterface->print_stats(exitstatus, OptimalSolution, AcceptableSolution, -1, MaxIter, LocalInfeas);

    /* Return */
    return exitstatus;
}

/** A wrapper around WORHPSolver::callSolve to load the static function from the dynamic library in runtime */
#ifdef __cplusplus
extern "C" {
#endif
MINOS_EXPORT_API int callSolve(
    OCPInterface* ocpInterface
) {
    return WORHPSolver::callSolve(ocpInterface);
}
#ifdef __cplusplus
} /* extern "C" */
#endif
