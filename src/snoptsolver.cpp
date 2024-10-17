/*  File snoptsolver.cpp
    NLP solver interface for SNOPT
    Copyright (C) 2024 Stefano Lovato
*/

#include "snoptsolver.h" // base class with all method implementations
#include "macros.h" // for nlp macros
#include <time.h> // for timing

void SNOPTSolver::nlpfun(
    int *Status, int *n,   double x[],
    int *needF,  int *nF,  double F[],
    int *needG,  int *neG, double G[],
    char   cu[], int   *lencu,
    int    iu[], int   *leniu,
    double ru[], int   *lenru
) {
    clock_t start_time = clock(); // start time
    OCPInterface* ocpInterface = CASTOCPINTERFACE(iu);
    if (*needF > 0) { // calc obj and constr
        double* obj = &F[0]; // obj
        double* g  = &F[1]; // constr
        // eval obj
        ocpInterface->eval_obj(x, obj);
        // eval constr
        ocpInterface->eval_constr(x, g);
    }
    if (*needG > 0) { // calc obj grad and constr jac
        double* grad = &G[0]; // obj grad
        double* jac = &G[*n]; // constr jacs
        // eval obj grad
        ocpInterface->eval_obj_grad(x, grad);
        // eval constr jac
        ocpInterface->eval_constr_jac(x, jac);
    }
    // update evaluation time
    ocpInterface->tcpu_eval += ((double) clock() - start_time) / CLOCKS_PER_SEC;
}

void SNOPTSolver::iter_callback(
    int *iAbort, int KTcond[], int *MjrPrt, int *minimz,
    int *m, int *maxS, int *n, int *nb,
    int *nnCon0, int *nnCon, int *nnObj0, int *nnObj, int *nS,
    int *itn, int *nMajor, int *nMinor, int *nSwap,
    double *condHz, int *iObj, double *sclObj, double *ObjAdd,
    double *fObj, double *fMrt, double PenNrm[], double *step,
    double *prInf, double *duInf, double *vimax, double *virel, int hs[],
    int *ne, int *nlocJ, int locJ[], int indJ[], double Jcol[], int *negCon,
    double Ascale[], double bl[], double bu[],
    double Fx[], double fCon[], double gCon[], double gObj[],
    double yCon[], double pi[], double rc[], double rg[], double x[],
    char cu[], int *lencu, int iu[], int *leniu, double ru[], int *lenru,
    char cw[], int *lencw, int iw[], int *leniw, double rw[], int *lenrw 
) {
    OCPInterface* ocpInterface = CASTOCPINTERFACE(iu);
    if (*nMajor<=ocpInterface->num_iter) { return; } // avoid repeating iter
    // eval current solution data
    if ((ocpInterface->print_itersol>0) && ((*nMajor % ocpInterface->print_itersol) == 0)) {
        double* lamz = ru + 0; // lamz are first nz entries of ru
        double* lamg = ru + (*n+1); // lamg are last ng entries of ru
        memcpy(ocpInterface->z_opt, x, ocpInterface->nz*sizeof(double));
        // below seems not working yet
        memcpy(ocpInterface->lamz_opt, lamz, ocpInterface->nz*sizeof(double));
        memcpy(ocpInterface->lamg_opt, lamg, ocpInterface->ng*sizeof(double));
        ocpInterface->eval_constr(ocpInterface->z_opt, ocpInterface->g_opt);
    }
    // print iter
    bool exit = ocpInterface->print_iter(*nMajor, *fObj, *prInf, *duInf, false);
    *iAbort = (exit==false);
}

void SNOPTSolver::finalize_solution(
    int status,
    double* z,
    double* F,
    double* lamz,
    double* lamF,
    OCPInterface* ocpInterface
) {
    memcpy(ocpInterface->z_opt,     z,      ocpInterface->nz*sizeof(double)); // nlp vars
    memcpy(ocpInterface->lamz_opt,  lamz,   ocpInterface->nz*sizeof(double)); // nlp mult
    memcpy(ocpInterface->lamg_opt,  lamF+1, ocpInterface->ng*sizeof(double)); // nlp constr mult
    // eval nlp functions at solution
    ocpInterface->eval_ocpfuncs();
    // solution data
    // num_iter handled separately
    // no need to save mu_curr b/c SQP
    // set guess solution to optimal if success
    if (status < 10) // <10 for success (SNOPT manual Sec. 3.4, page 21)
        ocpInterface->set_optsol_as_guess();
    // save cpu time
    // tcpu_eval and tcpu_tot already saved
}

int SNOPTSolver::callSolve(
    OCPInterface *ocpInterface
) {

    /* Init CPU times and num iter */
    ocpInterface->tcpu_tot = 0;
    ocpInterface->tcpu_eval = 0;
    ocpInterface->num_iter = -1;

    /* Set log file */
    int iprint = 0;
    if (ocpInterface->logfile != "none") {
        iprint = 1; // enable log
        //if (ocpInterface->logfile.empty()) ocpInterface->logfile = "snopt.log";
    } // else // no log file

    /* Create an instance of snProblem and init */
    snProblem nlp;
    snInitX(&nlp, (char*) "snoptnlp", (char*) ocpInterface->logfile.c_str(), iprint, (char*) "", 0);

    /* Get dims */
    int nz = ocpInterface->nz;
    int ng = ocpInterface->ng;
    int nnzj = ocpInterface->nnzj;
    //int nnzh = ocpInterface->nnzh;

    /* Get NLP bounds and init 
    The NLP objective and constraint fuctions are collected in a unique function F,
    which has 1+ng components. Here F[0] is the objective, and F[1..end] are the
    constraints. We need to define guess, bounds and lagrange multipliers accordingly.
    */
    double *z = NULL, *lamz = NULL, *lamF = NULL;
    double *lbz = NULL, *ubz = NULL, *lbF = NULL, *ubF = NULL;
    double *F = NULL; // F[0] is obj, F[1..end] is constr
    double *lamzF = new double[1+nz+ng] { 0 }; // collect lam for z , obj, and constr in same array
    z = new double[1+nz];
    lbz = new double[1+nz];
    ubz = new double[1+nz];
    lamz = lamzF + 0; // lamz starts at lamzF[0]
    F = new double[1+ng] { 0 }; // for both obj and constr
    lbF = new double[1+ng]; // for both obj and constr
    ubF = new double[1+ng]; // for both obj and constr
    lamF = lamzF + nz; // lamF starts at lamzF[nz], for both obj and constr
    ocpInterface->get_nlp_init(z != NULL, z, false, NULL, false, NULL); // init only nlp vars
    lbF[0] = -1.0e20; ubF[0] = +1.0e20; // no bounds for obj
    ocpInterface->get_nlp_bounds(lbz, ubz, lbF+1, ubF+1); // first val is obj
    if (ocpInterface->check_lambda_guess()) {
        // init multiplites
        ocpInterface->get_nlp_init(false, NULL, true, lamz, true, lamF+1); // first val is obj
        lamF[0] = 0; // mult for obj, should be 0 b/c no bounded
    }

    /* Create zstate and Fstate */
    int *zstate =  new int[1+nz] { 0 }; // always 0 (?)
    int *Fstate = new int[1+ng] { 0 };  // always 0 (?)

    /* Get sparse patterns 
    The NLP objective gradient and constraint jacobian are collected in a unique function G,
    which has nz (for the full gradient) + nnzj (for the sparse jacobian) non-zero entries.
    G[0..nz-1] is the objective gradient, which is full and has row index 0 and column
    indexes 0, 1, ..., nz-1.
    G[nz...end] is the constraint jacobian, which is sparse and has row and column indexes
    obtained from get_pattern_jac, with row indexes increases by 1 to account for objective
    in the first row.
    */
    int nnzG = nz+nnzj; 
    int *irG = new int[1+nnzG]; 
    int *jcG = new int[1+nnzG];
    // obj gradient (first row), which is full (element from 0 to nz-1)
    for (int i = 0; i < nz; ++i) {
        irG[i] = 0; // ir: 0, 0, 0, ..., 0
        jcG[i] = i; // jc: 0, 1, 2, ..., nz-1
    }
    // constr jacobian, which is sparse and starts from element nz
    ocpInterface->get_pattern_jac(irG+nz, jcG+nz);
    // increase irG by 1, b/c F[0] is obj
    for (int i = nz; i < nz+nnzj; ++i) ++irG[i];
    
    /* Allocates and set user workspace */
    setUserI(&nlp, (int*) ocpInterface, 1); // save ocp object as (int*)
    setUserR(&nlp, lamzF, 1); // save lamzF to have access in iter_callback


    /* SNOPT options */
    if (ocpInterface->max_iter<0) 
        ocpInterface->max_iter = MAX_DEFAULT_ITER; // use default if <0
    setIntParameter(&nlp, (char*) "Major iterations limit", ocpInterface->max_iter);
    setParameter(&nlp, (char*) "Solution No"); // do not print sol into snopt.log

    int start = 0; // 0 for cold start
    if (ocpInterface->check_lambda_guess()) start = 2; // 2 for warm start

    /* Load SNOPT option file if any */
    char optfile[] = "snopt.opt";
    if (ocpInterface->check_nlpopt_file(optfile))
        setSpecsfile(&nlp, optfile);

    /* Init history */
    ocpInterface->init_hist_mem(ocpInterface->max_iter);

    /* Set iter callback function
    Here either setSTOP or setLog may be employed to print info at every
    major iter. However, using setLog results in turning off the iter printing 
    into the log file (snopt.log), thus setSTOP is preferred.
    */
    setSTOP(&nlp, SNOPTSolver::iter_callback);
    
    /* Solve the problem */
    int    nS, nInf;
    double sInf;
    clock_t start_time = clock(); // start time
    int status = solveA( // solveA for general NLP with no assumption
            &nlp, start, // start [0/2] for cold/warm start
            1+ng, nz, // F has 1+ng comps (F[0] is obj), nlp has nz vars
            0, 0, // objAdd=0 (not neccesary) and objRow=0, b/c objective is F[0]
            SNOPTSolver::nlpfun, // NLP functin eval
            0, NULL, NULL, NULL, // generally no constant terms in the jacobian
            nnzG, irG, jcG, // gradient+jacobian pattern
            lbz, ubz, lbF, ubF, // z and F bounds
            z, zstate, lamz, // guess on z and lamz
            F, Fstate, lamF, // guess on F and lamF
            &nS, &nInf, &sInf);
    ocpInterface->tcpu_tot = ((double) (clock() - start_time)) / CLOCKS_PER_SEC;

    /* Finalize the solution */
    SNOPTSolver::finalize_solution(status, z, F, lamz, lamF, ocpInterface);

    /* Print last iter */
    int num_iter = ocpInterface->num_iter;
    double obj = F[0];
    double inf_pr = ocpInterface->infpr_history[num_iter];
    double inf_du = ocpInterface->infdu_history[num_iter];
    ocpInterface->print_iter(num_iter, obj, inf_pr, inf_du, true);

    /* Print exit status */
    if (ocpInterface->display)
        ocpInterface->print_stats(status, 1, 3, 74, 32, 13); // exit code from SNOPT manual Sec. 3.4, page 21

    /* Free mem */
    deleteSNOPT(&nlp);
    delete[] lamzF; 
    delete[] z; delete[] lbz; delete[] ubz; 
    delete[] F; delete[] lbF; delete[] ubF; 
    delete[] zstate; delete[] Fstate; 
    delete[] irG; delete[] jcG;

    return status;
}

/** A wrapper around SNOPTSolver::callSolve to load the static function from the dynamic library in runtime */
#ifdef __cplusplus
extern "C" {
#endif
MINOS_EXPORT_API int callSolve(
    OCPInterface* ocpInterface
) {
    return SNOPTSolver::callSolve(ocpInterface);
}
#ifdef __cplusplus
} /* extern "C" */
#endif
