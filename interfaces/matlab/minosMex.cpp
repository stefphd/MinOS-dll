/*  File minosMex.cpp 
    Mex interface for MinOS.
    Copyright (C) 2024 Stefano Lovato
*/

/* Macros */
#define ASSERTIN(x)             if (x) mexErrMsgIdAndTxt("minosMex:wrongNumInputs", "Expecting one input argument for 'minosMex'.")
#define ASSERTOUT(x)            if (x) mexErrMsgIdAndTxt("minosMex:wrongNumOutputs", "Expecting one output argument for 'minosMex'.")
#define ASSERTSTRUCT(x,p)       if ( (x) && (!mxIsStruct(x))) { if (p) delete p; mexErrMsgIdAndTxt("minosMex:wrongInputType", "Expecting structure for '%s'.", #x); }
#define ASSERTPOSINTSCALAR(x,p) if ( (x) && (!mxIsDouble(x) || mxIsComplex(x) || mxGetNumberOfElements(x)!=1 || mxGetScalar(x)<=0 || \
                                    ((((int) mxGetScalar(x)) != mxGetScalar(x)) != 0) ) ) \
                                    { if (p) delete p; mexErrMsgIdAndTxt("minosMex:notScalar","Expecting positive integer scalar for '%s'.", #x); }
#define ASSERTNONEGINTSCALAR(x,p) if ( (x) && (!mxIsDouble(x) || mxIsComplex(x) || mxGetNumberOfElements(x)!=1 || mxGetScalar(x)<0 || \
                                    ((((int) mxGetScalar(x)) != mxGetScalar(x)) != 0) ) ) \
                                    {if (p) delete p; mexErrMsgIdAndTxt("minosMex:notScalar","Expecting non-negative integer scalar for '%s'.", #x); }
#define ASSERTPOSSCALAR(x,p)    if ( (x) && (!mxIsDouble(x) || mxIsComplex(x) || mxGetNumberOfElements(x)!=1 || mxGetScalar(x)<=0 ) ) \
                                    {if (p) delete p; mexErrMsgIdAndTxt("minosMex:notScalar","Expecting positive scalar for '%s'.", #x); }
#define ASSERTSCALAR(x,p)       if ( (x) && ( !mxIsDouble(x) || mxIsComplex(x) || mxGetNumberOfElements(x)!=1 ) ) \
                                    { if (p) delete p; mexErrMsgIdAndTxt("minosMex:notScalar","Expecting scalar for '%s'.", #x); }
#define ASSERTDOUBLE(x,p)       if  ( (x) && (!mxIsDouble(x) || mxIsComplex(x)) ) \
                                    { if (p) delete p; mexErrMsgIdAndTxt("minosMex:notDouble","Expecting double for '%s'.", #x); }
#define ASSERTLOGICAL(x,p)      if ( (x) && (!mxIsLogicalScalar(x)) ) \
                                    { if (p) delete p; mexErrMsgIdAndTxt("minosMex:notLogical","Expecting logical scalar for '%s'.", #x); }
#define ASSERTSIZE(x,n,m,p)     if ( (x) && ((mxGetM(x) != n) || (mxGetN(x) != m)) ) \
                                    { if (p) delete p; mexErrMsgIdAndTxt("minosMex:wrongSize","Expecting dimension %d-by-%d for '%s' (found %d-by-%d).",n,m,#x,mxGetM(x),mxGetN(x)); }
#define ASSERTNUMEL(x,n,p)      if ( (x) && (mxGetNumberOfElements(x) != n) ) \
                                    { if (p) delete p; mexErrMsgIdAndTxt("minosMex:wrongSize","Expecting %d elements for '%s' (found %d).",n,#x,mxGetNumberOfElements(x)); }
#define ASSERTSTRING(x,n,p)     if ( (x) && (!mxIsChar(x) || mxGetM(x) < n) ) \
                                    { if (p) delete p; mexErrMsgIdAndTxt("minosMex:notString","Expecting string for '%s'.", #x); }
#define ASSERTOPTNUMEL(x,n,p)   if ( !(x) && (n > 0) ) \
                                    { if (p) delete p; mexErrMsgIdAndTxt("minosMex:missingField", "Expecting field '%s' with %d.", #x, n); } \
                                else if ( (x) && (mxGetNumberOfElements(x) != n) ) \
                                    { if (p) delete p; mexErrMsgIdAndTxt("minosMex:wrongSize","Expecting %d elements for '%s' (found %d).",n,#x,mxGetNumberOfElements(x)); }
#define ASSERT(x,p)             if ( !(x) ) { if (p) delete p; mexErrMsgIdAndTxt("minosMex:falseReturn", "Error: '%s' return false.", #x); }

#define GETFIELD(par,name,opt,p)    mxArray* name = (par) ? mxGetField(par,0,#name) : NULL; \
                                        if ( (name == NULL) & !opt ) {if (p) delete p; mexErrMsgIdAndTxt("minosMex:missingField", "Expecting field '%s' in '%s'.", #name, #par); }
#define GETVAL(x)   ((x) ? mxGetScalar(x)   : 0)
#define GETPTR(x)   ((x) ? mxGetPr(x)       : NULL)
#define GETIR(x)    ((x) ? mxGetIr(x)       : NULL)
#define GETJC(x)    ((x) ? mxGetJc(x)       : NULL)

#define CREATEDOUBLE(x, n, m) mxArray* x; { \
                                mwSize x ## _dims[2] = {static_cast<size_t>(n), static_cast<size_t>(m)}; \
                                x = mxCreateNumericArray(2, x ## _dims, mxDOUBLE_CLASS, mxREAL); }
#define CREATEINT32(x, n, m)  mxArray* x; { \
                                mwSize x ## _dims[2] = {static_cast<size_t>(n), static_cast<size_t>(m)}; \
                                x = mxCreateNumericArray(2, x ## _dims, mxINT32_CLASS, mxREAL); }
#define CREATESTRUCT(x,n,m,k,fields) mxArray* x; { \
                                        mwSize x ## _dims[2] = {static_cast<size_t>(n), static_cast<size_t>(m)}; \
                                        x = mxCreateStructArray(2, x ## _dims, k, fields); }
#define CREATESPARSE(x, n, m, nnz) mxArray* x = mxCreateSparse(static_cast<size_t>(n), static_cast<size_t>(m), static_cast<size_t>(nnz), mxREAL)

#define SETFIELDS(p, fields, vars)  for (int i = 0; i < sizeof(fields)/sizeof(void*); ++i) \
                                        if (p) mxSetField(p, 0, fields[i], vars[i]);
#define ASSERTSUM(x,s,e,p)  if (x) { \
                                double sum = 0; \
                                double* xp = GETPTR(x); \
                                for (int i = 0; i < mxGetNumberOfElements(x); ++i) sum += xp[i]; \
                                if (((sum-s) > e) || ((sum-s) < -e)) \
                                { if (p) delete p; mexErrMsgIdAndTxt("minosMex:invalidField", "Expecting sum to %f for '%s' (sum to %f).\n", s, #x, sum); } \
                            }

#define OUTPUT(x,c) plhs[c++] = x;

/* Includes */
#include "mex.h" // C-MEX file
#include "minos.h" // MinOS interface
#include <fstream> // for write output file
#include <stdarg.h> // for va_list

/* Externs */
#ifdef __cplusplus
    /* These functions are to allow interrupt from MATLAB using CTRL+C (see 
    https://undocumentedmatlab.com/articles/mex-ctrl-c-interrupt) 
    You need to link libut to employ these functions.
    */
    extern "C" bool utIsInterruptPending();
    extern "C" bool utSetInterruptPending(bool);
#else
    extern bool utIsInterruptPending();
    extern bool utSetInterruptPending(bool);
#endif

/* Function to print using mexPrintf with workaround to print immediately */
int print(const char *fmt, ...) {
    // handle variable args
    va_list ap, apcp;
    va_start(ap, fmt); va_copy(apcp, ap); 
    int len = vsnprintf(NULL, 0, fmt, apcp);
    va_end(apcp);
    if (len == -1) { return -1; }
    char *str = (char*) malloc((size_t) len + 1);
    if (!str) { return -1; }
    len = vsprintf(str, fmt, ap);
    va_end(ap);
    if (len == -1) {
        free(str);
        return -1;
    }
    // call mexPrintf and pause(0) to print immediately
    len = mexPrintf("%s", str); mexEvalString("pause(0);");
    free(str);
    return len;
}

/* Function to handle MATLAB CTRL+C */
bool handle_interrupt() {
    /* Handle interrupt using MATLAB functions*/
    if (utIsInterruptPending()) { // Handle MATLAB interrupt using CTRL+C
        utSetInterruptPending(false); // consume event and returning to MATLAB w/o interrupting main script
        return true; // stop IPOPT
    }
    return false;
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) {

    /* Assert number of input elements */
    ASSERTIN(nrhs!=1);

    /* Get input */
    const mxArray* input = prhs[0];
    
    /* Assert input type */
    ASSERTSTRUCT(input, (OCPInterface*) NULL);

    /* Assert input dimensions */
    ASSERTSIZE(input, 1, 1, (OCPInterface*) NULL);

    /* Get required fields */
    GETFIELD(input, name, false, (OCPInterface*) NULL);
    GETFIELD(input, N, false, (OCPInterface*) NULL);
    GETFIELD(input, ti, false, (OCPInterface*) NULL);
    GETFIELD(input, tf, false, (OCPInterface*) NULL);
    GETFIELD(input, bounds, false, (OCPInterface*) NULL);
    GETFIELD(input, guess, false, (OCPInterface*) NULL);

    /* Get optional fields */
    GETFIELD(input, options, true, (OCPInterface*) NULL);
    GETFIELD(input, auxdata, true, (OCPInterface*) NULL);
    GETFIELD(input, mesh, true, (OCPInterface*) NULL);

    /* Check name, N, ti, tf */
    ASSERTSTRING(name, 1, (OCPInterface*) NULL);
    //ASSERTPOSINTSCALAR(N, (OCPInterface*) NULL); // checked in constructor
    ASSERTSCALAR(ti, (OCPInterface*) NULL);
    ASSERTSCALAR(tf, (OCPInterface*) NULL);
    std::string namestr(mxArrayToString(name));

    /* Create OCPInterface object */
    OCPInterface* ocp;
    try {
        // this throws an excpetion if load of library fails 
        ocp = new OCPInterface(namestr, (int) GETVAL(N), GETVAL(ti), GETVAL(tf));
    } catch (const std::exception& e) {
        // terminate if errors
        mexErrMsgIdAndTxt("minosMex:loadLibFailed","%s", e.what());
    }

    /* Set the printing function to mexPrintf with workaround to print immediately instead of default one (printf) */
    ocp->set_printfun(&print);

    /* Set the interrupt handling function */
    ocp->set_interruptfun(&handle_interrupt);

    /* Get dimensions */
    int nx, nu, np, nc, nb, nq, nz, ng, nnzj, nnzh, na;
    ocp->get_dims(&nx, &nu, &np, &nc, &nb, &nq, &nz, &ng, &nnzj, &nnzh, &na);

    /* Check sub-struct */
    ASSERTSTRUCT(bounds, ocp);
    ASSERTSTRUCT(guess, ocp);
    ASSERTSTRUCT(options, ocp);
    ASSERTDOUBLE(auxdata, ocp);
    ASSERTSIZE(bounds, 1, 1, ocp);
    ASSERTSIZE(guess, 1, 1, ocp);
    ASSERTSIZE(options, 1, 1, ocp);
    ASSERTOPTNUMEL(auxdata, na, ocp);
    ASSERTDOUBLE(mesh, ocp); ASSERTNUMEL(mesh, (((int) GETVAL(N))-1), ocp); ASSERTSUM(mesh, 1.0, 1.0e-9, ocp);
    
    /* Get required fields */
    GETFIELD(bounds, lbx, false, ocp);
    GETFIELD(bounds, ubx, false, ocp);
    GETFIELD(bounds, lbu, false, ocp);
    GETFIELD(bounds, ubu, false, ocp);
    GETFIELD(bounds, lbp, false, ocp);
    GETFIELD(bounds, ubp, false, ocp);
    GETFIELD(bounds, lbc, false, ocp);
    GETFIELD(bounds, ubc, false, ocp);
    GETFIELD(bounds, lbb, false, ocp);
    GETFIELD(bounds, ubb, false, ocp);
    GETFIELD(bounds, lbq, false, ocp);
    GETFIELD(bounds, ubq, false, ocp);
    GETFIELD(guess, x, false, ocp);
    GETFIELD(guess, u, false, ocp);
    GETFIELD(guess, p, false, ocp);

    /* Get optional fields */
    GETFIELD(guess, lam_x, true, ocp);
    GETFIELD(guess, lam_u, true, ocp);
    GETFIELD(guess, lam_p, true, ocp);
    GETFIELD(guess, lam_f, true, ocp);
    GETFIELD(guess, lam_c, true, ocp);
    GETFIELD(guess, lam_b, true, ocp);
    GETFIELD(guess, lam_q, true, ocp);
    GETFIELD(options, max_iter, true, ocp);
    GETFIELD(options, flag_hessian, true, ocp);
    GETFIELD(options, mu_init, true, ocp);
    GETFIELD(options, display, true, ocp);
    GETFIELD(options, print_itersol, true, ocp);
    GETFIELD(options, outfile, true, ocp);
    GETFIELD(options, logfile, true, ocp);
    GETFIELD(options, nlpsolver, true, ocp);
    
    /* Check sub-fields */
    ASSERTDOUBLE(lbx, ocp); ASSERTNUMEL(lbx, nx, ocp);
    ASSERTDOUBLE(ubx, ocp); ASSERTNUMEL(ubx, nx, ocp);
    ASSERTDOUBLE(lbu, ocp); ASSERTNUMEL(lbu, nu, ocp);
    ASSERTDOUBLE(ubu, ocp); ASSERTNUMEL(ubu, nu, ocp);
    ASSERTDOUBLE(lbp, ocp); ASSERTNUMEL(lbp, np, ocp);
    ASSERTDOUBLE(ubp, ocp); ASSERTNUMEL(ubp, np, ocp);
    ASSERTDOUBLE(lbc, ocp); ASSERTNUMEL(lbc, nc, ocp);
    ASSERTDOUBLE(ubc, ocp); ASSERTNUMEL(ubc, nc, ocp);
    ASSERTDOUBLE(lbb, ocp); ASSERTNUMEL(lbb, nb, ocp);
    ASSERTDOUBLE(ubb, ocp); ASSERTNUMEL(ubb, nb, ocp);
    ASSERTDOUBLE(lbq, ocp); ASSERTNUMEL(lbq, nq, ocp);
    ASSERTDOUBLE(ubq, ocp); ASSERTNUMEL(ubq, nq, ocp);
    ASSERTDOUBLE(x, ocp); ASSERTSIZE(x, nx, ((int) GETVAL(N)), ocp);
    ASSERTDOUBLE(u, ocp); ASSERTSIZE(u, nu, ((int) GETVAL(N)), ocp);
    ASSERTDOUBLE(p, ocp); ASSERTNUMEL(p, np, ocp);
    ASSERTDOUBLE(lam_x, ocp); ASSERTSIZE(lam_x, nx, ((int) GETVAL(N)), ocp);
    ASSERTDOUBLE(lam_u, ocp); ASSERTSIZE(lam_u, nu, ((int) GETVAL(N)), ocp);
    ASSERTDOUBLE(lam_p, ocp); ASSERTNUMEL(lam_p, np, ocp);
    ASSERTDOUBLE(lam_f, ocp); ASSERTSIZE(lam_f, nx, (static_cast<size_t>((int) GETVAL(N)) - 1), ocp);
    ASSERTDOUBLE(lam_c, ocp); ASSERTSIZE(lam_c, nc, ((int) GETVAL(N)), ocp);
    ASSERTDOUBLE(lam_b, ocp); ASSERTNUMEL(lam_b, nb, ocp);
    ASSERTDOUBLE(lam_q, ocp); ASSERTNUMEL(lam_q, nq, ocp);
    ASSERTPOSINTSCALAR(max_iter, ocp);
    if (GETVAL(mu_init)!= -1) { 
        ASSERTPOSSCALAR(mu_init, ocp);
    }
    else {
        mu_init = NULL; // for default when -1  
    };   
    ASSERTLOGICAL(flag_hessian, ocp);
    ASSERTLOGICAL(display, ocp);
    ASSERTNONEGINTSCALAR(print_itersol, ocp);
    ASSERTSTRING(outfile, 1, ocp);
    ASSERTSTRING(logfile, 0, ocp);
    ASSERTSTRING(nlpsolver, 1, ocp);

    /* Set bounds */
    ocp->set_bounds(GETPTR(lbx), GETPTR(ubx),
                    GETPTR(lbu), GETPTR(ubu), 
                    GETPTR(lbp), GETPTR(ubp), 
                    GETPTR(lbc), GETPTR(ubc), 
                    GETPTR(lbb), GETPTR(ubb), 
                    GETPTR(lbq), GETPTR(ubq)
                    );

    /* Set guess */
    ocp->set_guess(GETPTR(x), GETPTR(u), GETPTR(p),
                   GETPTR(lam_x), GETPTR(lam_u), GETPTR(lam_p),
                   GETPTR(lam_f), GETPTR(lam_c), 
                   GETPTR(lam_b), GETPTR(lam_q)
                   );

    /* Set auxdata */
    ocp->set_auxdata(GETPTR(auxdata));

    /* Set mesh */
    ocp->set_mesh(GETPTR(mesh));

    /* Set options */
    if (max_iter)       ocp->set_option(OCPInterface::MAX_ITER, GETVAL(max_iter));
    if (mu_init)        ocp->set_option(OCPInterface::MU_INIT, GETVAL(mu_init));
    if (flag_hessian)   ocp->set_option(OCPInterface::FLAG_HESSIAN, GETVAL(flag_hessian));
    if (display)        ocp->set_option(OCPInterface::DISPLAY, GETVAL(display));
    if (print_itersol)  ocp->set_option(OCPInterface::PRINT_ITERSOL, GETVAL(print_itersol));
    if (logfile)        ocp->set_option(OCPInterface::LOGFILE, mxArrayToString(logfile));
    if (nlpsolver)      ocp->set_option(OCPInterface::NLPSOLVER, mxArrayToString(nlpsolver));

    /* Create result output matrices */
    CREATEDOUBLE(objval, 1, 1);
    CREATEDOUBLE(t, 1, ((int) GETVAL(N)));
    CREATEDOUBLE(xopt, nx, ((int) GETVAL(N)));
    CREATEDOUBLE(uopt, nu, ((int) GETVAL(N)));
    CREATEDOUBLE(popt, np, 1);
    CREATEDOUBLE(lam_xopt, nx, ((int) GETVAL(N)));
    CREATEDOUBLE(lam_uopt, nu, ((int) GETVAL(N)));
    CREATEDOUBLE(lam_popt, np, 1);
    CREATEDOUBLE(lam_fopt, nx, ( static_cast<size_t>((int) GETVAL(N)) - 1));
    CREATEDOUBLE(lam_copt, nc, ((int) GETVAL(N)));
    CREATEDOUBLE(lam_bopt, nb, 1);
    CREATEDOUBLE(lam_qopt, nq, 1);
    CREATEDOUBLE(fopt, nx, ( static_cast<size_t>((int) GETVAL(N)) - 1));
    CREATEDOUBLE(copt, nc, ((int) GETVAL(N)));
    CREATEDOUBLE(bopt, nb, 1);
    CREATEDOUBLE(qopt, nq, 1);
    CREATEDOUBLE(lopt, 1, ( static_cast<size_t>((int) GETVAL(N)) - 1));
    CREATEDOUBLE(mopt, 1, 1);
    CREATEDOUBLE(meshopt, 1, ( static_cast<size_t>((int) GETVAL(N)) - 1));
    CREATEDOUBLE(vj, nnzj, 1); CREATEINT32(irj, nnzj, 1); CREATEINT32(jcj, nnzj, 1);
    CREATEDOUBLE(vh, nnzh, 1); CREATEINT32(irh, nnzj, 1); CREATEINT32(jch, nnzj, 1);
    

    /* Call NLP solver */
    int status;
    try {
        // this throws an error if NLp lib not found
        status = ocp->solve();
    } catch (const std::exception& e) {
        // terminate if errors
        delete ocp;
        mexErrMsgIdAndTxt("minosMex:loadLibFailed","%s", e.what());
    }
    
    /* Get solution */
    ocp->get_sol(GETPTR(objval), GETPTR(t), 
                 GETPTR(xopt), GETPTR(uopt), GETPTR(popt),
                 GETPTR(lam_xopt), GETPTR(lam_uopt), GETPTR(lam_popt),
                 GETPTR(lam_fopt), GETPTR(lam_copt), GETPTR(lam_bopt), GETPTR(lam_qopt),
                 GETPTR(fopt), GETPTR(copt), GETPTR(bopt), GETPTR(qopt), 
                 GETPTR(lopt), GETPTR(mopt),
                 NULL, NULL, NULL,
                 (int*) GETPTR(irj), (int*) GETPTR(jcj), GETPTR(vj),
                 (int*) GETPTR(irh), (int*) GETPTR(jch), GETPTR(vh)
                 );
    
    /* Create sparse matrices */
    CREATESPARSE(jac, ng, nz, nnzj); 
    CREATESPARSE(hess, nz, nz, nnzh); 
    // convert COO to CSC format
    OCPInterface::coo2csc(ng, nz, nnzj, 
                         (int*) GETPTR(irj), (int*) GETPTR(jcj), GETPTR(vj),
                         GETIR(jac), GETJC(jac), GETPTR(jac));
    OCPInterface::coo2csc(nz, nz, nnzh, 
                         (int*) GETPTR(irh), (int*) GETPTR(jch), GETPTR(vh),
                         GETIR(hess), GETJC(hess), GETPTR(hess));

    /* Get mesh */
    ocp->get_mesh(GETPTR(meshopt));
    
    /* Get and print elapsed time */
    double ttot, talg, teval;
    ocp->get_cpu_time(ttot, talg, teval);

    /* Get and print num of iterations */
    int num_iter = ocp->get_num_iter();

    /* Get current barrier parameter */
    double mu_curr = ocp->get_mu_curr();
    
    /* Outfile */
    if (outfile) {
        std::ofstream file;
        char* buf = mxArrayToString(outfile);
        file.open(buf);
        if ( file.is_open() ) {
            file << ocp;
            file.close();
        } else {
            mexWarnMsgIdAndTxt("minosMex:failOpenFile", "Failed to create file %s.", buf);
        }
    }

    /* Get convergence history */
    CREATEDOUBLE(obj_history, static_cast<size_t>(num_iter)+1, 1);
    CREATEDOUBLE(infpr_history, static_cast<size_t>(num_iter)+1, 1);
    CREATEDOUBLE(infdu_history, static_cast<size_t>(num_iter)+1, 1);
    ocp->get_history(GETPTR(obj_history), GETPTR(infpr_history), GETPTR(infdu_history));

    /* Get employed options */
    if (!max_iter) {
        max_iter = mxCreateDoubleScalar(0);
        ocp->get_option(OCPInterface::MAX_ITER, GETPTR(max_iter));
    }
    if (!flag_hessian) {
        double val;
        ocp->get_option(OCPInterface::FLAG_HESSIAN, &val);
        flag_hessian = mxCreateLogicalScalar(val>0);
    }
    if (!display) {
        double val;
        ocp->get_option(OCPInterface::DISPLAY, &val);
        display = mxCreateLogicalScalar(val>0);
    }
    if (!print_itersol) {
        print_itersol = mxCreateDoubleScalar(0);
        ocp->get_option(OCPInterface::PRINT_ITERSOL, GETPTR(print_itersol));
    }
    if (!logfile) {
        std::string str;
        ocp->get_option(OCPInterface::LOGFILE, str);
        logfile = mxCreateString(str.c_str());
    }
    if (!nlpsolver) {
        std::string str;
        ocp->get_option(OCPInterface::NLPSOLVER, str);
        nlpsolver = mxCreateString(str.c_str());
    }

    /* Create stats struct */
    const char* stats_fields[] = {"num_iter", "obj_history", "infpr_history", "infdu_history", "mu_curr", 
                                  "nz", "ng", "nnzj", "nnzh", "jac", "hess",
                                  "ttot", "talg", "teval"};
    mxArray* stats_vars[] = {mxCreateDoubleScalar(num_iter), obj_history, infpr_history, infdu_history, mxCreateDoubleScalar(mu_curr), 
                             mxCreateDoubleScalar(nz), mxCreateDoubleScalar(ng), mxCreateDoubleScalar(nnzj), mxCreateDoubleScalar(nnzh),
                             jac, hess, 
                             mxCreateDoubleScalar(ttot), mxCreateDoubleScalar(talg), mxCreateDoubleScalar(teval)};
    CREATESTRUCT(stats, 1, 1, sizeof(stats_fields)/sizeof(void*), stats_fields);
    SETFIELDS(stats, stats_fields, stats_vars);
    
    /* Create next_problem struct */
    mxArray* next_bounds = mxDuplicateArray(bounds); // use bounds of input argument
    const char* opt_fields[] = {"nlpsolver", "max_iter", "flag_hessian", "mu_init", "display", "print_itersol", "logfile"};
    mxArray* opt_vars[] = {mxDuplicateArray(nlpsolver), mxDuplicateArray(max_iter), mxDuplicateArray(flag_hessian), mxCreateDoubleScalar(mu_curr), mxDuplicateArray(display), mxDuplicateArray(print_itersol), mxDuplicateArray(logfile)};
    CREATESTRUCT(next_options, 1, 1, sizeof(opt_fields)/sizeof(void*), opt_fields);
    SETFIELDS(next_options, opt_fields, opt_vars);
    if (outfile) {
        mxAddField(next_options, "outfile");
        mxSetField(next_options, 0, "outfile", mxDuplicateArray(outfile));
    }
    const char *guess_fields[] = {"x", "u", "p", 
                                  "lam_x", "lam_u", "lam_p", 
                                  "lam_f", "lam_c", "lam_b"
                                  };
    mxArray* guess_vars[] = {mxDuplicateArray(xopt), mxDuplicateArray(uopt), mxDuplicateArray(popt), 
                             mxDuplicateArray(lam_xopt), mxDuplicateArray(lam_uopt), mxDuplicateArray(lam_popt), 
                             mxDuplicateArray(lam_fopt), mxDuplicateArray(lam_copt), mxDuplicateArray(lam_bopt), mxDuplicateArray(lam_qopt)
                            };
    CREATESTRUCT(next_guess, 1, 1, sizeof(guess_fields)/sizeof(void*), guess_fields);
    SETFIELDS(next_guess, guess_fields, guess_vars);
    mwSize auxdata_dims[2] = {static_cast<size_t>(na), 1};
    mxArray* next_auxdata = (auxdata) ? mxDuplicateArray(auxdata) : mxCreateNumericArray(2, auxdata_dims, mxDOUBLE_CLASS, mxREAL);
    const char* prob_fields[] = {"name", "N", "ti", "tf", 
                                 "guess", "bounds", "auxdata", "mesh", "options"};
    mxArray* prob_vars[] = {mxDuplicateArray(name), mxDuplicateArray(N), mxDuplicateArray(ti), mxDuplicateArray(tf),
                            next_guess, next_bounds, next_auxdata, mxDuplicateArray(meshopt), next_options};
    CREATESTRUCT(next_problem, 1, 1, sizeof(prob_fields)/sizeof(void*), prob_fields);
    SETFIELDS(next_problem, prob_fields, prob_vars);

    /* Create output struct */
    const char *out_fields[] = {"exitflag", "objval", "t",
                                "x", "u", "p", 
                                "lam_x", "lam_u", "lam_p", 
                                "lam_f", "lam_c", "lam_b", "lam_q",
                                "f", "c", "b", "q",
                                "l", "m",
                                "stats", "next_problem",
                                };
    mxArray* out_vars[] = {mxCreateDoubleScalar(status), objval, t,
                          xopt, uopt, popt, 
                          lam_xopt, lam_uopt, lam_popt, 
                          lam_fopt, lam_copt, lam_bopt, lam_qopt, 
                          fopt, copt, bopt, qopt, 
                          lopt, mopt,
                          stats, next_problem
                          };
    CREATESTRUCT(out, 1, 1, sizeof(out_fields)/sizeof(void*), out_fields);
    SETFIELDS(out, out_fields, out_vars);
    
    /* Set output */
    int c = 0; // output counter
    OUTPUT(out, c);

    /* Assert number of output elements */
    ASSERTOUT(nlhs>c);

    /* Free mem */
    delete ocp;
}
