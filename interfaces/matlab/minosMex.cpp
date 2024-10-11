/*  File minosMex.cpp 
    Mex interface for MinOS.
    Copyright (C) 2024 Stefano Lovato
*/

/* Macros */
#define ASSERTIN(x) (x) ? mexErrMsgIdAndTxt("MinOS:wrongNumInputs", "One input argument required.\n") : (void(0))
#define ASSERTOUT(x) (x) ? mexErrMsgIdAndTxt("MinOS:wrongNumout_vars", "Too many output arguments.\n") : (void(0))
#define ASSERTSTRUCT(x) ( (x) && (!mxIsStruct(x))) ? mexErrMsgIdAndTxt("ocpSOlver:wrongInputType", "Input must be a struct\n") : (void(0))
#define ASSERTPOSINTSCALAR(x) ( (x) && (!mxIsDouble(x) || \
                            mxIsComplex(x) || \
                            mxGetNumberOfElements(x)!=1 || \
                            mxGetScalar(x)<=0 || \
                            ((((int) mxGetScalar(x)) != mxGetScalar(x)) != 0) ) ) ? \
                            mexErrMsgIdAndTxt("MinOS:notScalar","%s must be a positive integer scalar.", #x) : (void(0))
#define ASSERTNONEGINTSCALAR(x) ( (x) && (!mxIsDouble(x) || \
                            mxIsComplex(x) || \
                            mxGetNumberOfElements(x)!=1 || \
                            mxGetScalar(x)<0 || \
                            ((((int) mxGetScalar(x)) != mxGetScalar(x)) != 0) ) ) ? \
                            mexErrMsgIdAndTxt("MinOS:notScalar","%s must be a positive integer scalar.", #x) : (void(0))
#define ASSERTPOSSCALAR(x) ( (x) && (!mxIsDouble(x) || \
                            mxIsComplex(x) || \
                            mxGetNumberOfElements(x)!=1 || \
                            mxGetScalar(x)<=0 ) ) ? \
                            mexErrMsgIdAndTxt("MinOS:notScalar","%s must be a positive scalar.", #x) : (void(0))
#define ASSERTSCALAR(x) ( (x) && ( !mxIsDouble(x) || \
                            mxIsComplex(x) || \
                            mxGetNumberOfElements(x)!=1 ) ) ? \
                            mexErrMsgIdAndTxt("MinOS:notScalar","%s must be a scalar.", #x) : (void(0))
#define ASSERTGREATER(x1,x2) ( mxGetScalar(x1) <= mxGetScalar(x2) ) ? \
                            mexErrMsgIdAndTxt("MinOS:notGreater","%s must be greater than %s.", #x1, #x2) : (void(0))
#define ASSERTDOUBLE(x)   ( (x) && (!mxIsDouble(x) || mxIsComplex(x)) ) ? \
                            mexErrMsgIdAndTxt("MinOS:notDouble","%s must be of type double.", #x) : (void(0))
#define ASSERTLOGICAL(x) ( (x) && (!mxIsLogicalScalar(x)) ) ? \
                            mexErrMsgIdAndTxt("MinOS:notLogical","%s must be a logical scalar.", #x) : (void(0))
#define ASSERTSIZE(x,n,m) ( (x) && \
                            ((mxGetM(x) != n) || (mxGetN(x) != m)) ) ? \
                            mexErrMsgIdAndTxt("MinOS:wrongSize","%s must be a have dimension %d-by-%d (found %d-by-%d).",#x,n,m,mxGetM(x),mxGetN(x)) : (void(0))
#define ASSERTNUMEL(x,n)  ( (x) && \
                            (mxGetNumberOfElements(x) != n) ) ? \
                            mexErrMsgIdAndTxt("MinOS:wrongSize","%s must be have %d elements (found %d).", #x,n,mxGetNumberOfElements(x)) : (void(0))
#define ASSERTSTRING(x,n)   ( (x) && (!mxIsChar(x) || mxGetM(x) != n) ) ? \
                            mexErrMsgIdAndTxt("MinOS:notString","%s must be a string.", #x) : (void(0))
#define ASSERTOPTIONALNUMEL(x,n)  ( !(x) && (n > 0) ) ? \
                            mexErrMsgIdAndTxt("MinOS:missingField", "Field %s with %d elements is required.", #x, n) : \
                            ( (x) && \
                            (mxGetNumberOfElements(x) != n) ) ? \
                            mexErrMsgIdAndTxt("MinOS:wrongSize","%s must be have %d elements (found %d).", #x,n,mxGetNumberOfElements(x)) : (void(0))
#define ASSERT(x)  ( ! (x) ) ? \
                   m        exErrMsgIdAndTxt("MinOS:falseReturn", "%s return false.", #x) : (void(0))

#define GETFIELD(p,name,opt) mxArray* name = (p) ? mxGetField(p,0,#name) : NULL; \
                            ( (name == NULL) & !opt ) ? mexErrMsgIdAndTxt("MinOS:missingField", "Missing field %s in %s.", #name, #p) : (void(0))
#define GETVAL(x) ((x) ? mxGetScalar(x) : 0)
#define GETPTR(x) ((x) ? mxGetPr(x) : NULL)
#define GETIR(x) ((x) ? mxGetIr(x) : NULL)
#define GETJC(x) ((x) ? mxGetJc(x) : NULL)

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
#define ASSERTSUM(x,s,e) if (x) { \
                            double sum = 0; \
                            double* xp = GETPTR(x); \
                            for (int i = 0; i < mxGetNumberOfElements(x); ++i) sum += xp[i]; \
                            if (((sum-s) > e) || ((sum-s) < -e)) mexErrMsgIdAndTxt("MinOS:invalidField", "Field %s must sum to %f (sum to %f).\n", #x, s, sum); \
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
int print(const char *fmt,
         ...
) {
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
    ASSERTSTRUCT(input);

    /* Assert input dimensions */
    ASSERTSIZE(input, 1, 1);

    /* Get required fields */
    GETFIELD(input, name, false);
    GETFIELD(input, N, false);
    GETFIELD(input, ti, false);
    GETFIELD(input, tf, false);
    GETFIELD(input, bounds, false);
    GETFIELD(input, guess, false);

    /* Get optional fields */
    GETFIELD(input, options, true);
    GETFIELD(input, auxdata, true);
    GETFIELD(input, mesh, true);

    /* Check name, N, ti, tf */
    ASSERTSTRING(name,1);
    ASSERTPOSINTSCALAR(N);
    ASSERTSCALAR(ti);
    ASSERTSCALAR(tf);
    ASSERTGREATER(tf,ti);
    std::string namestr(mxArrayToString(name));

    /* Create OCPInterface object */
    OCPInterface* ocp;
    try {
        // this throws an excpetion if load of library fails 
        ocp = new OCPInterface(namestr, (int) GETVAL(N), GETVAL(ti), GETVAL(tf));
    } catch (const std::exception& e) {
        // terminate if errors
        mexErrMsgIdAndTxt("MinOS:loadLibFailed","%s", e.what());
    }

    /* Set the printing function to mexPrintf with workaround to print immediately instead of default one (printf) */
    ocp->set_printfun(&print);

    /* Set the interrupt handling function */
    ocp->set_interruptfun(&handle_interrupt);

    /* Get dimensions */
    int nx, nu, np, nc, nb, nq, nz, ng, nnzj, nnzh, na;
    ocp->get_dims(&nx, &nu, &np, &nc, &nb, &nq, &nz, &ng, &nnzj, &nnzh, &na);

    /* Check sub-struct */
    ASSERTSTRUCT(bounds);
    ASSERTSTRUCT(guess);
    ASSERTSTRUCT(options);
    ASSERTDOUBLE(auxdata);
    ASSERTSIZE(bounds, 1, 1);
    ASSERTSIZE(guess, 1, 1);
    ASSERTSIZE(options, 1, 1);
    ASSERTOPTIONALNUMEL(auxdata, na);
    ASSERTDOUBLE(mesh); ASSERTNUMEL(mesh, (((int) GETVAL(N))-1)); ASSERTSUM(mesh, 1.0, 1.0e-9);
    
    /* Get required fields */
    GETFIELD(bounds, lbx, false);
    GETFIELD(bounds, ubx, false);
    GETFIELD(bounds, lbu, false);
    GETFIELD(bounds, ubu, false);
    GETFIELD(bounds, lbp, false);
    GETFIELD(bounds, ubp, false);
    GETFIELD(bounds, lbc, false);
    GETFIELD(bounds, ubc, false);
    GETFIELD(bounds, lbb, false);
    GETFIELD(bounds, ubb, false);
    GETFIELD(bounds, lbq, false);
    GETFIELD(bounds, ubq, false);
    GETFIELD(guess, x, false);
    GETFIELD(guess, u, false);
    GETFIELD(guess, p, false);

    /* Get optional fields */
    GETFIELD(guess, lam_x, true);
    GETFIELD(guess, lam_u, true);
    GETFIELD(guess, lam_p, true);
    GETFIELD(guess, lam_f, true);
    GETFIELD(guess, lam_c, true);
    GETFIELD(guess, lam_b, true);
    GETFIELD(guess, lam_q, true);
    GETFIELD(options, max_iter, true);
    GETFIELD(options, flag_hessian, true);
    GETFIELD(options, mu_init, true);
    GETFIELD(options, display, true);
    GETFIELD(options, print_itersol, true);
    GETFIELD(options, outfile, true);
    GETFIELD(options, logfile, true);
    GETFIELD(options, nlpsolver, true);
    
    /* Check sub-fields */
    ASSERTDOUBLE(lbx); ASSERTNUMEL(lbx, nx);
    ASSERTDOUBLE(ubx); ASSERTNUMEL(ubx, nx);
    ASSERTDOUBLE(lbu); ASSERTNUMEL(lbu, nu);
    ASSERTDOUBLE(ubu); ASSERTNUMEL(ubu, nu);
    ASSERTDOUBLE(lbp); ASSERTNUMEL(lbp, np);
    ASSERTDOUBLE(ubp); ASSERTNUMEL(ubp, np);
    ASSERTDOUBLE(lbc); ASSERTNUMEL(lbc, nc);
    ASSERTDOUBLE(ubc); ASSERTNUMEL(ubc, nc);
    ASSERTDOUBLE(lbb); ASSERTNUMEL(lbb, nb);
    ASSERTDOUBLE(ubb); ASSERTNUMEL(ubb, nb);
    ASSERTDOUBLE(lbq); ASSERTNUMEL(lbq, nq);
    ASSERTDOUBLE(ubq); ASSERTNUMEL(ubq, nq);
    ASSERTDOUBLE(x); ASSERTSIZE(x, nx, ((int) GETVAL(N)));
    ASSERTDOUBLE(u); ASSERTSIZE(u, nu, ((int) GETVAL(N)));
    ASSERTDOUBLE(p); ASSERTNUMEL(p, np);
    ASSERTDOUBLE(lam_x); ASSERTSIZE(lam_x, nx, ((int) GETVAL(N)));
    ASSERTDOUBLE(lam_u); ASSERTSIZE(lam_u, nu, ((int) GETVAL(N)));
    ASSERTDOUBLE(lam_p); ASSERTNUMEL(lam_p, np);
    ASSERTDOUBLE(lam_f); ASSERTSIZE(lam_f, nx, (static_cast<size_t>((int) GETVAL(N)) - 1));
    ASSERTDOUBLE(lam_c); ASSERTSIZE(lam_c, nc, ((int) GETVAL(N)));
    ASSERTDOUBLE(lam_b); ASSERTNUMEL(lam_b, nb);
    ASSERTDOUBLE(lam_q); ASSERTNUMEL(lam_q, nq);
    ASSERTPOSINTSCALAR(max_iter);
    if (GETVAL(mu_init)!= -1) ASSERTPOSSCALAR(mu_init);
    else mu_init = NULL; // for default when -1    
    ASSERTLOGICAL(flag_hessian);
    ASSERTLOGICAL(display);
    ASSERTNONEGINTSCALAR(print_itersol);
    ASSERTSTRING(outfile,1);
    ASSERTSTRING(logfile,0);
    ASSERTSTRING(nlpsolver,1);

    /* Set defaults */
    if (!max_iter) max_iter = mxCreateDoubleScalar(3000); // set default
    if (!mu_init) mu_init = mxCreateDoubleScalar(-1); // set default
    if (!flag_hessian) flag_hessian = mxCreateLogicalScalar(false); // set default
    if (!display) display = mxCreateLogicalScalar(true); // set default
    if (!print_itersol) print_itersol = mxCreateDoubleScalar(0); // set default
    if (!logfile) logfile = mxCreateString(""); // set default
    if (!nlpsolver) nlpsolver = mxCreateString("ipopt"); // set default

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
    ocp->set_option(OCPInterface::MAX_ITER, GETVAL(max_iter));
    ocp->set_option(OCPInterface::MU_INIT, GETVAL(mu_init));
    ocp->set_option(OCPInterface::FLAG_HESSIAN, GETVAL(flag_hessian));
    ocp->set_option(OCPInterface::DISPLAY, GETVAL(display));
    ocp->set_option(OCPInterface::PRINT_ITERSOL, GETVAL(print_itersol));
    ocp->set_option(OCPInterface::LOGFILE, mxArrayToString(logfile));
    ocp->set_option(OCPInterface::NLPSOLVER, mxArrayToString(nlpsolver));

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

    /* Call NLP solver */
    int status;
    try {
        // this throws an error if NLp lib not found
        status = ocp->solve();
    } catch (const std::exception& e) {
        // terminate if errors
        delete ocp;
        mexErrMsgIdAndTxt("MinOS:loadLibFailed","%s", e.what());
    }

    /* Get solution */
    ocp->get_sol(GETPTR(objval), GETPTR(t), 
                 GETPTR(xopt), GETPTR(uopt), GETPTR(popt),
                 GETPTR(lam_xopt), GETPTR(lam_uopt), GETPTR(lam_popt),
                 GETPTR(lam_fopt), GETPTR(lam_copt), GETPTR(lam_bopt), GETPTR(lam_qopt),
                 GETPTR(fopt), GETPTR(copt), GETPTR(bopt), GETPTR(qopt), 
                 GETPTR(lopt), GETPTR(mopt)
                 );

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
            mexWarnMsgIdAndTxt("MinOS:failOpenFile", "Failed to create file %s.", buf);
        }
    }

    /* Get convergence history */
    CREATEDOUBLE(obj_history, static_cast<size_t>(num_iter)+1, 1);
    CREATEDOUBLE(infpr_history, static_cast<size_t>(num_iter)+1, 1);
    CREATEDOUBLE(infdu_history, static_cast<size_t>(num_iter)+1, 1);
    ocp->get_history(GETPTR(obj_history), GETPTR(infpr_history), GETPTR(infdu_history));

    /* Create stats struct */
    const char* stats_fields[] = {"num_iter", "obj_history", "infpr_history", "infdu_history", "mu_curr", 
                                  "nz", "ng", "nnzj", "nnzh", 
                                  "ttot", "talg", "teval"};
    mxArray* stats_vars[] = {mxCreateDoubleScalar(num_iter), obj_history, infpr_history, infdu_history, mxCreateDoubleScalar(mu_curr), 
                             mxCreateDoubleScalar(nz), mxCreateDoubleScalar(ng), mxCreateDoubleScalar(nnzj), mxCreateDoubleScalar(nnzh), 
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
