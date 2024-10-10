/*  File minos_base.cpp 
    OCPInterface class implementation with NLP methods
    Copyright (C) 2024 Stefano Lovato
*/

#include "minos.h" 
#include "macros.h" // for nlp macros

/** Solve the NLP using the specified NLP solver */
int OCPInterface::solve(
){
    /* Print MinOS version */
    if (display)
        (*print_funptr)("%s\n", get_header().c_str());

    /* Reset is_interrupt_requested */
    OCPInterface::is_interrupt_requested = false;

    /* Load NLP solver library */
    void* nlp_lib = NULL;
    SolveFunc nlp_solve = NULL;
    nlp_lib = load_nlplib("minos-" + nlpsolver);
    nlp_solve = reinterpret_cast<SolveFunc>(import_nlpsolve(nlp_lib));
    if (!nlp_solve) { return -1; }
    
    /* Lof filename */
    if (logfile.empty()) this->logfile = name + ".log"; // default logfile name is <name>.log

    /* Call to solve */
    int status = (*nlp_solve)(this);

    /* Free NLP library */
    free_library(nlp_lib);

    /* Return status */
    return status;
}

/** Init constant arguments */
void OCPInterface::init_args() {
    argrc[0] = &t;
    argrcg[0] = &t;
    argd[0] = &t;
    argp[0] = &t;
    argq[0] = &t;
    argdj[0] = &t;
    argpj[0] = &t;  
    argqj[0] = &t;
    argrc[5] = &h;
    argrcg[5] = &h;
    argd[5] = &h;
    argdj[5] = &h;
    argq[5] = &h;
    argqj[5] = &h;
    if (ocp_hessb && ocp_hessi) {
        arghi[0] = &t;
        arghi[5] = &h;
        arghb[0] = &tf;
    }
}

/** Init cost gradient */
void OCPInterface::init_cost_gradient(
) {
    // get nnz of running cost
    get_sizes(SPOUT(ocp_runcost_grad), 0, NULL, NULL, &nnzrcg);
    // get nnz of boundary cost
    get_sizes(SPOUT(ocp_bcscost_grad), 0, NULL, NULL, &nnzbcg);
    // get pattern of running cost
    irrcg = new int[1+nnzrcg] { 0 };
    resrcg[0] = new double[1+nnzrcg] { 0 }; // init resrcg
    get_pattern(SPOUT(ocp_runcost_grad), 0, irrcg, NULL); // use only row indexes (gradient is column matrix)
    // create krcg
    krcg = new int[1+nnzrcg*(N-1)] { -1 };
    int ir; // store temporary index
    for (int k = 0; k < N-1; ++k) {
        for (int i = 0; i < nnzrcg; ++i) {
            ir = irrcg[i];
            if (ir < (nx+nu+nx)) ir += k*(nx+nu); //x1,u,x2
            else ir += -(nx+nu+nx) + N*(nx+nu); // p
            krcg[k*nnzrcg + i] = ir;
        }
    }
    // get pattern of boundary cost
    irbcg = new int[1+nnzbcg] { 0 };
    resbcg[0] = new double[1+nnzbcg] { 0 }; // init resbcg
    get_pattern(SPOUT(ocp_bcscost_grad), 0, irbcg, NULL); // use only row indexes (gradient is column matrix)
    // create kbcg
    kbcg = new int[1+nnzbcg] { -1 };
    for (int i = 0; i < nnzbcg; ++i) {
        ir = irbcg[i];
        if (ir >= (nx+nu)) ir += -(nx+nu) + (N-1)*(nx+nu); // xf, uf, p
        kbcg[i] = ir;
    }
}

/** Init constraint jacobian */
void OCPInterface::init_constr_jac() {
    // get nnz of dyn, path, bcs
    get_sizes(SPOUT(ocp_dyn_jac), 0, NULL, NULL, &nnzdj); // dyn w.r.t. x1, u, x2, p
    get_sizes(SPOUT(ocp_path_jac), 0, NULL, NULL, &nnzpj);  // path w.r.t. x, u, p
    get_sizes(SPOUT(ocp_bcs_jac), 0, NULL, NULL, &nnzbj); // bcs w.r.t. x0, u0, xn, un, p
    get_sizes(SPOUT(ocp_int_jac), 0, NULL, NULL, &nnzqj); // int w.r.t. x1, u, x2, p

    // calc total nnz of jacobian, not including integral constraints
    nnzj = (N-1)*(nnzdj+nnzpj) + nnzpj + nnzbj; // total nnz of the jacobian not including integral constraints

    // get ir,jc of dyn, path, bcs
    irdj = new int[1+nnzdj] { 0 };
    jcdj = new int[1+nnzdj] { 0 };
    irpj = new int[1+nnzpj] { 0 };
    jcpj = new int[1+nnzpj] { 0 };
    irbj = new int[1+nnzbj] { 0 };
    jcbj = new int[1+nnzbj] { 0 };
    irqj = new int[1+nnzqj] { 0 };
    jcqj = new int[1+nnzqj] { 0 };
    get_pattern(SPOUT(ocp_dyn_jac), 0, irdj, jcdj); // dyn w.r.t. x1, u, x2, p
    get_pattern(SPOUT(ocp_path_jac), 0, irpj, jcpj); // path w.r.t. x, u, p
    get_pattern(SPOUT(ocp_bcs_jac), 0, irbj, jcbj); // bcs w.r.t. x0, u0, xn, un, p
    get_pattern(SPOUT(ocp_int_jac), 0, irqj, jcqj); // int w.r.t. x1, u, x2, p

    // procedire to compute the nnz of the entire int jac and corresponding indexes
    resq[0] = new double[1+nq] { 0 }; // init also resq for ocp_int
    resqj[0] = new double[1+nnzqj] { 0 }; // init resqj
    kjq = new int[1+nnzqj*(N-1)] { -1 };
    kj = new int[1+nnzqj*(N-1)] { -1 }; // save 'rolled' indexes of nz, with size from upper estimation of nnz
    int c_nnz = 0; // counter for int jacobian nz
    int ir, jc, kk; // temporary row, column, and rolled indexes
    for (int k = 0; k < N-1; ++k) { // iterate over mesh intervals (k from 0 to N-2)
        for (int i = 0; i < nnzqj; ++i) {
            ir = irqj[i]; jc = jcqj[i]; // get current ir,jc
            // convert local to global indexes
            if (jc < (nx+nu+nx)) jc += k*(nx+nu); //x1,u,x2
            else jc += -(nx+nu+nx) + N*(nx+nu); // p
            ir += (N-1)*(nx+nc) + nc + nb; // q
            kk = ir + nz*jc;
            int idx = search_value(kj, c_nnz, kk); // search kk in kh, return index position
            if (idx < 0) { // value not found: add new nz
                kj[c_nnz] = kk;
                kjq[k*nnzqj + i] = nnzj + c_nnz;
                ++c_nnz;
            } else { // value found: not add new nz
                kjq[k*nnzqj + i] = nnzj + idx;
            }
        }
    }
    // update nnzj
    nnzj += c_nnz;
}

/** Init lagragian hessian */
void OCPInterface::init_lag_hessian() {
    if (ocp_hessb && ocp_hessi) {
        // get nnz of hessb and hessi
        get_sizes(SPOUT(ocp_hessb), 0, NULL, NULL, &nnzhb); 
        get_sizes(SPOUT(ocp_hessi), 0, NULL, NULL, &nnzhi); 
        // init patterns
        irhb = new int[1+nnzhb] { 0 };
        jchb = new int[1+nnzhb] { 0 };
        irhi = new int[1+nnzhi] { 0 };
        jchi = new int[1+nnzhi] { 0 };
        // init res
        reshb[0] = new double[1+nnzhb] { 0 };
        reshi[0] = new double[1+nnzhi] { 0 };
        // get patterns
        get_pattern(SPOUT(ocp_hessb), 0, irhb, jchb);
        get_pattern(SPOUT(ocp_hessi), 0, irhi, jchi);
        // init khb, khi
        khb = new int[1+nnzhb] { -1 };
        khi = new int[1+nnzhi*(N-1)] { -1 };
        // prcedure to calculate nnzh
        nnzh = nnzhb + (N-1)*nnzhi; // upper estimation of nnz for init
        kh = new int[1+nnzh]; // save 'rolled' indexes of nz
        int c_nnz = 0; // counter for hessian nz
        int ir, jc, kk; // temporary row, column, and rolled indexes
        // nnz related to hessb
        for (int i = 0; i < nnzhb; ++i) {
            ir = irhb[i]; jc = jchb[i]; // get current ir,jc
            if (jc>ir) { continue; } // skip b/c not lower triangular
            // convert local ir,jc to global
            if (ir >= (nx+nu)) ir += -(nx+nu) + (N-1)*(nx+nu); // xn,un,p
            if (jc >= (nx+nu)) jc += -(nx+nu) + (N-1)*(nx+nu); // xn,un,p
            kk = ir + nz*jc; // current kk
            int idx = search_value(kh, c_nnz, kk); // search kk in kh, return index position
            if (idx < 0) { // value not found: add new nz
                kh[c_nnz] = kk;
                khb[i] = c_nnz;
                ++c_nnz;
            } else { // value found: not add new nz
                khb[i] = idx;
            }
        }
        // nnz related to hessi
        for (int k = 0; k < N-1; ++k) { // iterate over mesh intervals (k from 0 to N-2)
            for (int i = 0; i < nnzhi; ++i) {
                ir = irhi[i]; jc = jchi[i]; // get current ir,jc
                if (jc>ir) { continue; } // skip b/c not lower triangular
                // convert local ir,jc to global
                if (ir < (nx+nu+nx)) ir += k*(nx+nu); //x1,u,x2
                else ir += -(nx+nu+nx) + N*(nx+nu); // p
                if (jc < (nx+nu+nx)) jc += k*(nx+nu); //x1,u,x2
                else jc += -(nx+nu+nx) + N*(nx+nu); // p
                kk = ir + nz*jc; // current kk
                int idx = search_value(kh, c_nnz, kk); // search kk in kh, return index position
                if (idx < 0) { // value not found: add new nz
                    kh[c_nnz] = kk;
                    khi[k*nnzhi + i] = c_nnz;
                    ++c_nnz;
                } else { // value found: not add new nz
                    khi[k*nnzhi + i] = idx;
                }
            }
        }
        // update nnzh
        nnzh = c_nnz;
    } else nnzh = 0;
}

/** Get variable and constraint boundaries */
void OCPInterface::get_nlp_bounds(
    double* lbz,
    double* ubz,
    double* lbg,
    double* ubg
) {
    // z bounds
    // init pointers to data
    double* lbxdest = lbz;
    double* ubxdest = ubz;
    double* lbudest = lbz + nx;
    double* ubudest = ubz + nx;
    // iterate over mesh points (k from 0 to N-1)
    for (int k = 0; k < N; ++k) {
        // copy lbx, ubx, lbu, ubu
        if (lbz && lbx) memcpy(lbxdest, lbx, nx*sizeof(double));
        if (ubz && ubx) memcpy(ubxdest, ubx, nx*sizeof(double));
        if (lbz && lbu) memcpy(ubudest, ubu, nu*sizeof(double));
        if (ubz && ubu) memcpy(lbudest, lbu, nu*sizeof(double));
        // move pointers forwards
        lbxdest += nx+nu;
        lbudest += nx+nu;
        ubxdest += nx+nu;
        ubudest += nx+nu;
    }
    // copy lbp, ubp
    double* lbpdest = lbz + N*(nx+nu);
    double* ubpdest = ubz + N*(nx+nu);
    if (lbz && lbp) memcpy(lbpdest, lbp, np*sizeof(double));
    if (ubz && ubp) memcpy(ubpdest, ubp, np*sizeof(double));

    // g bounds
    // init pointers to data
    double* lbfdest = lbg;
    double* lbcdest = lbg + nx;
    double* ubfdest = ubg;
    double* ubcdest = ubg + nx;
    // iterate over mesh intervals (k from 0 to N-2)
    for (int k = 0; k < N-1; ++k) {
        // set lbf, ubf to zero
        if (lbg) memset((void*) lbfdest, 0, nx*sizeof(double));
        if (ubg) memset((void*) ubfdest, 0, nx*sizeof(double));
        // copy lbc, ubc
        if (lbg && lbc) memcpy(lbcdest, lbc, nc*sizeof(double));
        if (ubg && ubc) memcpy(ubcdest, ubc, nc*sizeof(double));
        // move pointers forwards
        lbfdest += nx+nc;
        lbcdest += nx+nc;
        ubfdest += nx+nc;
        ubcdest += nx+nc;
    }
    // copy lbc, ubc for last point
    lbcdest = lbg + (N-1)*(nx+nc);
    ubcdest = ubg + (N-1)*(nx+nc);
    if (lbg && lbc) memcpy(lbcdest, lbc, nc*sizeof(double));
    if (ubg && ubc) memcpy(ubcdest, ubc, nc*sizeof(double));
    // copy lbb, ubb
    double *lbbdest = lbg + (N-1)*(nx+nc) + nc;
    double *ubbdest = ubg + (N-1)*(nx+nc) + nc;
    if (lbg && lbb) memcpy(lbbdest, lbb, nb*sizeof(double));
    if (ubg && ubb) memcpy(ubbdest, ubb, nb*sizeof(double));
    // copy lbq, ubq
    double *lbqdest = lbg + (N-1)*(nx+nc) + nc + nb;
    double *ubqdest = ubg + (N-1)*(nx+nc) + nc + nb;
    if (lbg && lbq) memcpy(lbqdest, lbq, nq*sizeof(double));
    if (ubg && ubq) memcpy(ubqdest, ubq, nq*sizeof(double));
}

/** Get initial point */
void OCPInterface::get_nlp_init(
    bool init_z,
    double* z,
    bool init_lamz,
    double* lamz,
    bool init_lamg,
    double* lamg
) {
    // init z
    if (init_z) {
        // init pointers to data
        double* xdest = z;
        double* udest = z + nx;
        double* xsrc = x0;
        double* usrc = u0;
        // iterate over mesh points (k from 0 to N-1)
        for (int k = 0; k < N; ++k) {
            // copy x0[:,k] and u0[:,k]
            if (x0) memcpy(xdest, xsrc, nx*sizeof(double));
            if (u0) memcpy(udest, usrc, nu*sizeof(double)); 
            // move pointers forwards
            xdest += nx+nu;
            udest += nx+nu;
            xsrc += nx;
            usrc += nu;
        }
        // copy p0
        double* pdest = z + N*(nx+nu);
        double* psrc = p0;
        if (p0) memcpy(pdest, psrc, np*sizeof(double)); 
    }

    // init lamzl & lamzu
    if (init_lamz) {
        // init pointers to data
        double* lamxdest = lamz;
        double* lamudest = lamz + nx;
        double* lamxsrc = lam_x0;
        double* lamusrc = lam_u0;
        // iterate over mesh points (k from 0 to N-1)
        for (int k = 0; k < N; ++k) {
            // copy lamx0[:,k] and lamu0[:,k]
            if (lam_x0 && flag_lamx) memcpy(lamxdest, lamxsrc, nx*sizeof(double));
            if (lam_u0 && flag_lamu) memcpy(lamudest, lamusrc, nu*sizeof(double)); 
            // move pointers forwards
            lamxdest += nx+nu;
            lamudest += nx+nu;
            lamxsrc += nx;
            lamusrc += nu;
        }
        // calc lampl,lampu
        double* lampdest = lamz + N*(nx+nu);
        double* lampsrc = lam_p0;
        // copy lamp0
        if (lam_p0 && flag_lamp) memcpy(lampdest, lampsrc, np*sizeof(double));
    }

    // init lamg
    if (init_lamg) {
        // init pointers to data
        double* lamfdest = lamg;
        double* lamcdest = lamg + nx;
        const double* lamfsrc = lam_f0;
        const double* lamcsrc = lam_c0;
        // iterate over mesh intervals (k from 0 to N-2)
        for (int k = 0; k < N-1; ++k) {
            // copy lamf0[:,k] and lamc0[:,k]
            if (lam_f0 && flag_lamf) memcpy(lamfdest, lamfsrc, nx*sizeof(double));
            if (lam_c0 && flag_lamc) memcpy(lamcdest, lamcsrc, nc*sizeof(double)); 
            // move pointers forwards
            lamfdest += nx+nc;
            lamcdest += nx+nc;
            lamfsrc += nx;
            lamcsrc += nc;
        }
        // copy lamc0[:,end]
        lamcdest = lamg + (N-1)*(nx+nc);
        if (lam_c0  && flag_lamc) memcpy(lamcdest, lamcsrc, nc*sizeof(double)); 
        // copy lamb0
        double* lambdest = lamg + (N-1)*(nx+nc) + nc;
        const double* lambsrc = lam_b0;
        if (lam_b0  && flag_lamb) memcpy(lambdest, lambsrc, nb*sizeof(double));
        // copy lamq0
        double* lamqdest = lamg + (N-1)*(nx+nc) + nc + nb;
        const double* lamqsrc = lam_q0;
        if (lam_q0  && flag_lamq) memcpy(lamqdest, lamqsrc, nq*sizeof(double));
    }
}

/** Get pattern of lagragian hessian */
int OCPInterface::get_pattern_hess(
    int* irh,
    int* jch
) {
    if (ocp_hessb && ocp_hessi) {
        for (int i = 0; i < nnzh; ++i) {
            jch[i] = kh[i]/nz;
            irh[i] = kh[i]%nz;
        }
    }
    return nnzh;
}

/** Get pattern of lagragian hessian - 64bit int version */
long int OCPInterface::get_pattern_hess(
    long int* irh,
    long int* jch
) {
    if (ocp_hessb && ocp_hessi) {
        for (long int i = 0; i < nnzh; ++i) {
            jch[i] = kh[i]/nz;
            irh[i] = kh[i]%nz;
        }
    }
    return nnzh;
}

/** Get pattern of constrain jacobian */
int OCPInterface::get_pattern_jac(
    int* irj,
    int* jcj
) {
    // iterate over mesh intervals (k from 0 to N-2)
    int c = 0; // nnz counter
    int ir, jc; // temporary row and column indexes
    for (int k = 0; k < N-1; ++k) {
        // add indexes associated with dyn (w.r.t. x1,u,x2,p)
        for (int i = 0; i < nnzdj; ++i) {
            ir = irdj[i]; jc = jcdj[i]; // get current ir,jc
            // convert local ir,jc to global
            ir += k*(nx+nc);
            if (jc < (nx+nu+nx)) jc += k*(nx+nu); //x1,u,x2
            else jc += -(nx+nu+nx) + N*(nx+nu); // p
            // update irj,jcj
            irj[c] = ir;
            jcj[c] = jc;
            // update counter
            ++c;
        }
        // add indexes associated with path (w.r.t. x1,u,p)
        for (int i = 0; i < nnzpj; ++i) {
            ir = irpj[i]; jc = jcpj[i]; // get current ir,jc
            // convert local ir,jc to global
            ir += k*(nx+nc) + nx;
            if (jc < (nx+nu)) jc += k*(nx+nu); //x1,u
            else jc += -(nx+nu) + N*(nx+nu); // p
            // update irj,jcj
            irj[c] = ir;
            jcj[c] = jc;
            // update counter
            ++c;
        }
    }

    // add indexes associated with path for last point
    for (int i = 0; i < nnzpj; ++i) {
        ir = irpj[i]; jc = jcpj[i]; // get current ir,jc
        // convert local ir,jc to global
        ir += (N-1)*(nx+nc);
        if (jc < (nx+nu)) jc += (N-1)*(nx+nu); //x1,u
        else jc += -(nx+nu) + N*(nx+nu); // p
        // update irj,jcj
        irj[c] = ir;
        jcj[c] = jc;
        // update counter
        ++c;
    }

    // add indexes associated with bcs (w.r.t. x0,u0,xn,un)
    for (int i = 0; i < nnzbj; ++i) {
        ir = irbj[i]; jc = jcbj[i]; // get current ir,jc
        // convert local ir,jc to global
        ir += (N-1)*(nx+nc) + nc;
        if (jc >= (nx+nu)) jc += -(nx+nu) + (N-1)*(nx+nu); // xn,un,p
        // update irj,jcj
        irj[c] = ir;
        jcj[c] = jc;
        // update counter
        ++c;
    }

    // add indexes associated with int
    int imax = nnzj-c;
    for (int i = 0; i < imax; ++i) {
        jcj[c] = kj[i]/nz;
        irj[c] = kj[i]%nz;
        ++c;
    }

    // return c, which is the nnz (for internal checks)
    return c;
}

/** Get pattern of constrain jacobian - 64bit int version */
long int OCPInterface::get_pattern_jac(
    long int* irj,
    long int* jcj
) {
    // iterate over mesh intervals (k from 0 to N-2)
    long int c = 0; // nnz counter
    long int ir, jc; // temporary row and column indexes
    for (long int k = 0; k < N-1; ++k) {
        // add indexes associated with dyn (w.r.t. x1,u,x2,p)
        for (long int i = 0; i < nnzdj; ++i) {
            ir = irdj[i]; jc = jcdj[i]; // get current ir,jc
            // convert local ir,jc to global
            ir += k*(nx+nc);
            if (jc < (nx+nu+nx)) jc += k*(nx+nu); //x1,u,x2
            else jc += -(nx+nu+nx) + N*(nx+nu); // p
            // update irj,jcj
            irj[c] = ir;
            jcj[c] = jc;
            // update counter
            ++c;
        }
        // add indexes associated with path (w.r.t. x1,u,p)
        for (long int i = 0; i < nnzpj; ++i) {
            ir = irpj[i]; jc = jcpj[i]; // get current ir,jc
            // convert local ir,jc to global
            ir += k*(nx+nc) + nx;
            if (jc < (nx+nu)) jc += k*(nx+nu); //x1,u
            else jc += -(nx+nu) + N*(nx+nu); // p
            // update irj,jcj
            irj[c] = ir;
            jcj[c] = jc;
            // update counter
            ++c;
        }
    }

    // add indexes associated with path for last point
    for (long int i = 0; i < nnzpj; ++i) {
        ir = irpj[i]; jc = jcpj[i]; // get current ir,jc
        // convert local ir,jc to global
        ir += (N-1)*(nx+nc);
        if (jc < (nx+nu)) jc += (N-1)*(nx+nu); //x1,u
        else jc += -(nx+nu) + N*(nx+nu); // p
        // update irj,jcj
        irj[c] = ir;
        jcj[c] = jc;
        // update counter
        ++c;
    }

    // add indexes associated with bcs (w.r.t. x0,u0,xn,un)
    for (long int i = 0; i < nnzbj; ++i) {
        ir = irbj[i]; jc = jcbj[i]; // get current ir,jc
        // convert local ir,jc to global
        ir += (N-1)*(nx+nc) + nc;
        if (jc >= (nx+nu)) jc += -(nx+nu) + (N-1)*(nx+nu); // xn,un,p
        // update irj,jcj
        irj[c] = ir;
        jcj[c] = jc;
        // update counter
        ++c;
    }

    // add indexes associated with int
    long int imax = nnzj-c;
    for (long int i = 0; i < imax; ++i) {
        jcj[c] = kj[i]/nz;
        irj[c] = kj[i]%nz;
        ++c;
    }

    // return c, which is the nnz (for internal checks)
    return c;
}


int OCPInterface::eval_obj(
    const double* z, 
    double *obj
) {
    // init obj value to 0
    *obj = 0; 

    // exitflag
    int exit;

    // objective result
    double res;
    resrc[0] = &res; // set res

    // init pointers to x1,u,x2,p
    argrc[1] = z + 0; // x1
    argrc[2] = z + nx; // u
    argrc[3] = z + (nx+nu); // x2
    argrc[4] = z + N*(nx+nu); // p
    t = ti; // t

    // iterate over mesh intervals (k from 0 to N-2)
    for (int k = 0; k < N-1; ++k) {
        // update time step
        h = mesh[k] * (tf-ti);
        // call to ocp_runcost
        exit = (*ocp_runcost)(argrc, resrc, iwrc, wrc, memrc);
        // update obj
        *obj += res;
        // update pointers (move forward by nx+nu)
        argrc[1] += nx+nu;
        argrc[2] += nx+nu;
        argrc[3] += nx+nu;
        // update current time
        t += h; 
    }

    // eval boundary cost
    // init pointers to x0,u0,xn,un,p
    argbc[0] = z; // x0
    argbc[1] = z + nx; // u0
    argbc[2] = z + (N-1)*(nx+nu); // xn
    argbc[3] = z + (N-1)*(nx+nu) + nx; // un
    argbc[4] = z + N*(nx+nu); // p
    // set result pointer
    resbc[0] = &res; // set res
    // call to ocp_bcscost
    exit = (*ocp_bcscost)(argbc, resbc, iwbc, wbc, membc);
    // update obj
    *obj += res;
    // return
    return exit;
}

int OCPInterface::eval_obj_grad(
    const double* z, 
    double* grad
) {
    // exitflag
    int exit;

    // init grad to 0
    memset((void*) grad, 0, nz*sizeof(double));

    // init pointers to x1,u,x2,p
    argrcg[1] = z + 0; // x1
    argrcg[2] = z + nx; // u
    argrcg[3] = z + (nx+nu); // x2
    argrcg[4] = z + N*(nx+nu); // p
    t = ti; // t

    // iterate over mesh intervals (k from 0 to N-2)
    for (int k = 0; k < N-1; ++k) {
        // update time step
        h = mesh[k] * (tf-ti);
        // call to ocp_runcost_grad 
        exit = (*ocp_runcost_grad)(argrcg, resrcg, iwrcg, wrcg, memrcg);
        // update gradient
        for (int i = 0; i < nnzrcg; ++i) {
            int idx = krcg[k*nnzrcg + i];
            grad[idx] += resrcg[0][i];
        }
        // update pointers (move forward by nx+nu)
        argrcg[1] += nx+nu;
        argrcg[2] += nx+nu;
        argrcg[3] += nx+nu;
        // update current time
        t += h; 
    }

    // eval boundary cost gradient
    // init pointers to x0,u0,xn,un,p
    argbcg[0] = z; // x0
    argbcg[1] = z + nx; // u0
    argbcg[2] = z + (N-1)*(nx+nu); // xn
    argbcg[3] = z + (N-1)*(nx+nu) + nx; // un
    argbcg[4] = z + N*(nx+nu); // p
    // call to ocp_bcscost_grad 
    exit = (*ocp_bcscost_grad)(argbcg, resbcg, iwbcg, wbcg, membcg);
    // update gradient
    for (int i = 0; i < nnzbcg; ++i) {
        int idx = kbcg[i];
        grad[idx] += resbcg[0][i];
    }
    
    return exit;
}

int OCPInterface::eval_constr(
    const double* z, 
    double* g
){
    
    // exitflag
    int exit;

    // init g to 0
    memset((void*) g, 0, ng*sizeof(double));

    // init pointers to x1,u,x2,p
    argd[1] = z + 0; // x1
    argd[2] = z + nx; // u
    argd[3] = z + (nx+nu); // x2
    argd[4] = z + N*(nx+nu); // p
    argp[1] = z + 0; // x1
    argp[2] = z + nx; // u
    argp[3] = z + N*(nx+nu); // p
    argq[1] = z + 0; // x1
    argq[2] = z + nx; // u
    argq[3] = z + (nx+nu); // x2
    argq[4] = z + N*(nx+nu); // p
    t = ti; // t

    // init pointers to f,c,q
    resd[0] = g; // f
    resp[0] = g + nx; // c
    double *resq0 = g +(N-1)*(nx+nc)+nc+nb; // q

    // iterate over mesh intervals (k from 0 to N-2)
    for (int k = 0; k < N-1; ++k) {
        // update time step
        h = mesh[k] * (tf-ti);
        // call to ocp_dyn
        exit = (*ocp_dyn)(argd, resd, iwd, wd, memd);
        // call to ocp_path
        exit = (*ocp_path)(argp, resp, iwp, wp, memp);
        // call to ocp_int
        exit = (*ocp_int)(argq, resq, iwq, wq, memq);
        for (int i = 0; i < nq; ++i) resq0[i] += resq[0][i];
        // update pointers (move forward)
        argd[1] += nx+nu;
        argd[2] += nx+nu;
        argd[3] += nx+nu;
        argp[1] += nx+nu;
        argp[2] += nx+nu;
        argq[1] += nx+nu;
        argq[2] += nx+nu;
        argq[3] += nx+nu;
        resd[0] += nx+nc;
        resp[0] += nx+nc;
        // update current time
        t += h; 
    }
    // call to ocp_path for last point
    t = tf;
    resp[0] = g + (N-1)*(nx+nc); // set result pointer correctly
    exit = (*ocp_path)(argp, resp, iwp, wp, memp);
    // init pointers to x0,u0,xn,un,p
    argb[0] = z; // x0
    argb[1] = z + nx; // u0
    argb[2] = z + (N-1)*(nx+nu); // xn
    argb[3] = z + (N-1)*(nx+nu) + nx; // un
    argb[4] = z + N*(nx+nu); // p
    // init pointer to b
    resb[0] = g + (N-1)*(nx+nc) + nc; // b
    // call to ocp_bcs
    exit = (*ocp_bcs)(argb, resb, iwb, wb, memb);
 
    return exit;
}

int OCPInterface::eval_constr_jac(
    const double* z, 
    double* jac
){

    // exitflag
    int exit;

    // init jac to 0
    memset((void*) jac, 0, nnzj*sizeof(double));

    // init pointers to x1,u,x2,p
    argdj[1] = z + 0; // x1
    argdj[2] = z + nx; // u
    argdj[3] = z + (nx+nu); // x2
    argdj[4] = z + N*(nx+nu); // p
    argpj[1] = z + 0; // x1
    argpj[2] = z + nx; // u
    argpj[3] = z + N*(nx+nu); // p
    argqj[1] = z + 0; // x1
    argqj[2] = z + nx; // u
    argqj[3] = z + (nx+nu); // x2
    argqj[4] = z + N*(nx+nu); // p
    t = ti; // t

    // init pointers to f,c
    resdj[0] = jac + 0; // w.r.t x1,u,x2,p
    respj[0] = jac + nnzdj; // w.r.t. x1,u,p

    // iterate over mesh intervals (k from 0 to N-2)
    for (int k = 0; k < N-1; ++k) {
        // update time step
        h = mesh[k] * (tf-ti);
        // call to ocp_dyn_jac
        exit = (*ocp_dyn_jac)(argdj, resdj, iwdj, wdj, memdj);
        // call to ocp_path_jac
        exit = (*ocp_path_jac)(argpj, respj, iwpj, wpj, mempj);   
        // call to ocp_int_jac
        exit = (*ocp_int_jac)(argqj, resqj, iwqj, wqj, memqj);   
        for (int i = 0; i < nnzqj; ++i) {
            int idx = kjq[k*nnzqj+i]; // nz index
            if (idx>=0) jac[idx] += resqj[0][i];
        }
        // update pointers (move forward)
        argdj[1] += nx+nu;
        argdj[2] += nx+nu;
        argdj[3] += nx+nu;
        argpj[1] += nx+nu;
        argpj[2] += nx+nu;
        argqj[1] += nx+nu;
        argqj[2] += nx+nu;
        argqj[3] += nx+nu;
        resdj[0] += nnzdj+nnzpj;
        respj[0] += nnzdj+nnzpj;
        // update current time
        t += h; 
    }
    // call to ocp_path_jac for last point
    t = tf;
    respj[0] = jac + (N-1)*(nnzdj+nnzpj); // w.r.t x1,u,p
    exit = (*ocp_path_jac)(argpj, respj, iwpj, wpj, mempj); 
    // init pointers to x0,u0,xn,un,p
    argbj[0] = z; // x0
    argbj[1] = z + nx; // u0
    argbj[2] = z + (N-1)*(nx+nu); // xn
    argbj[3] = z + (N-1)*(nx+nu) + nx; // un
    argbj[4] = z + N*(nx+nu); // p
    // init pointer to b
    resbj[0] = jac + (N-1)*(nnzdj+nnzpj) + nnzpj; // w.r.t x0,u0,xn,un,p
    // call to ocp_bcs_jac
    exit = (*ocp_bcs_jac)(argbj, resbj, iwbj, wbj, membj);

    return exit;
}

int OCPInterface::eval_hessian(
    const double* z, 
    const double sigma, 
    const double* lamg, 
    double* hess
){
    if (ocp_hessb && ocp_hessi) {
        // init hess to 0
        memset((void*) hess, 0, nnzh*sizeof(double));

        // exitflag
        int exit;

        // init pointers to x1,u,x2,p,sigma,lamc,lamf
        arghi[1] = z + 0; // x1
        arghi[2] = z + nx; // u
        arghi[3] = z + (nx+nu); // x2
        arghi[4] = z + N*(nx+nu); // p
        arghi[6] = &sigma; // sigma
        arghi[7] = lamg + nx; // lamc
        arghi[8] = lamg; // lamf
        arghi[9] = lamg + (N-1)*(nx+nc) + nc + nb; // lamq
        t = ti; // t

        // iterate over mesh intervals (k from 0 to N-2)
        for (int k = 0; k < N-1; ++k) {
            // update time step
            h = mesh[k] * (tf-ti);
            // reset reshi
            memset((void*) reshi[0], 0, nnzhi*sizeof(double));
            // call to ocp_hessi
            exit = (*ocp_hessi)(arghi, reshi, iwhi, whi, memhi);
            // assign values
            for (int i = 0; i < nnzhi; ++i) {
                int idx = khi[k*nnzhi + i]; // nz index
                if (idx>=0) hess[idx] += reshi[0][i];
            }
            // update pointers (move forward)
            arghi[1] += nx+nu;
            arghi[2] += nx+nu;
            arghi[3] += nx+nu;
            arghi[7] += nx+nc;
            arghi[8] += nx+nc;
            // update current time
            t += h; 
        }

        /** Eval OCP functions */
        void eval_ocpfuncs(
        );
        arghb[5] = z + N*(nx+nu); // p
        arghb[6] = &sigma; // sigma
        arghb[7] = lamg + (N-1)*(nx+nc); // lamc
        arghb[8] = lamg + (N-1)*(nx+nc) + nc; // lamb
        // reset reshb
        memset((void*) reshb[0], 0, nnzhb*sizeof(double));
        // call to ocp_hessb
        exit = (*ocp_hessb)(arghb, reshb, iwhb, whb, memhb);
        // assign values
        for (int i = 0; i < nnzhb; ++i) {
            int idx = khb[i];
            if (idx>=0) hess[idx] += reshb[0][i];
        }
        return exit;
    } else {
        return 1;
    }
}

/** Eval OCP functions */
void OCPInterface::eval_ocpfuncs() {
    eval_obj_grad(z_opt,  grad_opt);
    eval_constr(z_opt, g_opt);
    eval_constr_jac(z_opt, jac_opt);
    eval_hessian(z_opt, 1, lamg_opt,  hess_opt);
}