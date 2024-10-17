/*  File minos_cinterface.cpp 
    C interface to OCPInterface class
    Copyright (C) 2024 Stefano Lovato
*/

#include "minos.h"
#include "minosc.h"
#include "macros.h"

void MINOSC_PREFIX(new)(
    OCP_t** ocp_ptr, 
    const char* name,
    int N,
    double ti,
    double tf
) {
    // create new object and assign pointer
    (*ocp_ptr) = new OCP_t;
    (*ocp_ptr)->ptr = NULL;
    (*ocp_ptr)->exitval = 0;
    (*ocp_ptr)->exitmsg = "";
    try {
        (*ocp_ptr)->ptr = (void*) new OCPInterface(name, N, ti, tf);
    } catch (const std::exception& e) {
        (*ocp_ptr)->exitval = 1;
        (*ocp_ptr)->exitmsg = e.what();
    }
}

OCP_t* MINOSC_PREFIX(new2)(
    const char* name,
    int N,
    double ti,
    double tf
) {
    // create new object and assign pointer
    OCP_t* ocp;
    MINOSC_PREFIX(new)(&ocp, name, N, ti, tf);
    return ocp;
}

void MINOSC_PREFIX(free)(
    OCP_t** ocp_ptr
) {
    delete CASTOCPINTERFACE((*ocp_ptr)->ptr); // call destructor
    delete (*ocp_ptr);
}

void MINOSC_PREFIX(get_dims)(
    const OCP_t* ocp, 
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
) {
    CASTOCPINTERFACE(ocp->ptr)->get_dims(nx, nu, np, nc, nb, nq, nz, ng, nnzj, nnzh, na);
}

int MINOSC_PREFIX(get_N)(
    const OCP_t* ocp
) {
    return CASTOCPINTERFACE(ocp->ptr)->get_N();
}

void MINOSC_PREFIX(set_guess)(
    const OCP_t* ocp,
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
) {
    CASTOCPINTERFACE(ocp->ptr)->set_guess(x0, u0, p0, lam_x0, lam_u0, lam_p0, lam_f0, lam_c0, lam_b0, lam_q0);
}

void MINOSC_PREFIX(set_bounds)(
    const OCP_t* ocp,
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
) {
    CASTOCPINTERFACE(ocp->ptr)->set_bounds(lbx, ubx, lbu, ubu, lbp, ubp, lbc, ubc, lbb, ubb, lbq, ubq);
}

int MINOSC_PREFIX(set_option_val)(
    const OCP_t* ocp,
    int optkey,
    double val
) {
    return CASTOCPINTERFACE(ocp->ptr)->set_option(optkey, val) ? 0 : 1;
}

int MINOSC_PREFIX(set_option)(
    const OCP_t* ocp,
    int optkey,
    const char* str
) {
    return CASTOCPINTERFACE(ocp->ptr)->set_option(optkey, str) ? 0 : 1;
}

void MINOSC_PREFIX(set_auxdata)(
    const OCP_t* ocp,
    double *auxdata
) {
    CASTOCPINTERFACE(ocp->ptr)->set_auxdata(auxdata);
}

int MINOSC_PREFIX(solve)(
    OCP_t* ocp
)  {
    ocp->exitval = 0;
    ocp->exitmsg = "";
    int status = -1;
    try {
        status = CASTOCPINTERFACE(ocp->ptr)->solve();
    } catch (const std::exception& e) {
        ocp->exitval = 1;
        ocp->exitmsg = e.what();
        return status;
    }
    return status;
}

void MINOSC_PREFIX(get_sol)(
    const OCP_t* ocp,
    double *J_opt, double *t,
    double *x_opt, double *u_opt, double *p_opt,
    double *lamx_opt, double *lamu_opt, double *lamp_opt,
    double *lamf_opt, double *lamc_opt, double *lamb_opt, double *lamq_opt,
    double *f_opt, double *c_opt, double *b_opt, double *q_opt,
    double *l_opt, double *m_opt,
    double *gradx_opt, double *gradu_opt, double *gradp_opt,
    int *irj, int *jcj, double *jac,
    int *irh, int *jch, double *hess
) {
    CASTOCPINTERFACE(ocp->ptr)->get_sol(J_opt, t, 
                                    x_opt, u_opt, p_opt,
                                    lamx_opt, lamu_opt, lamp_opt,
                                    lamf_opt, lamc_opt, lamb_opt, lamq_opt,
                                    f_opt, c_opt, b_opt, q_opt,
                                    l_opt, m_opt,
                                    gradx_opt, gradu_opt, gradp_opt,
                                    irj, jcj, jac,
                                    irh, jch, hess);
}

const char* MINOSC_PREFIX(print_sol)(
    const OCP_t* ocp
) {
    std::string str = CASTOCPINTERFACE(ocp->ptr)->toString();
    int len = str.length();
    char* cstr = (char*) malloc((len+1)*sizeof(char));
    std::strcpy(cstr, str.c_str());
    return cstr;    
}

void MINOSC_PREFIX(set_mesh)(
    const OCP_t* ocp,
    double *mesh
) {
    CASTOCPINTERFACE(ocp->ptr)->set_mesh(mesh);
}

void MINOSC_PREFIX(get_cpu_time)(
    const OCP_t* ocp,
    double* tcpu_tot,
    double* tcpu_alg,
    double* tcpu_eval
) {
    CASTOCPINTERFACE(ocp->ptr)->get_cpu_time(tcpu_tot, tcpu_alg, tcpu_eval);
}

double MINOSC_PREFIX(get_mu_curr)(
    const OCP_t* ocp
) {
    return CASTOCPINTERFACE(ocp->ptr)->get_mu_curr();
}

int MINOSC_PREFIX(get_num_iter)(
    const OCP_t* ocp
) {
    return CASTOCPINTERFACE(ocp->ptr)->get_num_iter();
}

void MINOSC_PREFIX(get_history)(
    const OCP_t* ocp,
    double* obj,
    double* inf_pr,
    double* inf_du
) {
    CASTOCPINTERFACE(ocp->ptr)->get_history(obj, inf_pr, inf_du);
}

void MINOSC_PREFIX(set_printfun)(
    const OCP_t* ocp,
    int (*ext_print_funptr)(const char *fmt, ...)
) {
    CASTOCPINTERFACE(ocp->ptr)->set_printfun(ext_print_funptr);
}

void MINOSC_PREFIX(set_interruptfun)(
    const OCP_t* ocp,
    bool (*ext_int_funptr)(void)
) {
    CASTOCPINTERFACE(ocp->ptr)->set_interruptfun(ext_int_funptr);
}