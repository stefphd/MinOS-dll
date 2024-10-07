/*  File minos_cinterface.cpp 
    C interface to OCPInterface class
    Copyright (C) 2024 Stefano Lovato
*/

#define MAKE_MINOS

#include "minos.h"
#include "minosc.h"
#include "macros.h"

int MINOSC_PREFIX(new)(
    OCP_t* ocp_ptr, 
    const char* name,
    int N,
    double ti,
    double tf) {
    // create new object and assign pointer
    try {
        *ocp_ptr = (void*) new OCPInterface(name, N, ti, tf);
    } catch (const std::exception& e) {
        *ocp_ptr = NULL;
        return 1;
    }
    return 0;
}

void MINOSC_PREFIX(free)(
    OCP_t* ocp_ptr
) {
    delete CASTOCPINTERFACE(*ocp_ptr); // call destructor
    *ocp_ptr = NULL; // reset to null
}

void MINOSC_PREFIX(get_dims)(
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
) {
    CASTOCPINTERFACE(ocp)->get_dims(nx, nu, np, nc, nb, nq, nz, ng, nnzj, nnzh, na);
}

int MINOSC_PREFIX(get_N)(
    const OCP_t ocp
) {
    return CASTOCPINTERFACE(ocp)->get_N();
}

void MINOSC_PREFIX(set_guess)(
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
) {
    CASTOCPINTERFACE(ocp)->set_guess(x0, u0, p0, lam_x0, lam_u0, lam_p0, lam_f0, lam_c0, lam_b0, lam_q0);
}

void MINOSC_PREFIX(set_bounds)(
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
) {
    CASTOCPINTERFACE(ocp)->set_bounds(lbx, ubx, lbu, ubu, lbp, ubp, lbc, ubc, lbb, ubb, lbq, ubq);
}

int MINOSC_PREFIX(set_option_val)(
    const OCP_t ocp,
    int optkey,
    double val
) {
    return CASTOCPINTERFACE(ocp)->set_option(optkey, val) ? 0 : 1;
}

int MINOSC_PREFIX(set_option)(
    const OCP_t ocp,
    int optkey,
    const char* str
) {
    return CASTOCPINTERFACE(ocp)->set_option(optkey, str) ? 0 : 1;
}

void MINOSC_PREFIX(set_auxdata)(
    const OCP_t ocp,
    double *auxdata
) {
    CASTOCPINTERFACE(ocp)->set_auxdata(auxdata);
}

int MINOSC_PREFIX(solve)(
    const OCP_t ocp
)  {
    return CASTOCPINTERFACE(ocp)->solve();
}

void MINOSC_PREFIX(get_sol)(
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
) {
    CASTOCPINTERFACE(ocp)->get_sol(J_opt, t, 
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
    const OCP_t ocp
) {
    std::string str = CASTOCPINTERFACE(ocp)->toString();
    int len = str.length();
    char* cstr = (char*) malloc((len+1)*sizeof(char));
    std::strcpy(cstr, str.c_str());
    return cstr;    
}

void MINOSC_PREFIX(set_mesh)(
    const OCP_t ocp,
    double *mesh
) {
    CASTOCPINTERFACE(ocp)->set_mesh(mesh);
}

void MINOSC_PREFIX(get_cpu_time)(
    const OCP_t ocp,
    double* tcpu_tot,
    double* tcpu_alg,
    double* tcpu_eval
) {
    CASTOCPINTERFACE(ocp)->get_cpu_time(tcpu_tot, tcpu_alg, tcpu_eval);
}

double MINOSC_PREFIX(get_mu_curr)(
    const OCP_t ocp
) {
    return CASTOCPINTERFACE(ocp)->get_mu_curr();
}

int MINOSC_PREFIX(get_num_iter)(
    const OCP_t ocp
) {
    return CASTOCPINTERFACE(ocp)->get_num_iter();
}

void MINOSC_PREFIX(get_history)(
    const OCP_t ocp,
    double* obj,
    double* inf_pr,
    double* inf_du
) {
    CASTOCPINTERFACE(ocp)->get_history(obj, inf_pr, inf_du);
}

const char* MINOSC_PREFIX(get_version)(
) {
    return OCPInterface::get_version();
}

void MINOSC_PREFIX(set_printfun)(
    const OCP_t ocp,
    int (*ext_print_funptr)(const char *fmt, ...)
) {
    CASTOCPINTERFACE(ocp)->set_printfun(ext_print_funptr);
}

void MINOSC_PREFIX(set_interruptfun)(
    const OCP_t ocp,
    bool (*ext_int_funptr)(void)
) {
    CASTOCPINTERFACE(ocp)->set_interruptfun(ext_int_funptr);
}