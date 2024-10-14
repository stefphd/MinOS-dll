/*  File minos_base.cpp 
    OCPInterface class implementation with basic methods
    Copyright (C) 2024 Stefano Lovato
*/

#include "minos.h" 
#include "macros.h" // for nlp macros
#include <signal.h> // for CTRL+C catch using signal

#ifdef WITH_HMACLIC
#include "hmaclic.h" // include for HMACLIC library
#include <stdexcept>  // For std::runtime_error
#ifndef MINOS_PRIVATE_KEY
#define MINOS_PRIVATE_KEY "0000000000"
#endif
#endif

/* Forward declarations */
std::atomic_bool OCPInterface::is_interrupt_requested;

/* Globals to handle interrupt with CTRL+C.*/
void signal_callback_handler(int signum) {
    // Set state of OCPInterface::is_interrupt_requested to true when CTRL+C pressed
    // This stops all instances of OCPInterface, i.e. acts globally for multi-tread applications.
    OCPInterface::is_interrupt_requested = true;
}

/** 
* Iterface class implementation 
*/
OCPInterface::OCPInterface(
    std::string name,
    int N,
    double ti,
    double tf
)  {
#ifdef WITH_HMACLIC
    // Validate license
    char licmsg[HMACLIC_MAXPATH];
    // get machine current hostname and MAC
    char* hostname = get_hostname();
    char* mac = get_mac();
    if (!mac) {
        throw std::runtime_error("Failed to retrieve MAC address for license validation");
        return;
    }
    // look for minos-<hostname>.lic in current directory + specified environment variable
    char lic_filename[HMACLIC_MAXPATH];
    sprintf(lic_filename, "minos-%s.lic", hostname);
    const char* search_envs[] = { "MINOS_LICENSE",
                                "USERPROFILE",
                                "HOME",
                                "PATH"
                                }; // search paths
    char* lic_filename_full = find_lic_file(lic_filename, search_envs, sizeof(search_envs)/sizeof(char*));
    // check if license file found
    if (!lic_filename_full) {
        sprintf(licmsg, "Failed to find license file: %s", lic_filename);
        throw std::runtime_error(licmsg);      
        free(hostname); free(mac);
        return;
    }
    // get license key
    char* license_key, *exp_date;
    if (read_lic_key(lic_filename_full, &license_key, &exp_date)) {
        sprintf(licmsg, "Failed to retrieve license from file: %s", lic_filename_full);
        throw std::runtime_error(licmsg);   
        free(hostname); free(mac);
        free(lic_filename_full);   
        return;
    }
    // validate license key
    int licexit = validate_lic(mac, exp_date, MINOS_PRIVATE_KEY, license_key);
    switch (licexit) {
        case EXIT_EXPIRED:
            sprintf(licmsg, "Expired license from file: %s\nExpiration date: %s", lic_filename_full, exp_date);
            break;
        case EXIT_UNVALID:
            sprintf(licmsg, "Unvalid license key from file: %s", lic_filename_full);
            break;
    }
    free(hostname); free(mac);
    free(lic_filename_full);
    free(license_key); free(exp_date);
    if (licexit != EXIT_VALID) {
        throw std::runtime_error(licmsg);   
        return;
    }
    // OK, valid license
#endif
    // Problem name 
    this->name = name;

    // Init internal variables
    this->N = N; // num of mesh points
    this->ti = ti; // initial time
    this->tf = tf; // final time

    // Load dynamic library <name> with OCP functions
    ocp_lib = load_ocplib(name);
    if (!ocp_lib) { return; }

    // Create default mesh consisting equally-space mesh points
    this->mesh = new double[N-1];
    for (int i = 0; i < N-1; ++i)
        this->mesh[i] = 1.0/( (double) N - 1.0);

    // Get dimensions of OCP
    get_sizes(ocp_path.spin,  1, &nx, NULL, NULL); // num of states
    get_sizes(ocp_path.spin,  2, &nu, NULL, NULL); // num of controls
    get_sizes(ocp_path.spin,  3, &np, NULL, NULL); // num of parameters
    get_sizes(ocp_path.spout, 0, &nc, NULL, NULL); // num of path constr
    get_sizes(ocp_bcs.spout,  0, &nb, NULL, NULL); // num of bcs
    get_sizes(ocp_int.spout,  0, &nq, NULL, NULL); // num of int constr

    // Check for auxdata
    na = -1; // init to -1, i.e. no auxdata
    if ((*ocp_path.nin)() > 4) get_sizes(ocp_path.spin, 4, &na, NULL, NULL); // user ocp_path, but with could also use one of the other functions

    // Calc NLP dims
    nz = N*(nx+nu) + np; // num of NLP vars
    ng = (N-1)*(nx+nc) + nc + nb + nq; // num of NLP constr
    
    // Register signal and signal handler
    // This catches CTRL+C and makes all instances of OCPInterface to stop
    signal(SIGINT, signal_callback_handler);

    // Allocate the memory
    alloc_mem();

    // Init cost gradient
    init_cost_gradient();

    // Init constraint jacobian
    init_constr_jac();

    // Init lagragian hessian
    init_lag_hessian();
    
    // Set constant arg pointers for CASADI functions
    init_args();

    // Init memory to store optimal sol 
    init_data_mem();

    // Set hessian flag
    if (!(ocp_hessb.eval && ocp_hessi.eval))
        this->flag_hessian = true;
    else flag_hessian = false;
}

/** Destructor */
OCPInterface::~OCPInterface() {
}

void OCPInterface::get_dims(
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
    if (nx) *nx = this->nx;
    if (nu) *nu = this->nu;
    if (np) *np = this->np;
    if (nc) *nc = this->nc;
    if (nb) *nb = this->nb;
    if (nq) *nq = this->nq;
    if (nz) *nz = this->nz;
    if (ng) *ng = this->ng;
    if (nnzj) *nnzj = this->nnzj;
    if (nnzh) *nnzh = this->nnzh;
    if (na) *na = (this->na>=0) ? this->na : 0;
}

int OCPInterface::get_N() {
     return this->N;
}

void OCPInterface::set_guess(
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
    // get data: store local copy
    if (x0 && this->x0) memcpy(this->x0, x0, nx*N*sizeof(double));
    if (u0 && this->u0) memcpy(this->u0, u0, nu*N*sizeof(double));
    if (p0 && this->p0) memcpy(this->p0, p0, np*sizeof(double));
    if (lam_x0 && this->lam_x0) memcpy(this->lam_x0, lam_x0, nx*N*sizeof(double));
    if (lam_u0 && this->lam_u0) memcpy(this->lam_u0, lam_u0, nu*N*sizeof(double));
    if (lam_p0 && this->lam_p0) memcpy(this->lam_p0, lam_p0, np*sizeof(double));
    if (lam_f0 && this->lam_f0) memcpy(this->lam_f0, lam_f0, nx*(N-1)*sizeof(double));
    if (lam_c0 && this->lam_c0) memcpy(this->lam_c0, lam_c0, nc*N*sizeof(double));
    if (lam_b0 && this->lam_b0) memcpy(this->lam_b0, lam_b0, nb*sizeof(double));   
    if (lam_q0 && this->lam_q0) memcpy(this->lam_q0, lam_q0, nq*sizeof(double));   
    // set lambda guess flags
    flag_lamx = (bool) lam_x0 || nx==0;
    flag_lamu = (bool) lam_u0 || nu==0;
    flag_lamp = (bool) lam_p0 || np==0;
    flag_lamf = (bool) lam_f0 || nx==0;
    flag_lamc = (bool) lam_c0 || nc==0;
    flag_lamb = (bool) lam_b0 || nb==0;
    flag_lamq = (bool) lam_q0 || nq==0;
}

void OCPInterface::set_optsol_as_guess() {
    // init solution vectors
    std::vector<double> x(nx*N), u(nu*N), p(np),
                        lam_x(nx*N), lam_u(nu*N), lam_p(np),
                        lam_f(nx*(N-1)), lam_c(nc*N), lam_b(nb), lam_q(nq);
    // call get_sol
    get_sol(NULL, NULL,
            x.data(), u.data(), p.data(),
            lam_x.data(), lam_u.data(), lam_p.data(),
            lam_f.data(), lam_c.data(), lam_b.data(), lam_q.data());

    // call set_guess
    set_guess(x.data(), u.data(), p.data(), 
              lam_x.data(), lam_u.data(), lam_p.data(), 
              lam_f.data(), lam_c.data(), lam_b.data(), lam_q.data());
    // set mu_init
    set_option(MU_INIT, mu_curr);
}

void OCPInterface::set_bounds(
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
    // get data: store local copy
    if (lbx && this->lbx) memcpy(this->lbx,  lbx, nx*sizeof(double));
    if (ubx && this->ubx) memcpy(this->ubx,  ubx, nx*sizeof(double));
    if (lbu && this->lbu) memcpy(this->lbu,  lbu, nu*sizeof(double));
    if (ubu && this->ubu) memcpy(this->ubu,  ubu, nu*sizeof(double));
    if (lbp && this->lbp) memcpy(this->lbp,  lbp, np*sizeof(double));
    if (ubp && this->ubp) memcpy(this->ubp,  ubp, np*sizeof(double));
    if (lbc && this->lbc) memcpy(this->lbc,  lbc, nc*sizeof(double));
    if (ubc && this->ubc) memcpy(this->ubc,  ubc, nc*sizeof(double));
    if (lbb && this->lbb) memcpy(this->lbb,  lbb, nb*sizeof(double));
    if (ubb && this->ubb) memcpy(this->ubb,  ubb, nb*sizeof(double));
    if (lbb && this->lbq) memcpy(this->lbq,  lbq, nq*sizeof(double));
    if (ubb && this->ubq) memcpy(this->ubq,  ubq, nq*sizeof(double));
}

void OCPInterface::set_auxdata(
    double* auxdata
) {
    if (auxdata) 
        memcpy(this->auxdata, auxdata, na*sizeof(double));
    if (na >= 0) {
        ocp_runcost.arg[6] = this->auxdata;
        ocp_bcscost.arg[5] = this->auxdata;
        ocp_dyn.arg[6] = this->auxdata;
        ocp_path.arg[4] = this->auxdata;
        ocp_bcs.arg[5] = this->auxdata;
        ocp_int.arg[6] = this->auxdata;
        ocp_runcost_grad.arg[6] = this->auxdata;
        ocp_bcscost_grad.arg[5] = this->auxdata;
        ocp_dyn_jac.arg[6] = this->auxdata;
        ocp_path_jac.arg[4] = this->auxdata;
        ocp_bcs_jac.arg[5] = this->auxdata;
        ocp_int_jac.arg[6] = this->auxdata;
    if (ocp_hessb.eval && ocp_hessi.eval) {
            ocp_hessi.arg[10] = this->auxdata;
            ocp_hessb.arg[9] = this->auxdata;
    }
    }
}

void OCPInterface::set_mesh(
    double *mesh
) {
    if (mesh) memcpy(this->mesh, mesh, (N-1)*sizeof(double));
    else { // equally spaced
        for (int i = 0; i < N-1; ++i)
        this->mesh[i] = 1.0/( (double) N - 1.0);
    }
}

void OCPInterface::get_mesh(
    double *mesh
) {
    if (mesh) memcpy(mesh, this->mesh, (N-1)*sizeof(double));
}

void OCPInterface::get_sol(
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
    
    // objective value
    if (J_opt) *J_opt = this->J_opt;
    
    // time
    if (t) {
        t[0] = ti; // init to ti
        // iterate over mesh intervals (k from 0 to N-2)
        for (int k = 0; k < N-1; ++k)
            t[k+1] = t[k] + mesh[k]*(tf-ti);
    }

    // x, u, p
    if (z_opt) {
        // init pointers to data
        double* xdest = x_opt;
        double* udest = u_opt;
        const double* xsrc = z_opt;
        const double* usrc = z_opt + nx;
        // iterate over mesh points (k from 0 to N-1)
        for (int k = 0; k < N; ++k) {
            // copy x[:,k] and u[:,k]
            if (x_opt) memcpy(xdest, xsrc, nx*sizeof(double));
            if (u_opt) memcpy(udest, usrc, nu*sizeof(double)); 
            // move pointers forwards
            xdest += nx;
            udest += nu;
            xsrc += nx+nu;
            usrc += nx+nu;
        }
        // copy p
        double* pdest = p_opt;
        const double* psrc = z_opt + N*(nx+nu);
        if (p_opt) memcpy(pdest, psrc, np*sizeof(double)); 
    }

    // lamx, lamu, lamp
    if (lamz_opt) {
        // init pointers to data
        double* lamxdest = lamx_opt;
        double* lamudest = lamu_opt;
        const double* lamxsrc = lamz_opt;
        const double* lamusrc = lamz_opt + nx;
        // iterate over mesh points (k from 0 to N-1)
        for (int k = 0; k < N; ++k) {
            // copy lamx[:,k] and lamu[:,k]
            if (lamx_opt) memcpy(lamxdest, lamxsrc, nx*sizeof(double));
            if (lamu_opt) memcpy(lamudest, lamusrc, nu*sizeof(double)); 
            // move pointers forwards
            lamxdest += nx;
            lamudest += nu;
            lamxsrc += nx+nu;
            lamusrc += nx+nu;
        }
        // copy p
        double* lampdest = lamp_opt;
        const double* lampsrc = lamz_opt + N*(nx+nu);
        if (lamp_opt) memcpy(lampdest, lampsrc, np*sizeof(double));
    }

    // lamf, lamc, lamb, lamq
    if (lamg_opt) {
        // init pointers to data
        double* lamfdest = lamf_opt;
        double* lamcdest = lamc_opt;
        const double* lamfsrc = lamg_opt;
        const double* lamcsrc = lamg_opt + nx;
        // iterate over mesh intervals (k from 0 to N-2)
        for (int k = 0; k < N-1; ++k) {
            // copy lamf[:,k] and lamc[:,k]
            if (lamf_opt) memcpy(lamfdest, lamfsrc, nx*sizeof(double));
            if (lamc_opt) memcpy(lamcdest, lamcsrc, nc*sizeof(double)); 
            // move pointers forwards
            lamfdest += nx;
            lamcdest += nc;
            lamfsrc += nx+nc;
            lamcsrc += nx+nc;
        }
        // copy lamc[:,end]
        lamcsrc = lamg_opt + (N-1)*(nx+nc);
        if (lamc_opt) memcpy(lamcdest, lamcsrc, nc*sizeof(double));
        // copy lamb
        double* lambdest = lamb_opt;
        const double* lambsrc = lamg_opt + (N-1)*(nx+nc) + nc;
        if (lamb_opt) memcpy(lambdest, lambsrc, nb*sizeof(double));
        // copy lamq
        double* lamqdest = lamq_opt;
        const double* lamqsrc = lamg_opt + (N-1)*(nx+nc) + nc + nb;
        if (lamq_opt) memcpy(lamqdest, lamqsrc, nq*sizeof(double));
    }

    // f, c, b, q
    if (g_opt) {
        // init pointers to data
        double* fdest = f_opt;
        double* cdest = c_opt;
        const double* fsrc = g_opt;
        const double* csrc = g_opt + nx;
        // iterate over mesh intervals (k from 0 to N-2)
        for (int k = 0; k < N-1; ++k) {
            // copy f[:,k] and c[:,k]
            if (f_opt) memcpy(fdest, fsrc, nx*sizeof(double));
            if (c_opt) memcpy(cdest, csrc, nc*sizeof(double)); 
            // move pointers forwards
            fdest += nx;
            cdest += nc;
            fsrc += nx+nc;
            csrc += nx+nc;
        }
        // copy c[:,end]
        csrc = g_opt + (N-1)*(nx+nc);
        if (c_opt) memcpy(cdest, csrc, nc*sizeof(double));
        // copy b
        double* bdest = b_opt;
        const double* bsrc = g_opt + (N-1)*(nx+nc) + nc;
        if (b_opt) memcpy(bdest, bsrc, nb*sizeof(double));
        // copy q
        double* qdest = q_opt;
        const double* qsrc = g_opt + (N-1)*(nx+nc) + nc + nb;
        if (q_opt) memcpy(qdest, qsrc, nq*sizeof(double));
    }

    // l_opt
    if (z_opt && l_opt) {
        //set arg pointer
        ocp_runcost.arg[1] = z_opt + 0; // x1
        ocp_runcost.arg[2] = z_opt + nx; // u
        ocp_runcost.arg[3] = z_opt + (nx+nu); // x2
        ocp_runcost.arg[4] = z_opt + N*(nx+nu); // p
        this->t = ti; // t
        // set result pointer
        ocp_runcost.res[0] = l_opt;
        for (int k = 0; k < N-1; ++k) {
            // update time step
            this->h = mesh[k] * (tf-ti);
            // call to ocp_runcost
            ocp_runcost.eval(ocp_runcost.arg, ocp_runcost.res, ocp_runcost.iw, ocp_runcost.w, ocp_runcost.mem);
            // update pointers (move forward by nx+nu)
            ocp_runcost.arg[1] += nx+nu;
            ocp_runcost.arg[2] += nx+nu;
            ocp_runcost.arg[3] += nx+nu;
            ocp_runcost.res[0]++;
            // update current time
            this->t += this->h; 
        }
    }

    // m_opt
    if (z_opt && m_opt) {
        //set arg pointer
        ocp_bcscost.arg[0] = z_opt; // x0
        ocp_bcscost.arg[1] = z_opt + nx; // u0
        ocp_bcscost.arg[2] = z_opt + (N-1)*(nx+nu); // xn
        ocp_bcscost.arg[3] = z_opt + (N-1)*(nx+nu) + nx; // un
        ocp_bcscost.arg[4] = z_opt + N*(nx+nu); // p
        // set result pointer
        ocp_bcscost.res[0] = m_opt; // set res
        // call to ocp_bcscost
        ocp_bcscost.eval(ocp_bcscost.arg, ocp_bcscost.res, ocp_bcscost.iw, ocp_bcscost.w, ocp_bcscost.mem);
    }

    // grad_x, grad_u, grad_p
    if (grad_opt) {
        // init pointers to data
        double* gradxdest = gradx_opt;
        double* gradudest = gradu_opt;
        const double* gradxsrc = grad_opt;
        const double* gradusrc = grad_opt + nx;
        // iterate over mesh points (k from 0 to N-1)
        for (int k = 0; k < N; ++k) {
            // copy grad_x[:,k] and grad_u[:,k]
            if (gradx_opt) memcpy(gradxdest, gradxsrc, nx*sizeof(double));
            if (gradu_opt) memcpy(gradudest, gradusrc, nu*sizeof(double)); 
            // move pointers forwards
            gradxdest += nx;
            gradudest += nu;
            gradxsrc += nx+nu;
            gradusrc += nx+nu;
        }
        // copy grad_p
        double* gradpdest = gradp_opt;
        const double* gradpsrc = grad_opt + N*(nx+nu);
        if (gradp_opt) memcpy(gradpdest, gradpsrc, np*sizeof(double)); 
    }

    // return jac 'as is', but using CSC format (instead of COO used here)
    if (irj && jcj && jac) {
        // allocate memory for ir_coo and jc_coo
        int* ir_coo = new int[1+nnzj] { 0 };
        int* jc_coo = new int[1+nnzj] { 0 };
        // get ir_coo and jc_coo
        get_pattern_jac(ir_coo, jc_coo);
        // convert COO to CSC format
        if (jac_opt && ir_coo && jc_coo) coo2csc(ng, nz, nnzj,
                                                 ir_coo, jc_coo, jac_opt,
                                                 irj, jcj, jac);
        // delete
        delete[] ir_coo; delete[] jc_coo;
    }
    // return hess 'as is', but using CSC format (instead of COO used here)
    if (irh && jch && hess) {
        // allocate memory for ir_coo and jc_coo
        int* ir_coo = new int[1+nnzh] { 0 };
        int* jc_coo = new int[1+nnzh] { 0 };
        // get ir_coo and jc_coo
        get_pattern_hess(ir_coo, jc_coo);
        // convert COO to CSC format
        if (hess_opt && ir_coo && jc_coo) coo2csc(nz, nz, nnzh,
                                                 ir_coo, jc_coo, hess_opt,
                                                 irh, jch, hess);
        // delete
        delete[] ir_coo; delete[] jc_coo;
    }
}

bool OCPInterface::set_option(
    int optkey,
    double val
) {
    switch (optkey) {
        case OCPInterface::MAX_ITER:
            if (val<0) { return false; }
            this->max_iter = (int) val;
            break;
        case OCPInterface::MU_INIT:
            if (val<0) { return false; }
            this->mu_init = (double) val;
            break;
        case OCPInterface::FLAG_HESSIAN:
            if ((ocp_hessb.eval && ocp_hessi.eval)) this->flag_hessian = (bool) (val > 0);
            else this->flag_hessian = true;
            break;
        case OCPInterface::DISPLAY:
            this->display = (bool) (val > 0);
            break;
        case OCPInterface::PRINT_ITERSOL:
            if (val<0) { return false; }
            this->print_itersol = (int) val;
            break;
        default: // case not convered
            return false; 
    }
    return true;
}

bool OCPInterface::set_option(
    int optkey,
    std::string str
) {
    switch (optkey) {
        case OCPInterface::LOGFILE:
            logfile = str;
            break;
        case OCPInterface::NLPSOLVER:
            if (str.empty()) nlpsolver = "ipopt"; // set to default
            else nlpsolver = str;
            break;
        default: // case not convered
            return false; 
    }
    return true;
}

void OCPInterface::get_cpu_time(
    double& tcpu_tot,
    double& tcpu_alg,
    double& tcpu_eval
) {
    tcpu_tot = (this->tcpu_tot != NAN) ? this->tcpu_tot : 0;
    tcpu_alg = (this->tcpu_eval != NAN) ? (this->tcpu_tot - this->tcpu_eval) : this->tcpu_tot;
    tcpu_eval = (this->tcpu_eval != NAN) ? this->tcpu_eval : 0;
}
void OCPInterface::get_cpu_time(
    double* tcpu_tot,
    double* tcpu_alg,
    double* tcpu_eval
) {
    this->get_cpu_time(*tcpu_tot, *tcpu_alg, *tcpu_eval);
}

/** Get current value of barrier parameter */
double OCPInterface::get_mu_curr() {
    return mu_curr;
}

/** Get number of iterations */
int OCPInterface::get_num_iter() {
    return num_iter;
}

/** Get objective convergence history */
void OCPInterface::get_history(
    double* obj,
    double* inf_pr,
    double* inf_du
) {
    if (obj && this->obj_history && num_iter>=0) 
        memcpy(obj, this->obj_history, (num_iter+1)*sizeof(double)); 
    if (obj && this->infpr_history && num_iter>=0) 
        memcpy(inf_pr, this->infpr_history, (num_iter+1)*sizeof(double)); 
    if (obj && this->infdu_history && num_iter>=0) 
        memcpy(inf_du, this->infdu_history, (num_iter+1)*sizeof(double)); 
}

/** Check CTRL+C during NLP solution. */
bool OCPInterface::check_interrupt() {
    // do not reset is_interrupt_requested to false in order to stop all instances
    // reset is performed when calling again OCPInterface::solve.
    return OCPInterface::is_interrupt_requested;
}