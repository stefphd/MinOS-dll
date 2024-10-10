from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp cimport bool

cdef extern from "minos.h":
    cdef cppclass OCPInterface:
        # Constructor
        OCPInterface(string name, int N, double ti, double tf) except +
        # Methods
        void get_dims(int* nx, int* nu, int* np, int* nc, int* nb, int* nq, int* nz, int* ng, int* nnzj, int* nnzh, int* na)
        int get_N()
        void set_guess(double* x0, double* u0, double* p0, double* lam_x0, double* lam_u0, double* lam_p0, double* lam_f0, double* lam_c0, double* lam_b0, double* lam_q0)
        void set_bounds(double* lbx, double* ubx, double* lbu, double* ubu, double* lbp, double* ubp, double* lbc, double* ubc, double* lbb, double* ubb, double* lbq, double* ubq)
        void set_auxdata(double* auxdata)
        void set_mesh(double* mesh)
        int solve() except +
        void get_sol(double *J_opt, double *t, double *x_opt, double *u_opt, double *p_opt, double *lamx_opt, double *lamu_opt, double *lamp_opt, double *lamf_opt, double *lamc_opt, double *lamb_opt, double* lam_qopt, double *f_opt, double *c_opt, double *b_opt, double* q_opt, double *l_opt, double *m_opt, double *gradx_opt, double *gradu_opt, double *gradp_opt, size_t *irj, size_t *jcj, double *jac, size_t *irh, size_t *jch, double *hess)
        void get_mesh(double* mesh)
        string toString()
        void get_cpu_time(double* tcpu_tot, double* tcpu_alg, double* tcpu_eval)
        double get_mu_curr()
        int get_num_iter()
        void get_history(double* obj_history, double* infpr_history, double* infdu_history)
        bool set_option(int opt_id, double val)
        bool set_option(int opt_id, string str)