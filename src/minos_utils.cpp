/*  File minos_utils.cpp 
    OCPInterface class implementation with utilities
    Copyright (C) 2024 Stefano Lovato
*/

#include "minos.h" 
#include "macros.h" // for nlp macros

/* Check if file exists */
bool OCPInterface::exist_file(
    const char *filename
) {
    std::ifstream f(filename);
    return f.good();
}

/** Check if guess for lambda is given */
bool OCPInterface::check_lambda_guess() {
    if (ASSERT_LAMBDA_GUESS) { return true; }
    return false;
}

/** Get n, m, and nnz of function input/output */
void OCPInterface::get_sizes(
    const int* (*sparsity_fun) (int), 
    const int id,
    int *n,
    int *m,
    int *nnz
) {
    // pointer to sparsity infomation
    const int *sp_i = sparsity_fun(id);
    if (sp_i==0) return; // error

    // number of rows and columns
    int n0 = *sp_i++; // rows
    int m0 = *sp_i++; // columns
    if (n) *n = n0; // rows
    if (m) *m = m0; // columns
    
    // number of nonzeros
    if (nnz) *nnz = sp_i[m0]; 
}

/* Get pattern indexes of function output */
void OCPInterface::get_pattern(
    const int* (*sparsity_fun) (int), 
    const int id,
    int *ir,
    int *jc
) {
    // pointer to sparsity infomation
    const int *sp_i = sparsity_fun(id);
    //if (sp_i==0) return; // error

    // number of rows and columns
    int n = *sp_i++; // rows
    int m = *sp_i++; // columns
    
    // column and row data pointers
    const int *colind = sp_i; // column offsets
    const int *row = sp_i + m + 1; // row nonzero

    // row i, column j indexes of nonzeros
    int c = 0; // init nnz counter
    for(int cc=0; cc < m; ++cc) { // loop over columns
        for(int el = colind[cc]; el < colind[cc+1]; ++el) { // loop over the nonzeros entries of the column
            // save indexes
            if (ir) ir[c] = row[el]; 
            if (jc) jc[c] = cc;
            // update nnz counter
            ++c;
        }
    }
}

/** Convert sparse matric from COO (i.e. tripled) format to CSC (i.e. column compressed) format */
void OCPInterface::coo2csc(int n, int m, int nnz, 
    int *ir1, int *jc1, double *v1, // COO format
    size_t *ir2, size_t *jc2, double *v2  // CSC format
) {
    // Initialize column pointers
    for (int j = 0; j <= m; j++) {
        jc2[j] = 0;
    }
    
    // Count nonzeros in each column
    for (int k = 0; k < nnz; k++) {
        jc2[jc1[k] + 1]++;
    }
    
    // Cumulative sum to get column pointers
    for (int j = 0; j < m; j++) {
        jc2[j + 1] += jc2[j];
    }
    
    // Fill values and row indices
    for (int k = 0; k < nnz; k++) {
        int col = jc1[k];
        int dest = jc2[col];
        ir2[dest] = ir1[k];
        v2[dest] = v1[k];
        jc2[col]++;
    }
    
    // Restore column pointers
    for (int j = m - 1; j >= 0; j--) {
        jc2[j + 1] = jc2[j];
    }
    jc2[0] = 0;
}

int OCPInterface::search_value(
    int arr[], 
    int size, 
    int value
) {
    for (int i = size - 1; i >= 0; --i) {
        if (arr[i] == value) {
            return i; // Return the index if value is found
        }
    }
    return -1; // Return -1 if value is not found
}

/** Init memory */
void OCPInterface::alloc_mem() {
    // init CASADI memory
    alloc_casadi_mem(&ocp_runcost);
    alloc_casadi_mem(&ocp_bcscost);
    alloc_casadi_mem(&ocp_dyn);
    alloc_casadi_mem(&ocp_path);
    alloc_casadi_mem(&ocp_bcs);
    alloc_casadi_mem(&ocp_int);
    alloc_casadi_mem(&ocp_runcost_grad);
    alloc_casadi_mem(&ocp_bcscost_grad);
    alloc_casadi_mem(&ocp_dyn_jac);
    alloc_casadi_mem(&ocp_path_jac);
    alloc_casadi_mem(&ocp_bcs_jac);
    alloc_casadi_mem(&ocp_int_jac);
    if (ocp_hessb.eval && ocp_hessi.eval) {
        alloc_casadi_mem(&ocp_hessb);
        alloc_casadi_mem(&ocp_hessi);
    }
}

void OCPInterface::init_hist_mem(
    int len
) {
    // dealloc old history if any
    if (obj_history)   delete[] obj_history;
    if (infpr_history) delete[] infpr_history;
    if (infdu_history) delete[] infdu_history;
    // alloc history
    obj_history   = new double[len+1]; 
    infpr_history = new double[len+1]; 
    infdu_history = new double[len+1]; 
}

void OCPInterface::init_data_mem() {
    // alloc guess mem
    x0 = new double[1+nx*N];
    u0 = new double[1+nu*N];
    p0 = new double[1+np];
    lam_x0 = new double[1+nx*N];
    lam_u0 = new double[1+nu*N];
    lam_p0 = new double[1+np];
    lam_f0 = new double[1+nx*(N-1)];
    lam_c0 = new double[1+nc*N];
    lam_b0 = new double[1+nb];
    lam_q0 = new double[1+nq*(N-1)];
    // alloc bounds mem
    lbx = new double[1+nx];
    ubx = new double[1+nx];
    lbu = new double[1+nu];
    ubu = new double[1+nu];
    lbp = new double[1+np];
    ubp = new double[1+np];
    lbc = new double[1+nc];
    ubc = new double[1+nc];
    lbb = new double[1+nb];
    ubb = new double[1+nb];
    lbq = new double[1+nq];
    ubq = new double[1+nq];
    //alloc auxdata
    if(na>=0) auxdata = new double[1+na];
    //alloc opt sol
    z_opt = new double[1+nz];
    lamz_opt = new double[1+nz];
    lamg_opt = new double[1+ng];
    g_opt = new double[1+ng];
    grad_opt = new double[1+nz];
    jac_opt = new double[1+nnzj];
    if (nnzh > 0) hess_opt = new double[1+nnzh];
    else hess_opt = NULL;
}

/** Init CASADI memory */
void OCPInterface::alloc_casadi_mem(
    ocpfun_t *ocpfun
) {
    int sz_arg = 0, sz_res = 0, sz_iw = 0, sz_w = 0;
    ocpfun->mem = ocpfun->alloc();
    ocpfun->work(&sz_arg, &sz_res, &sz_iw, &sz_w);
    ocpfun->arg = new const double*[1+sz_arg];
    ocpfun->res = new double*[1+sz_res];
    ocpfun->iw = new int[1+sz_iw];
    ocpfun->w = new double[1+sz_w];
}

/** Delete CASADI memory */
void OCPInterface::dealloc_casadi_mem(
    ocpfun_t *ocpfun
) {
    ocpfun->free(ocpfun->mem);
    if (ocpfun->arg) {
        delete[] ocpfun->arg;
        delete[] ocpfun->res;
        delete[] ocpfun->iw;
        delete[] ocpfun->w;
    }
}

/** Check NLP option file */
bool OCPInterface::check_nlpopt_file(
    const char *filename
) {
    if (OCPInterface::exist_file(filename)) { // check for file exist
        // notify user
        if (display)
            (*print_funptr)("Found option file %s\n", filename);
        return true;
    }
    return false;
}

/** Set the printing function */
void OCPInterface::set_printfun(
    int (*ext_print_funptr)(const char *fmt, ...)
) {
    if (ext_print_funptr) this->print_funptr = ext_print_funptr;
    else this->print_funptr = &printf; // set to default if NULL
}

/** Set the interrupt handling function */
void OCPInterface::set_interruptfun(
    bool (ext_int_funptr)(void)
) {
    if (ext_int_funptr) this->int_funptr = ext_int_funptr;
    else this->int_funptr = &OCPInterface::check_interrupt; // set to default if NULL
}

/** Get mapping from unsorted to sorted indexes */
void OCPInterface::get_sorted_mapping(
    int *arr, 
    int *id, 
    int n
) {
    // Init id to natural indexes
    for (int i = 0; i < n; ++i) id[i] = i;
    // Perform insertion sort on the indices array *id
    for (int i = 1; i < n; ++i) {
        int key_id = id[i];
        int j = i - 1;

        // Move elements of id that are greater than arr[key_id] to one position ahead
        while (j >= 0 && arr[id[j]] > arr[key_id]) {
            id[j + 1] = id[j];
            j = j - 1;
        }
        id[j + 1] = key_id;
    }
}