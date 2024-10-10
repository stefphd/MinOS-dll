/*  File macros.h 
    Macros for MINOS.
    Copyright (C) 2024 Stefano Lovato
*/

#define SPOUT(x) x ## _sparsity_out // get name of CASADI sparsity out function
#define SPIN(x) x ## _sparsity_in // get name of CASADI sparsity in function
#define NIN(x) x ## _n_in // get number of input arguments of CASADI function
#define ALLOC(x) x ## _alloc_mem // get name of CASADI alloc memory function
#define DEALLOC(x) x ## _free_mem // get name of CASADI dealloc memory function
#define WORK(x) x ## _work // get name of CASADI work function
#define POSPART(x) ((x>0) ? (+x) : 0) // positive part of x
#define NEGPART(x) ((x<0) ? (-x) : 0) // positive part of x

#define ASSERT_LAMBDA_GUESS ((flag_lamx) && (flag_lamu) && (flag_lamp) && (flag_lamf) && (flag_lamc) && (flag_lamb))

#define PRINTVAL(x) std::cout << #x << ": " << x << std::endl

#define ITERPRINTINT 10 // Print iteration interval
#define HEADERPRINTINT 100 // Print header interval
#define MAX_DEFAULT_ITER 3000 // Default max number of iterations

#define CASTOCPINTERFACE(x) ((OCPInterface*) x) // Cast pointer to OCPInterface class