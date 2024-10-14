/*  File macros.h 
    Macros for MINOS.
    Copyright (C) 2024 Stefano Lovato
*/

#define POSPART(x) ((x>0) ? (+x) : 0) // positive part of x
#define NEGPART(x) ((x<0) ? (-x) : 0) // positive part of x

#define ASSERT_LAMBDA_GUESS ((flag_lamx) && (flag_lamu) && (flag_lamp) && (flag_lamf) && (flag_lamc) && (flag_lamb))

#define PRINTVAL(x) std::cout << #x << ": " << x << std::endl

#define ITERPRINTINT 10 // Print iteration interval
#define HEADERPRINTINT 100 // Print header interval
#define MAX_DEFAULT_ITER 3000 // Default max number of iterations

#define CASTOCPINTERFACE(x) ((OCPInterface*) x) // Cast pointer to OCPInterface class