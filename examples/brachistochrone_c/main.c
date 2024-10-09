/*  File main.c
    Main for brachistochron_c example
    Copyright (C) 2024 Stefano Lovato
*/

#include <stdio.h> // for FILE
#include <stdlib.h> // for malloc and free
#include "minosc.h" // MinOS C interface

#define PI 3.14159265358979323846

/** Entry-point function */
int main(int argc, char* argv[]) {
    /* Check number of input arguments */
    if (argc > 3) {
        fprintf(stderr, "Too many input arguments. Number of input arguments must be 1.\n");
        return 1;
    }

    /* Override nlpsolver if given in input argument */
    const char* nlpsolver = NULL;
    if (argc > 1) {
        nlpsolver = argv[1];
    }
    /* Override N if given in input argument */
    int N = 500;
    if (argc > 2) { 
        N = atoi(argv[2]);
        if (N<=0) {
            fprintf(stderr, "Invalid input arguments. Argument must be a positive number.\n");
            return 1;
        }
    }

    /* Initial and final time */
    double ti = 0;
    double tf = 1;

    /* Create OCP Interface */
    OCP_t ocp;
    minos_new(&ocp, "brachistochrone_c", N, ti, tf);
    if (ocp->exitval) {
        fprintf(stderr, "%s", ocp->exitmsg);
        return 1;
    }
    /*
    // Alternativelly
    ocp = minos_new2("brachistochrone_c", N, ti, tf);
    if (!ocp) {
        fprintf(stderr, "Unable to create OCP instance.\n");
        return 1;
    }
    */

    /* Get problem dimensions */
    int nx, nu, np, nc, nb, nq;
    minos_get_dims(ocp, &nx, &nu, &np, &nc, &nb, &nq, NULL, NULL, NULL, NULL, NULL);

    /* Bounds */
    double lbx[] = { 0, 0, -50 }; // dim nx
    double ubx[] = { 10, 10, 50 }; // dim nx
    double lbu[] = { -PI/2 }; // dim nu
    double ubu[] = { PI/2 }; // dim nu
    double lbp[] = { 0 }; // dim np
    double ubp[] = { 5 }; // dim np
    double *lbc = NULL; // dim nc
    double *ubc = NULL; // dim nc
    double lbb[5] = { 0, 0, 0, 2, 2 }; // dim nb
    double ubb[5] = { 0, 0, 0, 2, 2 }; // dim nb
    double *lbq = NULL; // dim nq
    double *ubq = NULL; // dim nq

    /* Guess */
    double *x0 = (double*) malloc(nx*N*sizeof(double)); // dim nx*N
    double *u0 = (double*) malloc(nu*N*sizeof(double)); // dim nu*N
    double p0[] = { 1 };
    for (int k = 0; k < N; ++k) {
        x0[3*k] = 2.0*((double) k / (double) (N-1));
        x0[3*k+1] = 2.0*((double) k / (double) (N-1));
        x0[3*k+2] = 0;
        u0[k] = 1.0*((double) k / (double) (N-1));
    }

    /* Set bounds */
    minos_set_bounds(ocp, lbx, ubx, lbu, ubu, lbp, ubp, lbc, ubc, lbb, ubb, lbq, ubq);
    
    /* Set guess */
    minos_set_guess(ocp, x0, u0, p0, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    /* Set NLP solver if any */
    if (nlpsolver) minos_set_option(ocp, NLPSOLVER, nlpsolver);

    /* Call to MinOS */
    int status = minos_solve(ocp);
    if (ocp->exitval) {
        fprintf(stderr, "%s", ocp->exitmsg);
        return 1;
    }
    
    /* Print solution to file */
    FILE *filePtr;
    filePtr = fopen("brachistochrone_c.txt", "w");
    const char* str = minos_print_sol(ocp);
    fprintf(filePtr, "%s", str);
    fclose(filePtr);

    /* Free mem */
    minos_free(&ocp);
    free(x0); free(u0);
    
    /* Return */
    return 0;
}