/*  File main.cpp 
    Main for direct collocation example
    Copyright (C) 2024 Stefano Lovato
*/

#include <iostream> // for std::cout etc.
#include <fstream> // for write output file
#include <cstdlib> // for atoi
#include "minos.h" // MinOS interface

/** Entry-point function */
int main(int argc, char* argv[]) {
    /* Check number of input arguments */
    if (argc > 3) {
        std::cerr << "Too many input arguments. Number of input arguments must be 1." << std::endl;
        return 1;
    }
    
    /* Override nlpsolver if given in input argument */
    std::string nlpsolver; // default empty
    if (argc > 1) {
        nlpsolver = std::string(argv[1]);
    }
    /* Override N if given in input argument */
    int N = 500;
    if (argc > 2) { 
        N = atoi(argv[2]);
        if (N<=0) {
            std::cerr << "Invalid input arguments. Argument must be a positive number." << std::endl;
            return 1;
        }
    }

    /* Initial and final time */
    double ti = 0;
    double tf = 10;

    /* Create OCP Interface */
    OCPInterface* ocp;
    try {
        ocp = new OCPInterface("directCollocation", N, ti, tf);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    /* Get problem dimensions */
    int nw, nu, np, nc, nb, nq;
    ocp->get_dims(&nw, &nu, &np, &nc, &nb, &nq);

    /* Bounds */
    double lbw[] = { -0.25, -1e9, -0.25, -1e9, -0.25, -1e9, -0.25, -1e9 }; // dim nw
    double ubw[] = { +1e9, +1e9, +1e9, +1e9, +1e9, +1e9, +1e9, +1e9 }; // dim nw
    double lbu[] = { -1 }; // dim nu
    double ubu[] = { +1 }; // dim nu
    double *lbp = NULL; // dim np
    double *ubp = NULL; // dim np
    double *lbc = NULL; // dim nc
    double *ubc = NULL; // dim nc
    double lbb[] = { 0, 1 }; // dim nb
    double ubb[] = { 0, 1 }; // dim nb
    double *lbq = NULL; // dim nq
    double *ubq = NULL; // dim nq

    /* Guess */
    double *w0 = new double[nw*N]; // dim nw*N
    double *u0 = new double[nu*N]; // dim nu*N
    double *p0 = NULL;
    for (int k = 0; k < N; ++k) {
        w0[nw*k] = 0;
        w0[nw*k+1] = +1.0 - 2.0*((double) k / (double) (N-1));
        w0[nw*k+2] = w0[nw*k];
        w0[nw*k+3] = w0[nw*k+1];
        w0[nw*k+4] = w0[nw*k];
        w0[nw*k+5] = w0[nw*k+1];
        w0[nw*k+6] = w0[nw*k];
        w0[nw*k+7] = w0[nw*k+1];
        u0[k] = 0;
    }

    /* Set bounds */
    ocp->set_bounds(lbw, ubw,
                   lbu, ubu,
                   lbp, ubp,
                   lbc, ubc,
                   lbb, ubb,
                   lbq, ubq);
    
    /* Set guess */
    ocp->set_guess(w0, u0, p0);

    /* Set NLP solver if any */
    if (!nlpsolver.empty()) ocp->set_option(OCPInterface::NLPSOLVER, nlpsolver);

    /* Call to MinOS */
    try {
        ocp->solve();
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        delete ocp;
        return 1;
    }
    
    /* Print solution to file */
    std::ofstream outfile;
    outfile.open("directCollocation.txt");
    outfile << *ocp;
    outfile.close();

    /* Free mem */
    delete ocp;
    delete[] w0; delete[] u0;

    /* Return */
    return 0;
}