/*  File main.cpp 
    Main for brachistochrone example
    Copyright (C) 2024 Stefano Lovato
*/

#include <iostream> // for std::cout etc.
#include <fstream> // for write output file
#include <cstdlib> // for atoi
#include "minos.h" // MinOS interface

#define PI 3.14159265358979323846

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
    double tf = 1;

    /* Create OCP Interface */
    OCPInterface* ocp;
    try {
        ocp = new OCPInterface("brachistochrone", N, ti, tf);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    /* Get problem dimensions */
    int nx, nu, np, nc, nb, nq;
    ocp->get_dims(&nx, &nu, &np, &nc, &nb, &nq);

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
    double *x0 = new double[nx*N]; // dim nx*N
    double *u0 = new double[nu*N]; // dim nu*N
    double p0[] = { 1 };
    for (int k = 0; k < N; ++k) {
        x0[3*k] = 2.0*((double) k / (double) (N-1));
        x0[3*k+1] = 2.0*((double) k / (double) (N-1));
        x0[3*k+2] = 0;
        u0[k] = 1.0*((double) k / (double) (N-1));
    }

    /* Set bounds */
    ocp->set_bounds(lbx, ubx,
                   lbu, ubu,
                   lbp, ubp,
                   lbc, ubc,
                   lbb, ubb,
                   lbq, ubq);
    
    /* Set guess */
    ocp->set_guess(x0, u0, p0);

    /* Set NLP solver if any */
    if (!nlpsolver.empty()) ocp->set_option(OCPInterface::NLPSOLVER, nlpsolver);

    /* Call to NLP solver */
    try {
        ocp->solve();
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        delete ocp;
        return 1;
    }
    
    /* Print solution to file */
    std::ofstream outfile;
    outfile.open("brachistochrone.txt");
    outfile << *ocp;
    outfile.close();

    /* Free mem */
    delete ocp;
    delete[] x0; delete[] u0;

    /* Return */
    return 0;
}