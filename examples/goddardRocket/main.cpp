/*  File main.cpp 
    Main for Goddard's rocket example
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
    double tf = 1;

    /* Create OCP Interface */
    OCPInterface* ocp;
    try {
        ocp = new OCPInterface("goddardRocket", N, ti, tf);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    /* Get problem dimensions */
    int nx, nu, np, nc, nb, nq;
    ocp->get_dims(&nx, &nu, &np, &nc, &nb, &nq);

    /* Bounds */
    double g0 = 32.174;
    double hi = 0, vi = 0, mi = 3, mf = 1;
    double hmin = 0, hmax = 30e3, vmin = -1.5e3, vmax = 1.5e3, mmin = 0.2*mi, mmax = mi;
    double Tfmin = 0, Tfmax = 500, Tmax = 2*mi*g0;
    double lbx[] = { hmin, vmin, mmin }; // dim nx
    double ubx[] = { hmax, vmax, mmax }; // dim nx
    double lbu[] = { 0 }; // dim nu
    double ubu[] = { Tmax }; // dim nu
    double lbp[] = { Tfmin }; // dim np
    double ubp[] = { Tfmax }; // dim np
    double *lbc = NULL; // dim nc
    double *ubc = NULL; // dim nc
    double lbb[] = { hi, vi, mi, 0 }; // dim nb
    double ubb[] = { hi, vi, mi, 0 }; // dim nb
    double lbq[] = { -(mi-mf) }; // dim nq
    double ubq[] = { -(mi-mf) }; // dim nq

    /* Guess */
    double *x0 = new double[3*N]; // dim nx*N
    double *u0 = new double[1*N]; // dim nu*N
    double p0[] = { Tfmax };
    for (int k = 0; k < N; ++k) {
        x0[3*k] = hi + (hmax-hi)*((double) k / (double) (N-1));
        x0[3*k+1] = vi;
        x0[3*k+2] = mi + (mf-mi)*((double) k / (double) (N-1));
        u0[k] = Tmax*((double) k / (double) (N-1));;
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
    outfile.open("goddardRocket.txt");
    outfile << *ocp;
    outfile.close();

    /* Free mem */
    delete ocp;
    delete[] x0; delete[] u0;

    /* Return */
    return 0;
}