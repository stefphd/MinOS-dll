/*  File worhpsolver.h 
    Declaration of NLP solver interface for WORHP
    Copyright (C) 2024 Stefano Lovato
*/

#ifndef _WORHPSOLVER_H
#define _WORHPSOLVER_H

#include "minos.h" // for OCPInterface
#include "worhp.h" // for WORHP application interface
#include <memory> // for std::unique_ptr

/** WORHPSolver class for C++ interface with WORHP */
class WORHPSolver {
    private:

    /**@name Internal members */
    //@{
    /** OCPInterface pointer */
    OCPInterface *ocpInterface = NULL;

    friend class OCPInterface;

    /** Worhp data structures */
    OptVar    opt;
    Workspace wsp;
    Params    par;
    Control   cnt;
    int       status;
    TimerType totTimer;
    TimerType evalTimer;
    double    inf_pr, inf_du;

    /** Work arrays */
    int *mapj = NULL; // unsort->sort mapping for jacobian
    int *maph = NULL; // unsort->sort mapping for hessian, -1 for elements not mapped (and thus 0)
    double *jac = NULL, *hess = NULL; // to save temp jacobian and hessian

    /* Init dimensions */
    void init_dim();

    /* Alloc memory */
    void alloc_mem();

    /* Set sparsity */
    void set_sparsity();

    /* Intermediate callback */
    bool int_callback();

    //@}

    public:

    /** Default constructor */
    WORHPSolver(
        OCPInterface *ocpInterface
    );
    
    /** Default destructor */
    virtual ~WORHPSolver();

    /** Solve method (static) */
    static int callSolve(
      OCPInterface *ocpInterface
    );

    /** Solve method */
    int solve();

    /** Finalize the solution */
    void finalize_solution();

    /** Printing function */
    static std::unique_ptr<std::ofstream> logfile;
    static void print_logfile(
        int mode, 
        const char msg[]
    );


};



#endif