/*  File knitrosolver.h 
    Declaration of NLP solver interface for KNITRO
    Copyright (C) 2024 Stefano Lovato
*/

#ifndef _KNITROSOLVER_H
#define _KNITROSOLVER_H

#include "minos.h" // for OCPInterface
#include "knitro.h" // for KNITRO application interface

class KNITROSolver {
public:
    /** @name Static methods for KNITRO callbacks */
    //@{
    /** Eval objective and constraints */
    static int callbackEvalFC(
        KN_context_ptr             kc,
        CB_context_ptr             cb,
        KN_eval_request_ptr const  evalRequest,
        KN_eval_result_ptr  const  evalResult,
        void              * const  ocpInterface
    );

    /** Eval gradient and jacobian */
    static int callbackEvalGA(
        KN_context_ptr             kc,
        CB_context_ptr             cb,
        KN_eval_request_ptr const  evalRequest,
        KN_eval_result_ptr  const  evalResult,
        void              * const  ocpInterface
    );
    
    /** Eval lagragian hessian */
    static int callbackEvalH(
        KN_context_ptr             kc,
        CB_context_ptr             cb,
        KN_eval_request_ptr const  evalRequest,
        KN_eval_result_ptr  const  evalResult,
        void              * const  ocpInterface
    );

    /** Callback for new point */
    static int callbackNewPoint(
        KN_context_ptr        kc,
        const double * const  x,
        const double * const  lambda,
        void         *        ocpInterface
    );
    
    /** Finalize the solution */
    static void finalize_solution(
        KN_context_ptr        kc,
        int                   status,
        OCPInterface*         ocpInterface
    );

    //@}

    /** Solve method (static) */
    static int callSolve(
      OCPInterface *ocpInterface
    );

};

#endif