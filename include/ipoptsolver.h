/*  File ipoptsolver.h 
    Declaration of NLP solver interface for IPOPT
    Copyright (C) 2024 Stefano Lovato
*/

#ifndef _IPOPTSOLVER_H
#define _IPOPTSOLVER_H

#define SET_IPOPT_STR_OPTION(obj, opt, str) obj->Options()->SetStringValue(#opt, str)
#define SET_IPOPT_INT_OPTION(obj, opt, val) obj->Options()->SetIntegerValue(#opt, val)
#define SET_IPOPT_NUM_OPTION(obj, opt, val) obj->Options()->SetNumericValue(#opt, val)

#include <iostream> // for std::cout
#include <ostream> // for std::ostream
#include <cstring> // for memset
#include <math.h> // for log10
#include "minos.h" // for OCPInterface

/*
* IPOPT C++ Interface
* see also https://coin-or.github.io/Ipopt/INTERFACES.html#INTERFACE_CPP
*/
#include <IpTNLP.hpp> // IPOPT C++ interface
#include <IpIpoptApplication.hpp>
#include <IpSolveStatistics.hpp>
#include <IpIpoptData.hpp>
#include <IpTimingStatistics.hpp>
#include <IpTimedTask.hpp>
#include <IpIpoptCalculatedQuantities.hpp>

using namespace Ipopt;

/** IPOPTSolver class for C++ interface with IPOPT */
class IPOPTSolver: public TNLP {

private:

    /**@name Internal variables */
    //@{

    /** OCPInterface pointer */
    OCPInterface *ocpInterface = NULL;

    double *lamzu = NULL, *lamzl = NULL;

    //@}

public:

    /** Default constructor */
    IPOPTSolver(
        OCPInterface *ocpInterface
    );
    
    /** Default destructor */
    virtual ~IPOPTSolver();

    /** Free mem */
    void free_mem();

    /** Solve method (static) */
    static int callSolve(
      OCPInterface *ocpInterface
    );

    /**@name Overloaded from TNLP */
    //@{
    /** Method to return some info about the nlp */
    virtual bool get_nlp_info(
      Index&          n,
      Index&          m,
      Index&          nnz_jac_g,
      Index&          nnz_h_lag,
      IndexStyleEnum& index_style
    );
    
    /** Method to return the bounds for my problem */
    virtual bool get_bounds_info(
      Index   n,
      Number* zl,
      Number* zu,
      Index   m,
      Number* gl,
      Number* gu
    );
    
    /** Method to return the starting point for the algorithm */
    virtual bool get_starting_point(
      Index   n,
      bool    init_z,
      Number* z,
      bool    init_lamz,
      Number* lamzl,
      Number* lamzu,
      Index   m,
      bool    init_lamg,
      Number* lamg
    );
    
    /** Method to return the objective value */
    virtual bool eval_f(
      Index         n,
      const Number* z,
      bool          new_z,
      Number&       obj_value
    );
    
    /** Method to return the gradient of the objective */
    virtual bool eval_grad_f(
      Index         n,
      const Number* z,
      bool          new_z,
      Number*       grad
    );
    
    /** Method to return the constraint residuals */
    virtual bool eval_g(
      Index         n,
      const Number* z,
      bool          new_z,
      Index         m,
      Number*       g
    );
    
    /** Method to return:
    *   1) The structure of the Jacobian (if "values" is NULL)
    *   2) The values of the Jacobian (if "values" is not NULL)
    */
    virtual bool eval_jac_g(
      Index         n,
      const Number* z,
      bool          new_z,
      Index         m,
      Index         nnz,
      Index*        ir,
      Index*        jc,
      Number*       jac
    );
    
    /** Method to return:
    *   1) The structure of the Hessian of the Lagrangian (if "values" is NULL)
    *   2) The values of the Hessian of the Lagrangian (if "values" is not NULL)
    */
    virtual bool eval_h(
      Index         n,
      const Number* z,
      bool          new_z,
      Number        sigma,
      Index         m,
      const Number* lamg,
      bool          new_lamg,
      Index         nnz,
      Index*        ir,
      Index*        jc,
      Number*       hess
    );
    
    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    virtual void finalize_solution(
      SolverReturn               status,
      Index                      n,
      const Number*              z,
      const Number*              lamzl,
      const Number*              lamu,
      Index                      m,
      const Number*              g,
      const Number*              lamg,
      Number                     obj_value,
      const IpoptData*           ip_data,
      IpoptCalculatedQuantities* ip_cq
    );
    
    virtual bool intermediate_callback(
      AlgorithmMode              mode,
      Index                      iter,
      Number                     obj_value,
      Number                     inf_pr,
      Number                     inf_du,
      Number                     mu,
      Number                     d_norm,
      Number                     regularization_size,
      Number                     alpha_du,
      Number                     alpha_pr,
      Index                      ls_trials,
      const IpoptData*           ip_data,
      IpoptCalculatedQuantities* ip_cq
    );
    //@}

private:
    /**@name Methods to block default compiler methods.
    *
    * The compiler automatically generates the following three methods.
    *  Since the default compiler implementation is generally not what
    *  you want (for all but the most simple classes), we usually
    *  put the declarations of these methods in the private section
    *  and never implement them. This prevents the compiler from
    *  implementing an incorrect "default" behavior without us
    *  knowing. (See Scott Meyers book, "Effective C++")
    */
   //@{
   IPOPTSolver(
      const IPOPTSolver&
   );

   IPOPTSolver& operator=(
      const IPOPTSolver&
   );
   //@}

};

#endif