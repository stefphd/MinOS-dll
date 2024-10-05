/*  File snoptsolver.h 
    Declaration of NLP solver interface for SNOPT
    Copyright (C) 2024 Stefano Lovato
*/

#ifndef _SNOPTSOLVER_H
#define _SNOPTSOLVER_H

#include "minos.h" // for OCPInterface
#ifdef __cplusplus
extern "C" {
#endif
#include "snopt_cwrap.h" // for SNOPT C interface
#ifdef __cplusplus
}
#endif

class SNOPTSolver {
public:
    /** NLP function 
     * Computes the nonlinear objective and constraint. 
     * F[0] is the objective, F[1:end] are the constraint
     * G[0:nz-1] is the objective gradient, G[nz:end] is the constraint jacobian
    */
    static void nlpfun(
        int *Status, int *n,   double x[],
        int *needF,  int *nF,  double F[],
        int *needG,  int *neG, double G[],
        char   cu[], int   *lencu,
        int    iu[], int   *leniu,
        double ru[], int   *lenru
    );

    /** Iter callback, run at every major iter */
    static void iter_callback(
        int *iAbort, int KTcond[], int *MjrPrt, int *minimz,
        int *m, int *maxS, int *n, int *nb,
        int *nnCon0, int *nnCon, int *nnObj0, int *nnObj, int *nS,
        int *itn, int *nMajor, int *nMinor, int *nSwap,
        double *condHz, int *iObj, double *sclObj, double *ObjAdd,
        double *fObj, double *fMrt, double PenNrm[], double *step,
        double *prInf, double *duInf, double *vimax, double *virel, int hs[],
        int *ne, int *nlocJ, int locJ[], int indJ[], double Jcol[], int *negCon,
        double Ascale[], double bl[], double bu[],
        double Fx[], double fCon[], double gCon[], double gObj[],
        double yCon[], double pi[], double rc[], double rg[], double x[],
        char cu[], int *lencu, int iu[], int *leniu, double ru[], int *lenru,
        char cw[], int *lencw, int iw[], int *leniw, double rw[], int *lenrw 
    );

    /** Finalize the solution  */
    static void finalize_solution(
        int status,
        double* z,
        double* F,
        double* lamz,
        double* lamF,
        OCPInterface* ocpInterface
    );

    /** Solve method (static) */
    static int callSolve(
      OCPInterface *ocpInterface
    );
};

#endif