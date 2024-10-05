# Stub file for OCPInterface
import numpy
from typing import Tuple
import casadi 

class Builder:
    """
    A class to generate and build the OCP Python module embedding the solver.
    Usage is 
        - b = minos.Builder(name: str, outdir: str = "./")
        - b.generate(ocp_runcost: casadi.Function, ocp_bcscost: casadi.Function, ocp_dyn: casadi.Function, ocp_path: casadi.Function, ocp_bcs: casadi.Function, skip_hessian: bool = False)
        - b.build(nlpsolver: str = "ipopt")
    
    Args:
        name (str): name for the OCP Python module to build.
        outdir (str): output directory. Should ends with "/".
    """
    def __init__(self, name: str, outdir: str = "./") -> None:
        """
        Initialize the builder object.

        Args:
            name (str): name for the OCP Python module to build.
            outdir (str): output directory. Should ends with "/".
        """

    def generate(self, ocp_runcost: casadi.Function, ocp_bcscost: casadi.Function, ocp_dyn: casadi.Function, ocp_path: casadi.Function, ocp_bcs: casadi.Function, ocp_int: casadi.Function, skip_hessian: bool = False) -> None:
        """
        Generate the C code for the OCP.

        Args:
            ocp_runcost (casadi.Function): CasADi function of OCP running cost with signature dl=ocp_runcost(t, x1, u, x2, p, h) (dl=ocp_runcost(t, x1, u, x2, p, h, auxdata) in the case of auxdata).
            ocp_bcscost (casadi.Function): CasADi function of OCP boundary cost with signature m=ocp_bcscost(xi, ui, xf, uf, p) (m=ocp_bcscost(xi, ui, xf, uf, p, auxdata) in the case of auxdata).
            ocp_dyn (casadi.Function): CasADi function of OCP dynamic constraints with signature dx=ocp_dyn(t, x, u, p, h) (dx=ocp_dyn(t, x, u, p, h, auxdata) in the case of auxdata).
            ocp_path (casadi.Function): CasADi function of OCP path constraints with signature c=ocp_path(t, x, u, p, h) (dl=ocp_path(t, x, u, p, h, auxdata) in the case of auxdata).
            ocp_bcs (casadi.Function): CasADi function of OCP boundary constraints with signature b=ocp_bcs(xi, ui, xf, uf, p) (m=ocp_bcs(xi, ui, xf, uf, p, auxdata) in the case of auxdata).
            ocp_int (casadi.Function): CasADi function of OCP integral constraints with signature dq=ocp_int(t, x, u, p, h) (dq=ocp_int(t, x, u, p, h, auxdata) in the case of auxdata).
            skip_hessian (bool): use approximated Hessian (skip exact Hessian calculation). This helps in reducing compiling time in the case exact Hessian needs not to be employed.
        
        Returns:
            None
        """
    
    def build(self, nlpsol: str = "ipopt") -> None:
        """
        Build the Python module.

        Args:
            nlpsol (str): NLP solver (either "ipopt" or "knitro").
        
        Returns:
            None
        """

class OCPInterface:
    """
    Class to solve the set and solve the OCP using the built Python module.
    Class methods "set_bounds" and "set_guess" are employed to set the problem bounds and initial guess.
    Solving is done using the class method "solve", while "get_sol" to obtain the solution.

    Args:
        N (int): number of discretization mesh points.
        ti (float): initial time.
        tf (float): final time.
    """

    def __init__(self, N: int, ti: float, tf: float):
        """
        Initialize the OCPInterfac object.
        
        Args:
            N (int): number of discretization mesh points.
            ti (float): initial time.
            tf (float): final time.
        """

    def get_dims(self) -> dict:
        """
        Get the problem dimensions.

        Returns:
            dict: dictionary contaning the problem dimensions (nx, nu, np, nc, nb, nq, nz, ng, nnzj, nnzh, na).
        """

    def get_N(self) -> int:
        """
        Get the number of discretization mesh points.

        Returns
            int: number of discretization points.
        """

    def set_guess(self, x0: numpy.ndarray = None, u0: numpy.ndarray = None, p0: numpy.ndarray = None, lam_x0: numpy.ndarray = None, lam_u0: numpy.ndarray = None, lam_p0: numpy.ndarray = None, lam_f0: numpy.ndarray = None, lam_c0: numpy.ndarray = None, lam_b0: numpy.ndarray = None, lam_q0: numpy.ndarray = None) -> None:
        """
        Set the OCP guess.

        Args:
            x0 (numpy.ndarray): guess for state (size nx-by-N).
            u0 (numpy.ndarray): guess for control (size nu-by-N).
            p0 (numpy.ndarray): guess for parameter (size np).
            lam_x0 (numpy.ndarray): guess for state multiplier (size nx-by-N).
            lam_u0 (numpy.ndarray): guess for control multiplier (size nu-by-N).
            lam_p0 (numpy.ndarray): guess for parameter multiplier (size np).
            lam_f0 (numpy.ndarray): guess for dynamic multiplier (size nx-by-(N-1)).
            lam_c0 (numpy.ndarray): guess for path multiplier (size nc-by-N).
            lam_b0 (numpy.ndarray): guess for boundary multiplier (size nb).
            lam_q0 (numpy.ndarray): guess for integral multiplier (size nq).

        Returns:
            None
        """

    def set_bounds(self, lbx: numpy.ndarray = None, ubx: numpy.ndarray = None, lbu: numpy.ndarray = None, ubu: numpy.ndarray = None, lbp: numpy.ndarray = None, ubp: numpy.ndarray = None, lbc: numpy.ndarray = None, ubc: numpy.ndarray = None, lbb: numpy.ndarray = None, ubb: numpy.ndarray = None, lbq: numpy.ndarray = None, ubq: numpy.ndarray = None) -> None:
        """
        Set the OCP bounds.

        Args:
            lbx (numpy.ndarray): lower bounds for state (size nx).
            ubx (numpy.ndarray): upper bounds for state (size nx).
            lbu (numpy.ndarray): lower bounds for control (size nu).
            ubu (numpy.ndarray): upper bounds for control (size nu).
            lbp (numpy.ndarray): lower bounds for parameters (size np).
            ubp (numpy.ndarray): upper bounds for parameters (size np).
            lbc (numpy.ndarray): lower bounds for path constraints (size nc).
            ubc (numpy.ndarray): upper bounds for path constraints (size nc).
            lbb (numpy.ndarray): lower bounds for boundary conditions (size nb).
            ubb (numpy.ndarray): upper bounds for boundary conditions (size nb).
            lbq (numpy.ndarray): lower bounds for integral constraints (size nq).
            ubq (numpy.ndarray): upper bounds for integral constraints (size nq).

        Returns:
            None
        """

    def set_auxdata(self, auxdata: numpy.ndarray = None) -> None:
        """
        Set the OCP auxdata (if any).

        Args:
            auxdata (numpy.ndarray): auxdata (size na).

        Returns:
            None
        """

    def set_auxdata(self, auxdata: numpy.ndarray = None) -> None:
        """
        Set the mesh fractions.

        Args:
            mesh (numpy.ndarray): mesh fractions that sum to unity (size N-1).

        Returns:
            None
        """

    def solve(self) -> None:
        """
        Solve the OCP.

        Returns:
            None
        """

    def get_sol(self) -> dict:
        """
        Get the OCP solution.

        Args:
            None

        Returns:
            dict: dictionary containing the OCP solution.
        """
    
    def get_mesh(self) -> numpy.ndarray:
        """
        Get the mesh fractions.

        Args:
            None

        Returns:
            numpy.ndarray: mesh fractions that sum to unity (size N-1).

        """

    def get_cpu_time(self) -> Tuple[float, float, float]:
        """
        Get the CPU times.

        Returns:
            Tuple[float, float, float]: A tuple containing:
                - total CPU time.
                - algorithm CPU time.
                - evaluation CPU time.
        """

    def get_mu_curr(self) -> float:
        """
        Get current barrier parameter. This applies only when using interior-point and may be usefull for solving a new problem starting from a previous solution.

        Args:
            None

        Returns:
            float: current barrier parameter.
        """

    def get_num_iter(self) -> int:
        """
        Get the number of iterations.

        Returns:
            int: number of iterations.
        """

    def get_history(self) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
        """
        Get the convergence history.

        Returns:
            Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]: A tuple containining:
                - objective history with size get_num_iter()+1
                - infinite norm of constraint violation history with size get_num_iter()+1
                - infinite norm of optimality history with size get_num_iter()+1
        """

    def set_option(self, **kwargs) -> None:
        """
        Set solver option.
        
        Args:
            **kwargs: keyword-value pair arguments. Possible keys are:
                - max_iter: maximum number of iterations.
                - mu_init: initial barrier parameter.  This applies only when using interior-point and may be usefull for solving a new problem starting from a previous solution.
                - flag_hessian: set to True to employ approximated Hessian (skip exact Hessian evaluation).
                - display: set to False to turn off printing in standard output.
                - print_itersol: print solution file every "print_itersol". 0 to disable.
                - logfile: name of NLP log file. Use "" for default (<NLPSolverName>.log) and "none" for no log file.

        Returns:
            None
        """

    def __str__(self) -> str:
        """
        Return the OCP solution as a string. This may be usefull to print the solution to a file, e.g. with
            with open("output.txt", "w") as file:
                file.write(str(ocp))

        Returns:
            str: string representation.
        """

    @staticmethod
    def get_version() -> str:
        """
        Get the version.

        Returns:
            str: The version string.
        """