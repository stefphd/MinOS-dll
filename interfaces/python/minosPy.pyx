#cython: c_string_type=str, c_string_encoding=ascii, language_level=3
from minosPy cimport OCPInterface as CObj
from cython import cast
from libcpp.string cimport string
from distutils import ccompiler
import casadi as ca
import numpy as npy
import os
import sys
import subprocess
import time
import shutil

# Builder class
class Builder:
    # Codegen directories
    __basedir__ = os.path.abspath(os.path.dirname(__file__)) + "/" # folder of this file
    __skip_hessian__ = False # skip hessian

    def __init__(self, name: str, outdir: str = "./") -> None:
        if not outdir.endswith("/"):
            outdir += "/"
        self.name = name
        self.outdir = outdir
        self.cfilename = self.name + ".c"

    def __add2path__(self, path) -> None:
        PATH = os.getenv('PATH') if 'PATH' in os.environ else ""
        if not PATH:
            PATH = ""
        # Add path to the PATH
        if path not in PATH:
            os.environ['PATH'] = path + os.pathsep + PATH
    
    def generate(self, ocp_runcost: ca.Function, ocp_bcscost: ca.Function, ocp_dyn: ca.Function, ocp_path: ca.Function, ocp_bcs: ca.Function, ocp_int: ca.Function, skip_hessian: bool = False) -> None:
        self.__skip_hessian__ = skip_hessian

        # Generate ocp functions
        nx = ocp_runcost.size1_in(1)
        nu = ocp_runcost.size1_in(2)
        np = ocp_runcost.size1_in(4)
        na = -1
        if ocp_runcost.n_in() > 6:
            na = ocp_runcost.size1_in(6)
        # Variables
        x = ca.MX.sym('x', nx)
        x1 = ca.MX.sym('x1', nx)
        x2 = ca.MX.sym('x2', nx)
        u = ca.MX.sym('u', nu); 
        xi = ca.MX.sym('xi', nx)
        ui = ca.MX.sym('ui', nu)
        xf = ca.MX.sym('xf', nx)
        uf = ca.MX.sym('uf', nu)
        p = ca.MX.sym('p', np)
        t = ca.MX.sym('t')
        h = ca.MX.sym('h')
        # Arguments
        args_runcost = [t, x1, u, x2, p, h]
        args_bcscost = [xi, ui, xf, uf, p]
        args_dyn = [t, x1, u, x2, p, h]
        args_path = [t, x, u, p]
        args_bcs = [xi, ui, xf, uf, p]
        args_int = [t, x1, u, x2, p, h]
        # Add auxdata
        if na >= 0:
            auxdata = ca.MX.sym('auxdata', na)
            args_runcost.append(auxdata)
            args_bcscost.append(auxdata)
            args_dyn.append(auxdata)
            args_path.append(auxdata)
            args_bcs.append(auxdata)
            args_int.append(auxdata)
        # Create gradients
        ocp_runcost_grad = ca.Function('ocp_runcost_grad', args_runcost, 
                                    [ca.gradient(ocp_runcost(*args_runcost), ca.vertcat(x1, u, x2, p))])
        ocp_bcscost_grad = ca.Function('ocp_bcscost_grad', args_bcscost, 
                                    [ca.gradient(ocp_bcscost(*args_bcscost), ca.vertcat(xi, ui, xf, uf, p))])
        # Create jacobians
        ocp_dyn_jac = ca.Function('ocp_dyn_jac', args_dyn, 
                                [ca.jacobian(ocp_dyn(*args_dyn), ca.vertcat(x1, u, x2, p))])
        ocp_path_jac = ca.Function('ocp_path_jac', args_path, 
                                [ca.jacobian(ocp_path(*args_path), ca.vertcat(x, u, p))])
        ocp_bcs_jac = ca.Function('ocp_bcs_jac', args_bcs, 
                                [ca.jacobian(ocp_bcs(*args_bcs), ca.vertcat(xi, ui, xf, uf, p))])
        ocp_int_jac = ca.Function('ocp_int_jac', args_int, 
                                [ca.jacobian(ocp_int(*args_int), ca.vertcat(x1, u, x2, p))])
        # Collect all function
        ocp_funcs = [ocp_dyn, ocp_path, ocp_bcs, ocp_int, ocp_runcost, ocp_bcscost, 
                    ocp_dyn_jac, ocp_path_jac, ocp_bcs_jac, ocp_int_jac, ocp_runcost_grad, ocp_bcscost_grad]   
        # Create hessians if not skip_hessian
        if not skip_hessian:
            # num of path and bcs
            nc = ocp_path.size1_out(0)
            nb = ocp_bcs.size1_out(0)
            nq = ocp_int.size1_out(0)
            # create ca variables for multipliers
            sigma = ca.MX.sym('sigma')  # cost multiplier
            lamf = ca.MX.sym('lamf',nx) # dynamic multiplier
            lamc = ca.MX.sym('lamc',nc) # path multiplier
            lamb = ca.MX.sym('lamb',nb) # bcs multiplier
            lamq = ca.MX.sym('lamq',nq) # int multiplier
            # build lagragians
            lagb = sigma*ocp_bcscost(*args_bcscost) # boundary lagragian
            lagi = sigma*ocp_runcost(*args_runcost) # internal lagragian
            if nb > 0: # add bcs lagragian
                lagb += ca.dot(lamb, ocp_bcs(*args_bcs))
            if nx > 0: # add dynamic lagragian
                lagi += ca.dot(lamf, ocp_dyn(*args_dyn))
            if nc > 0: # add path lagragian
                lagb += ca.dot(lamc, ocp_path(args_path[0], xf, uf, args_path[3:]))
                lagi += ca.dot(lamc, ocp_path(args_path[0], x1, u, args_path[3:]))
            if nq > 0: # add integral lagragian
                lagi += ca.dot(lamq, ocp_int(*args_int))
            # generate hessians of lagragians
            hessb, _ = ca.hessian(lagb, ca.vertcat(xi, ui, xf, uf, p))
            hessi, _ = ca.hessian(lagi, ca.vertcat(x1, u, x2, p))
            # create ca functions
            args_hessb = [t, xi, ui, xf, uf, p, sigma, lamc, lamb]
            args_hessi = [t, x1, u, x2, p, h, sigma, lamc, lamf, lamq]
            if na >= 0:
                args_hessb.append(auxdata)
                args_hessi.append(auxdata)
            ocp_hessb = ca.Function('ocp_hessb', args_hessb, [ca.tril(hessb)])
            ocp_hessi = ca.Function('ocp_hessi', args_hessi, [ca.tril(hessi)])
            # append to ocp_funcs
            ocp_funcs.append(ocp_hessb)
            ocp_funcs.append(ocp_hessi)

        # Create outdir if not exists
        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir)
        
        # Generate code
        cg = ca.CodeGenerator(self.cfilename, 
                             {'casadi_int': 'int'})
        for item in ocp_funcs:
            cg.add(item)
        print("Generating C code...")
        cg.generate(self.outdir)

    def build(self) -> None:
        # Check if C file exists
        csource = os.path.join(self.outdir, self.cfilename)
        if not os.path.exists(csource):
            raise FileNotFoundError(f"Unable to find source C file '{csource}'.")
        # Get current PATH
        oldpath = os.getenv('PATH') if 'PATH' in os.environ else ""
        # Select the C compiler
        cc = 'gcc'
        # Test the C compiler
        ccpath = shutil.which(cc)
        if sys.platform == "win32":
            if ccpath: # Global GCC found: give info message
                try:
                    subprocess.run([cc, '--version'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
                except subprocess.CalledProcessError as e:
                    # Throw error
                    raise RuntimeError(f"Unable to find C compiler '{cc}'.\n{e.stderr}")
                print(f"Using user C compiler at '{ccpath}'.")
            else: # Global GCC not found, try local GCC distribution
                self.__add2path__(os.path.join(self.__basedir__, "gcc/bin"))
                # Test GCC again
                ccpath = shutil.which(cc)
        if ccpath is None:
            # Reset default PATH
            if oldpath:
                os.environ['PATH'] = oldpath
            # Throw error
            raise RuntimeError(f"Unable to find C compiler '{cc}'.")       
        try:
            subprocess.run([cc, '--version'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        except subprocess.CalledProcessError as e:
             # Reset default PATH
            if oldpath:
                os.environ['PATH'] = oldpath
            # Throw error
            raise RuntimeError(f"Unable to find C compiler '{cc}'.\n{e.stderr}")
        # Determine the library extension
        if sys.platform == "win32":  # DLL for Windows
            libext = 'dll'
        else:  # SO for other platforms (Linux, macOS, etc.)
            libext = 'so'
        # Define the output library name
        libname = os.path.join(self.outdir, f"{self.name}.{libext}")
        # Build command
        cc_args = f"-shared -O1 -fPIC {csource} -o {libname}"
        cc_cmd = cc + " " + cc_args
        # Run the build process
        start_time = time.time()
        try:
            subprocess.run(cc_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        except subprocess.CalledProcessError as e:
            # Reset default PATH
            if oldpath:
                os.environ['PATH'] = oldpath
            # Throw error
            raise RuntimeError(f"Unable to build library '{libname}' from file '{csource}'.\n{e.stderr}")
        # Print end message
        compile_time = time.time() - start_time
        print(f"Library {libname} built in {compile_time:.2f} seconds.")
        # Reset default PATH
        if oldpath:
            os.environ['PATH'] = oldpath
        # Clean
        os.remove(csource)
        # For Win add <__basedir__> to PATH if not present
        if sys.platform == "win32":
            self.__add2path__(self.__basedir__)

## C++ class
cdef class OCP:
    cdef CObj* cobj # Hold the C++ object

    # Constructor and descrtuctor
    def __init__(self, string name, int N, double ti, double tf):
        self.cobj = new CObj(name, N, ti, tf)
    def __dealloc__(self):
        del self.cobj

    # Methods 
    def get_dims(self) -> dict:
        cdef int nx, nu, np, nc, nb, nq, nz, ng, nnzj, nnzh, na
        self.cobj.get_dims(&nx, &nu, &np, &nc, &nb, &nq, &nz, &ng, &nnzj, &nnzh, &na)
        dims = {"nx": nx,
                "nu": nu,
                "np": np,
                "nc": nc,
                "nb": nb,
                "nq": nq,
                "nz": nz,
                "ng": ng,
                "nnzj": nnzj,
                "nnzh": nnzh,
                "na": na,
                }
        return dims
    
    def get_N(self) -> int:
        cdef int N = self.cobj.get_N()
        return N

    def set_guess(self, x0=None, u0=None, p0=None, lam_x0=None, lam_u0=None, lam_p0=None, lam_f0=None, lam_c0=None, lam_b0=None, lam_q0=None) -> None:
        # Declare vars
        cdef int nx, nu, np, nc, nb, nq, N
        # Get dims
        N = self.cobj.get_N()
        self.cobj.get_dims(&nx, &nu, &np, &nc, &nb, &nq, NULL, NULL, NULL, NULL, NULL)
        # Get memviews, flatten with column-major
        cdef double[:] x0m = x0.flatten(order='F') if x0 is not None else None
        cdef double[:] u0m = u0.flatten(order='F') if u0 is not None else None
        cdef double[:] p0m = p0.flatten(order='F') if p0 is not None else None
        cdef double[:] lam_x0m = lam_x0.flatten(order='F') if lam_x0 is not None else None
        cdef double[:] lam_u0m = lam_u0.flatten(order='F') if lam_u0 is not None else None
        cdef double[:] lam_p0m = lam_p0.flatten(order='F') if lam_p0 is not None else None
        cdef double[:] lam_f0m = lam_f0.flatten(order='F') if lam_f0 is not None else None
        cdef double[:] lam_c0m = lam_c0.flatten(order='F') if lam_c0 is not None else None
        cdef double[:] lam_b0m = lam_b0.flatten(order='F') if lam_b0 is not None else None
        cdef double[:] lam_q0m = lam_q0.flatten(order='F') if lam_q0 is not None else None
        # Check dims
        if x0m is not None and x0m.shape[0]!=nx*N:
            raise ValueError("x0 must have dimension %d" % (nx*N))
        if u0m is not None and u0m.shape[0]!=nu*N:
            raise ValueError("u0 must have dimension %d" % (nu*N))
        if p0 is not None and p0m.shape[0]!=np:
            raise ValueError("p0 must have dimension %d" % (np))
        if lam_x0m is not None and lam_x0m.shape[0]!=nx*N:
            raise ValueError("lam_x0 must have dimension %d" % (nx*N))
        if lam_u0m is not None and lam_u0m.shape[0]!=nu*N:
            raise ValueError("lam_u0 must have dimension %d" % (nu*N))
        if lam_p0m is not None and lam_p0m.shape[0]!=np:
            raise ValueError("lam_p0 must have dimension %d" % (np))
        if lam_f0m is not None and lam_f0m.shape[0]!=nx*(N-1):
            raise ValueError("lam_f0 must have dimension %d" % (nx*(N-1)))
        if lam_c0m is not None and lam_c0m.shape[0]!=nc*N:
            raise ValueError("lam_c0 must have dimension %d" % (nc*N))
        if lam_b0m is not None and lam_b0m.shape[0]!=nb:
            raise ValueError("lam_b0 must have dimension %d" % (nb))
        if lam_q0m is not None and lam_q0m.shape[0]!=nq:
            raise ValueError("lam_q0 must have dimension %d" % (nq))
        # Assign pointers
        cdef double *x0p = &x0m[0] if nx>0 and x0m is not None else NULL
        cdef double *u0p = &u0m[0] if nu>0 and u0m is not None else NULL
        cdef double *p0p = &p0m[0] if np>0 and p0m is not None else NULL
        cdef double *lam_x0p = &lam_x0m[0] if nx>0 and lam_x0m is not None else NULL
        cdef double *lam_u0p = &lam_u0m[0] if nu>0 and lam_u0m is not None else NULL
        cdef double *lam_p0p = &lam_p0m[0] if np>0 and lam_p0m is not None else NULL
        cdef double *lam_f0p = &lam_f0m[0] if nx>0 and lam_f0m is not None else NULL
        cdef double *lam_c0p = &lam_c0m[0] if nc>0 and lam_c0m is not None else NULL
        cdef double *lam_b0p = &lam_b0m[0] if nb>0 and lam_b0m is not None else NULL
        cdef double *lam_q0p = &lam_q0m[0] if nq>0 and lam_q0m is not None else NULL
        # Call to set_guess
        self.cobj.set_guess(x0p, u0p, p0p, lam_x0p, lam_u0p, lam_p0p, lam_f0p, lam_c0p, lam_b0p, lam_q0p)

    def set_bounds(self, double[:] lbx=None, double[:] ubx=None, double[:] lbu=None, double[:] ubu=None, double[:] lbp=None, double[:] ubp=None, double[:] lbc=None, double[:] ubc=None, double[:] lbb=None, double[:] ubb=None, double[:] lbq=None, double[:] ubq=None) -> None:
        # Declare vars
        cdef int nx, nu, np, nc, nb, nq
        # Get dims
        self.cobj.get_dims(&nx, &nu, &np, &nc, &nb, &nq, NULL, NULL, NULL, NULL, NULL)
        # Check dims
        if lbx.shape[0]!=nx:
            raise ValueError("lbx must have dimension %d" % (nx))
        if ubx.shape[0]!=nx:
            raise ValueError("ubx must have dimension %d" % (nx))
        if lbu.shape[0]!=nu:
            raise ValueError("lbu must have dimension %d" % (nu))
        if ubu.shape[0]!=nu:
            raise ValueError("ubu must have dimension %d" % (nu))
        if lbp.shape[0]!=np:
            raise ValueError("lbp must have dimension %d" % (np))
        if ubp.shape[0]!=np:
            raise ValueError("ubp must have dimension %d" % (np))
        if lbc.shape[0]!=nc:
            raise ValueError("lbc must have dimension %d" % (nc))
        if ubc.shape[0]!=nc:
            raise ValueError("ubc must have dimension %d" % (nc))
        if lbb.shape[0]!=nb:
            raise ValueError("lbb must have dimension %d" % (nb))
        if ubb.shape[0]!=nb:
            raise ValueError("ubb must have dimension %d" % (nb))
        if lbq.shape[0]!=nq:
            raise ValueError("lbq must have dimension %d" % (nq))
        if ubq.shape[0]!=nq:
            raise ValueError("ubq must have dimension %d" % (nq))
        # Assign pointers
        cdef double *lbxp = &lbx[0] if nx>0 and lbx is not None else NULL
        cdef double *ubxp = &ubx[0] if nx>0 and ubx is not None else NULL
        cdef double *lbup = &lbu[0] if nu>0 and lbu is not None else NULL
        cdef double *ubup = &ubu[0] if nu>0 and ubu is not None else NULL
        cdef double *lbpp = &lbp[0] if np>0 and lbp is not None else NULL
        cdef double *ubpp = &ubp[0] if np>0 and ubp is not None else NULL
        cdef double *lbcp = &lbc[0] if nc>0 and lbc is not None else NULL
        cdef double *ubcp = &ubc[0] if nc>0 and ubc is not None else NULL
        cdef double *lbbp = &lbb[0] if nb>0 and lbb is not None else NULL
        cdef double *ubbp = &ubb[0] if nb>0 and ubb is not None else NULL
        cdef double *lbqp = &lbq[0] if nq>0 and lbq is not None else NULL
        cdef double *ubqp = &ubq[0] if nq>0 and ubq is not None else NULL
        # Call to set_bounds
        self.cobj.set_bounds(lbxp, ubxp, lbup, ubup, lbpp, ubpp, lbcp, ubcp, lbbp, ubbp, lbqp, ubqp)

    def set_auxdata(self, double[:] auxdata = None) -> None:
        # Declare vars
        cdef int na
        # Get dims
        self.cobj.get_dims(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &na)
        # Check dim 
        if auxdata.shape[0]!=na:
            raise ValueError("auxdata must have dimension %d" % (na))
        # Assign pointers
        cdef double *auxdatap = &auxdata[0] if na>0 and auxdata is not None else NULL
        # Call to set_auxdata
        self.cobj.set_auxdata(auxdatap)

    def set_mesh(self, double[:] mesh = None) -> None:
        # Declare vars
        cdef int N
        cdef double s
        # Get dims 
        N = self.cobj.get_N()
        # Check dim
        if mesh is not None and mesh.shape[0]!=(N-1):
            raise ValueError("mesh must have dimension %d" % (N-1))
        # Check sum to 1
        if mesh is not None and abs(sum(mesh)-1.0)>1.0e-9:
            raise ValueError("mesh must sum to 1")
        # Assign pointers
        cdef double *meshp = &mesh[0] if mesh is not None else NULL
        #Call to set_mesh
        self.cobj.set_mesh(meshp)

    def solve(self) -> None:
        return self.cobj.solve()

    def get_sol(self) -> dict:
        # Declare vars
        cdef int nx, nu, np, nc, nb, nq, N
        # Get dims
        N = self.cobj.get_N()
        self.cobj.get_dims(&nx, &nu, &np, &nc, &nb, &nq, NULL, NULL, NULL, NULL, NULL)
        # Init sol dictionary
        sol = {"objval": 0.0,
                "t": npy.zeros(N, dtype=npy.float64),
                "x": npy.zeros(nx*N, dtype=npy.float64),
                "u": npy.zeros(nu*N, dtype=npy.float64),
                "p": npy.zeros(np, dtype=npy.float64),
                "lamx": npy.zeros(nx*N, dtype=npy.float64),
                "lamu": npy.zeros(nu*N, dtype=npy.float64),
                "lamp": npy.zeros(np, dtype=npy.float64),
                "f": npy.zeros(nx*(N-1), dtype=npy.float64),
                "c": npy.zeros(nc*N, dtype=npy.float64),
                "b": npy.zeros(nb, dtype=npy.float64),
                "q": npy.zeros(nq, dtype=npy.float64),
                "lamf": npy.zeros(nx*(N-1), dtype=npy.float64),
                "lamc": npy.zeros(nc*N, dtype=npy.float64),
                "lamb": npy.zeros(nb, dtype=npy.float64),
                "lamq": npy.zeros(nq, dtype=npy.float64),
                "l": npy.zeros(N-1, dtype=npy.float64),
                "m": 0.0
            }
        # Get pointers
        cdef double J = 0.0
        cdef double[:] t = sol["t"]
        cdef double[:] x = sol["x"]
        cdef double[:] u = sol["u"]
        cdef double[:] p = sol["p"]
        cdef double[:] lamx = sol["lamx"]
        cdef double[:] lamu = sol["lamu"]
        cdef double[:] lamp = sol["lamp"]
        cdef double[:] f = sol["f"]
        cdef double[:] c = sol["c"]
        cdef double[:] b = sol["b"]
        cdef double[:] q = sol["q"]
        cdef double[:] lamf = sol["lamf"]
        cdef double[:] lamc = sol["lamc"]
        cdef double[:] lamb = sol["lamb"]
        cdef double[:] lamq = sol["lamq"]
        cdef double[:] l = sol["l"]
        cdef double m = 0.0
        # Call to get_sol

        self.cobj.get_sol(&J, &t[0], 
                        &x[0] if nx>0 else NULL, &u[0] if nu>0 else NULL, &p[0] if np>0 else NULL,
                        &lamx[0] if nx>0 else NULL, &lamu[0] if nu>0 else NULL, &lamp[0] if np>0 else NULL,
                        &lamf[0] if nx>0 else NULL, &lamc[0] if nc>0 else NULL, &lamb[0] if nb>0 else NULL, &lamq[0] if nq>0 else NULL,
                        &f[0] if nx>0 else NULL, &c[0] if nc>0 else NULL, &b[0] if nb>0 else NULL, &q[0] if nq>0 else NULL, 
                        &l[0], &m, 
                        NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
        # Assign scalars
        sol["objval"] = J
        sol["m"] = m
        # Reshape data with column-major
        sol["x"] = npy.reshape(sol["x"], (nx, N), order='F')
        sol["u"] = npy.reshape(sol["u"], (nu, N), order='F')
        sol["lamx"] = npy.reshape(sol["lamx"], (nx, N), order='F')
        sol["lamu"] = npy.reshape(sol["lamu"], (nu, N), order='F')
        sol["f"] = npy.reshape(sol["f"], (nx, N-1), order='F')
        sol["c"] = npy.reshape(sol["c"], (nc, N-1), order='F')
        sol["lamf"] = npy.reshape(sol["lamf"], (nx, N-1), order='F')
        sol["lamc"] = npy.reshape(sol["lamc"], (nc, N-1), order='F')
        # Return
        return sol

    def get_mesh(self):
        # Declare vars
        cdef int N
        # Get dims
        N = self.cobj.get_N()
        # Init
        mesh = npy.zeros(N-1, dtype=npy.float64)
        # Get pointers
        cdef double[:] meshp = mesh
        # Call to get_mesh
        self.cobj.get_mesh(&meshp[0])
        # Return
        return mesh

    def get_cpu_time(self):
        cdef double tcpu_tot, tcpu_alg, tcpu_eval
        self.cobj.get_cpu_time(&tcpu_tot, &tcpu_alg, &tcpu_eval)
        return tcpu_tot, tcpu_alg, tcpu_eval

    def get_mu_curr(self) -> float:
        return self.cobj.get_mu_curr()
    
    def get_num_iter(self) -> int:
        return self.cobj.get_num_iter()

    def get_history(self):
        cdef int num_iter = self.cobj.get_num_iter()
        cdef Py_ssize_t length = num_iter + 1
        obj_history = npy.zeros(length, dtype=npy.float64)
        infpr_history = npy.zeros(length, dtype=npy.float64)
        infdu_history = npy.zeros(length, dtype=npy.float64)
        cdef double[:] obj_history_mv = obj_history
        cdef double[:] infpr_history_mv = infpr_history
        cdef double[:] infdu_history_mv = infdu_history
        if length>0:
            self.cobj.get_history(&obj_history_mv[0], &infpr_history_mv[0], &infdu_history_mv[0])
        return obj_history, infpr_history, infdu_history

    def set_option(self, **kwargs) -> None:
        for key, value in kwargs.items():
            if key == "max_iter":
                self.cobj.set_option(0, cast(double, value))
            elif key == "mu_init":
                self.cobj.set_option(1, cast(double, value))
            elif key == "flag_hessian":
                self.cobj.set_option(2, cast(double, value))
            elif key == "display":
                self.cobj.set_option(3, cast(double, value))
            elif key == "print_itersol":
                self.cobj.set_option(4, cast(double, value))
            elif key == "logfile":
                self.cobj.set_option(5, cast(string, value))
            elif key == "nlpsolver":
                self.cobj.set_option(6, cast(string, value))
            else:
                raise ValueError(f"Unsupported option key: {key}")

    def __repr__(self):
        return self.cobj.toString()
        
    @staticmethod
    def get_version():
        return CObj.get_version()
