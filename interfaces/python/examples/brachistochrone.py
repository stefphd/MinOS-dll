#!/usr/bin/env python

# Brachistochrone problem
import sys
if sys.platform == "win32":
    sys.path.append("./../../bin")
else:
    sys.path.append("./../../lib")

import casadi as ca # import CASADI
import minosPy # import MINOS
import numpy as npy

# Solver options
name = "brachistochrone"
skip_hessian = False # use approximated hessian
N = 500
ti = 0
tf = 1

# Data
g = 10

# Create CASADI variables
nx = 3
nu = 1
np = 1
t = ca.SX.sym('t') # independent variable
x = ca.SX.sym('x', nx) # state
u = ca.SX.sym('u', nu) # control
p = ca.SX.sym('p', np) # parameter
xi = ca.SX.sym('xi', nx) # initial state
ui = ca.SX.sym('ui', nu) # initial control
xf = ca.SX.sym('xf', nx) # final state
uf = ca.SX.sym('uf', nu) # final control

# Model equations ( using time scaled by Tf=p[0] )
xdot = ca.vertcat(p[0]*x[2]*ca.sin(u[0]), p[0]*x[2]*ca.cos(u[0]), p[0]*g*ca.cos(u[0])) # state derivative
c = [] # path constraint
q = [] # integral constraint
l = 0 # running cost
m = p[0] # boundary cost

# Boundary conditions
b = ca.vertcat(xi[0], xi[1], xi[2], xf[0], xf[1])

# Integrator (trapezoidal)
h = ca.SX.sym('h')
f = ca.Function('f', [t, x, u, p], [xdot, l, q])
x1 = ca.SX.sym('x1', nx) # initial state of interval 
x2 = ca.SX.sym('x2', nx) # final state of interval 
f1, dl1, dq1 = f(t,   x1, u, p)
f2, dl2, dq2 = f(t+h, x2, u, p)
xnext = x1 + 0.5*h*(f1+f2)
dx = xnext - x2
dl = 0.5*h*(dl1+dl2)
dq = 0.5*h*(dq1+dq2)

# OCP functions
# OCP running cost
ocp_runcost = ca.Function('ocp_runcost', [t, x1, u, x2, p, h], [dl])
# OCP boundary cost
ocp_bcscost = ca.Function('ocp_bcscost', [xi, ui, xf, uf, p], [m])
# OCP dynamic
ocp_dyn = ca.Function('ocp_dyn', [t, x1, u, x2, p, h], [dx])
# OCP path
ocp_path = ca.Function('ocp_path', [t, x, u, p], [c])
# OCP boundary conditions
ocp_bcs = ca.Function('ocp_bcs', [xi, ui, xf, uf, p], [b])
# OCP integral constraints
ocp_int = ca.Function('ocp_int', [t, x1, u, x2, p, h], [dq])

# Build OCP
builder = minosPy.Builder(name)
builder.generate(ocp_runcost, ocp_bcscost, ocp_dyn, ocp_path, ocp_bcs, ocp_int, skip_hessian)
builder.build()

# Solve OCP
# create ocp object
ocp = minosPy.OCP(name, N, ti, tf) 
# get dims
dims = ocp.get_dims()
# set bounds
lbx = npy.array([0.0, 0.0, -50.0])
ubx = npy.array([10.0, 10.0, 50.0])
lbu = npy.array([-3.14/2])
ubu = npy.array([+3.14/2])
lbp = npy.array([0.0])
ubp = npy.array([5.0])
lbc = None
ubc = None
lbb = npy.array([0.0, 0.0, 0.0, 2.0, 2.0])
ubb = npy.array([0.0, 0.0, 0.0, 2.0, 2.0])
lbq = None
ubq = None
ocp.set_bounds(lbx, ubx, lbu, ubu, lbp, ubp, lbc, ubc, lbb, ubb, lbq, ubq)
# set guess
x0 = npy.vstack((npy.linspace(0, 2, N), 
                 npy.linspace(0, 2, N), 
                 npy.linspace(0, 0, N)))
u0 = npy.linspace(0, 1, N)
p0 = npy.array([1.0])
ocp.set_guess(x0, u0, p0)
# set mesh
# mesh = 1.0 / (N-1) * npy.ones(N-1)
# ocp.set_mesh(mesh)
# set option
# ocp.set_option(max_iter=30,display=False,nlpsolver="worhp")
# solve
ocp.solve()
# get solution dictionary
sol = ocp.get_sol()
# get CPU times
tcpu_tot, tcpu_alg, tcpu_eval = ocp.get_cpu_time()
# get convergence history
num_iter = ocp.get_num_iter()
obj, infpr, infdu = ocp.get_history()
iters = [i for i in range(num_iter+1)]
# print current barrier parameter
# mu_curr = ocp.get_mu_curr()
# print("Current barrier parameter: %.3e" % mu_curr)
# get mesh
# mesh = ocp.get_mesh()

# Print to file
with open(name + ".txt", "w") as file:
    file.write(str(ocp))
