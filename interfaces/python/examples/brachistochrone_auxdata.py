#!/usr/bin/env python

# Brachistochrone problem with auxdata
import sys
sys.path.append("./../../bin")

import casadi as ca # import CASADI
import minosPy # import MINOS
import matplotlib.pyplot as plt

# Solver options
name = "brachistochroneAuxdata"
skip_hessian = False # use approximated hessian
N = 500
ti = 0
tf = 1

# Data
# g = 10 # now this is an auxdata (set later)

# Create CASADI variables
nx = 3
nu = 1
np = 1
na = 1
t = ca.SX.sym('t') # independent variable
x = ca.SX.sym('x', nx) # state
u = ca.SX.sym('u', nu) # control
p = ca.SX.sym('p', np) # parameter
xi = ca.SX.sym('xi', nx) # initial state
ui = ca.SX.sym('ui', nu) # initial control
xf = ca.SX.sym('xf', nx) # final state
uf = ca.SX.sym('uf', nu) # final control
auxdata = ca.SX.sym('auxdata', na) # auxdata

# Model equations ( using time scaled by Tf=p[0] )
g = auxdata[0]
xdot = ca.vertcat(p[0]*x[2]*ca.sin(u[0]), p[0]*x[2]*ca.cos(u[0]), p[0]*g*ca.cos(u[0])) # state derivative
c = [] # path constraint
q = [] # integral constraint
l = 0 # running cost
m = p[0] # boundary cost

# Boundary conditions
b = ca.vertcat(xi[0], xi[1], xi[2], xf[0], xf[1])

# Integrator (trapezoidal)
h = ca.SX.sym('h')
f = ca.Function('f', [t, x, u, p, auxdata], [xdot, l, q])
x1 = ca.SX.sym('x1', nx) # initial state of interval 
x2 = ca.SX.sym('x2', nx) # final state of interval 
f1, dl1, dq1 = f(t,   x1, u, p, auxdata)
f2, dl2, dq2 = f(t+h, x2, u, p, auxdata)
xnext = x1 + 0.5*h*(f1+f2)
dx = xnext - x2
dl = 0.5*h*(dl1+dl2)
dq = 0.5*h*(dq1+dq2)

# OCP functions
# In this case the OCP functions have also 'auxdata' as the last argument.
# OCP running cost
ocp_runcost = ca.Function('ocp_runcost', [t, x1, u, x2, p, h, auxdata], [dl])
# OCP boundary cost
ocp_bcscost = ca.Function('ocp_bcscost', [xi, ui, xf, uf, p, auxdata], [m])
# OCP dynamic
ocp_dyn = ca.Function('ocp_dyn', [t, x1, u, x2, p, h, auxdata], [dx])
# OCP path
ocp_path = ca.Function('ocp_path', [t, x, u, p, auxdata], [c])
# OCP boundary conditions
ocp_bcs = ca.Function('ocp_bcs', [xi, ui, xf, uf, p, auxdata], [b])
# OCP integral constraints
ocp_int = ca.Function('ocp_int', [t, x1, u, x2, p, h, auxdata], [dq])

# Build OCP
builder = minosPy.Builder(name)
builder.generate(ocp_runcost, ocp_bcscost, ocp_dyn, ocp_path, ocp_bcs, ocp_int, skip_hessian)
builder.build()

# Solve OCP
import numpy as npy
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
# solver with different auxdata
auxdata_list = [npy.array([5.0]), 
                npy.array([10.0]), 
                npy.array([15.0])]
sol_list = []
for auxdata in auxdata_list:
    # set auxdata
    ocp.set_auxdata(auxdata)
    # call to solve
    ocp.solve()
    # get sol and append to sol_list
    sol = ocp.get_sol()
    sol_list.append(sol)
    mu_curr = ocp.get_mu_curr()
    # Set new guess
    ocp.set_guess(sol["x"], sol["u"], sol["p"],
                  sol["lamx"], sol["lamu"], sol["lamp"], 
                  sol["lamf"], sol["lamc"], sol["lamb"])
    ocp.set_option(mu_init=mu_curr)

# Print sol
k = 1
for sol in sol_list:
    print("Obj sol %d: %f" % (k, sol["objval"]))
    k += 1
