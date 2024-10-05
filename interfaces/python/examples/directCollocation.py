#!/usr/bin/env python

# Direct Collocation problem
import sys
sys.path.append("./../../bin")

import casadi as ca # import CASADI
import minosPy # import MINOS
import matplotlib.pyplot as plt

# Solver options
name = "directCollocation"
skip_hessian = False # use approximated hessian
N = 50
ti = 0
tf = 10

# Create CASADI variables
nx = 2
nu = 1
np = 0
t = ca.SX.sym('t') # independent variable
x = ca.SX.sym('x', nx) # state
u = ca.SX.sym('u', nu) # control
p = ca.SX.sym('p', np) # parameter
xi = ca.SX.sym('xi', nx) # initial state
ui = ca.SX.sym('ui', nu) # initial control
xf = ca.SX.sym('xf', nx) # final state
uf = ca.SX.sym('uf', nu) # final control

# Model equations
xdot = ca.vertcat((1-x[1]**2)*x[0]-x[1]+u[0], x[0]) # state derivative
c = [] # path constraint
q = [] # integral constraint
l = x[0]**2+x[1]**2+u[0]**2 # running cost
m = 0 # boundary cost

# Boundary conditions
b = xi

# Collocation points
d = 3 # degree of interpolating polynomial
tau_root = ca.collocation_points(d, "legendre")
C, D, B = ca.collocation_coeff(tau_root)

# Create augmented state including collocation points
# Augmented state is w = [x, x_1, x_2, ..., x_d], where
# x: state at the start of the interval
# x_1,x_2,...,x_d: states at the collocation points
w = [x] # augmented state
wi = [xi] # augmented initial state
wf = [xf] # augmented final state
# Loop over collocation points
for j in range(0,d):
    xj = ca.SX.sym('x_' + str(j+1), nx)
    xij = ca.SX.sym('xi_' + str(j+1), nx)
    xfj = ca.SX.sym('xf_' + str(j+1), nx)
    w.append(xj)
    wi.append(xij)
    wf.append(xfj)


# Integrator using collocation points
h = ca.SX.sym('h')
f = ca.Function('f', [t, x, u, p], [xdot, l, q])
x1 = ca.SX.sym('x1', nx) # initial state of interval 
x2 = ca.SX.sym('x2', nx) # final state of interval 
# Create auygmented states for x1, x2
w1 = [x1]
w2 = [x2]
for j in range(0,d):
    x1j = ca.SX.sym('x1_' + str(j+1), nx)
    x2j = ca.SX.sym('x2_' + str(j+1), nx)
    w1.append(x1j)
    w2.append(x2j)
# Loop over collocation points
dw = [] # collocation equations
dl = 0 # cost over one mesh interval
dq = 0 # integral over one mesh interval
xend = D[0]*x1 # contribution to the end state due to the start state
for j in range(0,d):
    xp = C[0,j]*w1[0] 
    for r in range(0,d):
        # Expression for the state derivative at the collocation point
        xp += C[r+1,j]*w1[r+1]
    # Append collocation equations
    fj, lj, qj = f(t+tau_root[j]*h, w1[j+1], u, p)
    dw.append(h*fj-xp)
    # Add contributon to the end state
    xend += D[j+1]*w1[j+1]
    # Add contributuon to the quadrature function
    dl += B[j]*lj*h
    dq += B[j]*qj*h
# Continuity constraint of the end state
dw.append(x2-xend)

# CAT list to obtain vectors
w = ca.vertcat(*w)
wi = ca.vertcat(*wi)
wf = ca.vertcat(*wf)
w1 = ca.vertcat(*w1)
w2 = ca.vertcat(*w2)
dw = ca.vertcat(*dw)
nw = w.size1()

# OCP functions
# OCP running cost
ocp_runcost = ca.Function('ocp_runcost', [t, w1, u, w2, p, h], [dl])
# OCP boundary cost
ocp_bcscost = ca.Function('ocp_bcscost', [wi, ui, wf, uf, p], [m])
# OCP dynamic
ocp_dyn = ca.Function('ocp_dyn', [t, w1, u, w2, p, h], [dw])
# OCP path
ocp_path = ca.Function('ocp_path', [t, w, u, p], [c])
# OCP boundary conditions
ocp_bcs = ca.Function('ocp_bcs', [wi, ui, wf, uf, p], [b])
# OCP integral constraints
ocp_int = ca.Function('ocp_int', [t, w1, u, w2, p, h], [dq])

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
lbw = npy.array([-0.25, -1.0e9, -0.25, -1.0e9, -0.25, -1.0e9, -0.25, -1.0e9])
ubw = npy.array([+1.0e9, +1.0e9, +1e9, +1.0e9, +1.0e9, +1.0e9, +1.0e9, +1.0e9])
lbu = npy.array([-1.0])
ubu = npy.array([+1.0])
lbp = None
ubp = None
lbc = None
ubc = None
lbb = npy.array([0.0, 1.0])
ubb = npy.array([0.0, 1.0])
lbq = None
ubq = None
ocp.set_bounds(lbw, ubw, lbu, ubu, lbp, ubp, lbc, ubc, lbb, ubb, lbq, ubq)
# set guess
w0 = npy.vstack((npy.linspace(0, 0, N), 
                 npy.linspace(1, -1, N))*(d+1))
u0 = npy.linspace(0, 0, N)
p0 = npy.array([])
ocp.set_guess(w0,u0,p0)
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

# Plot
plt.figure(1)
plt.plot(sol["t"], sol["x"][0,:], label='x1')
plt.plot(sol["t"], sol["x"][1,:], label='x2')
plt.xlabel("Time")
plt.ylabel("States")
plt.legend()
# plt.savefig("state.png", dpi=300)

plt.figure(2)
plt.plot(sol["t"], sol["u"][0,:], label='u')
plt.xlabel("Time")
plt.ylabel("Control")
plt.legend()
# plt.savefig("control.png", dpi=300)

plt.figure(3)
plt.subplot(2, 1, 1)
plt.plot(iters, obj, label='objective')
plt.xlabel("Iteration")
plt.ylabel("Objective")
plt.subplot(2, 1, 2)
plt.semilogy(iters, infpr, label='Feasiblity')
plt.semilogy(iters, infdu, label='Optimality')
plt.xlabel("Iteration")
plt.ylabel("Max norm")
plt.legend()
# plt.savefig("history.png", dpi=300)

plt.show()
