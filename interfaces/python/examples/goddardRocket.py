#!/usr/bin/env python

# Goddard-Rocket problem
import sys
sys.path.append("./../../bin")

import casadi as ca # import CASADI
import minosPy # import MINOS
import matplotlib.pyplot as plt

# Solver options
name = "goddardRocket"
skip_hessian = False # use approximated hessian
N = 500
ti = 0
tf = 1

# Data
hi = 0.0
vi = 0.0
mi = 3.0
mf = 1.0
g0 = 32.174
rho0 = 0.002378
H = 23800
hmin = 0.0
hmax = 30.0e3
vmin = -1.5e3
vmax = 1.5e3
mmin = 0.1*mi
mmax = mi

Tfmin = 0.0
Tfmax = 500.0
Tmax = 2.0*mi*g0
csqrd = 3.264*g0*H
cc  = csqrd**0.5
Tmax = 2.0*mi*g0
dragk = 0.7110*Tmax/csqrd

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
h = x[0]
v = x[1]
m = x[2]
T = u[0]
Tf = p[0]
D = dragk*(v**2)*ca.exp(-h/H)
xdot = ca.vertcat(Tf*v, Tf*((T-D)/m-g0), Tf*(-T/cc)) # state derivative
c = [] # path constraint
q = ca.vertcat(Tf*(-T/cc)) # integral constraint
l = 0 # running cost
m = -xf[0] # boundary cost

# Boundary conditions
b = ca.vertcat(xi[0], xi[1], xi[2], uf[0])

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
import numpy as npy
# create ocp object
ocp = minosPy.OCP(name, N, ti, tf) 
# get dims
dims = ocp.get_dims()
# set bounds
lbx = npy.array([hmin, vmin, mmin])
ubx = npy.array([hmax, vmax, mmax])
lbu = npy.array([0.0])
ubu = npy.array([Tmax])
lbp = npy.array([Tfmin])
ubp = npy.array([Tfmax])
lbc = None
ubc = None
lbb = npy.array([hi, vi, mi, 0])
ubb = npy.array([hi, vi, mi, 0])
lbq = npy.array([-(mi-mf)])
ubq = npy.array([-(mi-mf)])
ocp.set_bounds(lbx, ubx, lbu, ubu, lbp, ubp, lbc, ubc, lbb, ubb, lbq, ubq)
# set guess
x0 = npy.vstack((npy.linspace(hi, hmax, N), 
                 npy.linspace(vi, vi, N), 
                 npy.linspace(mi, mf, N)))
u0 = npy.linspace(0, Tmax, N)
p0 = npy.array([Tfmax])
ocp.set_guess(x0,u0,p0)
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
plt.plot(sol["t"], sol["x"][0,:], label='h')
plt.plot(sol["t"], sol["x"][1,:], label='v')
plt.plot(sol["t"], sol["x"][2,:], label='m')
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
