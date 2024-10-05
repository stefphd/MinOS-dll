# Generate goddardRocket problem

from casadi import * # import CASADI
import argparse

# Parse argument
parser = argparse.ArgumentParser()
parser.add_argument('cfilename', type = str, nargs='?', default='goddardRocket.c')
parser.add_argument('outdir', type = str, nargs='?', default='./')
parser.add_argument('--skip_hessian', dest='skip_hessian', action='store_true')
args = parser.parse_args()

# Args
cfilename =  args.cfilename
outdir = args.outdir
skip_hessian = args.skip_hessian

# Boundary conditions
ti = 0
tf = 1
hi = 0
vi = 0
mi = 3
mf = 1

# Data
g0 = 32.174
rho0 = 0.002378
H = 23800
csqrd = 3.264*g0*H
cc  = csqrd**0.5
Tmax = 2*mi*g0
dragk = 0.7110*Tmax/csqrd

# Create CASADI variables
nx = 3
nu = 1
np = 1
t = SX.sym('t') # independent variable
x = SX.sym('x', nx) # state
u = SX.sym('u', nu) # control
p = SX.sym('p', np) # parameter
xi = SX.sym('xi', nx) # initial state
ui = SX.sym('ui', nu) # initial control
xf = SX.sym('xf', nx) # final state
uf = SX.sym('uf', nu) # final control

# Model equations ( using time scaled by Tf=p[0] )
h = x[0]
v = x[1]
m = x[2]
T = u[0]
Tf = p[0]
D = dragk*(v**2)*exp(-h/H)
xdot = vertcat(Tf*v, Tf*((T-D)/m-g0), Tf*(-T/cc)) # state derivative
c = [] # path constraint
q = vertcat(Tf*(-T/cc)) # integral constraint
l = 0 # running cost
m = -xf[0] # boundary cost

# Boundary conditions
b = vertcat(xi[0], xi[1], xi[2], uf[0])

# Integrator (trapezoidal)
h = SX.sym('h')
f = Function('f', [t, x, u, p], [xdot, l, q])
x1 = SX.sym('x1', nx) # initial state of interval 
x2 = SX.sym('x2', nx) # final state of interval 
f1, dl1, dq1 = f(t,   x1, u, p)
f2, dl2, dq2 = f(t+h, x2, u, p)
xnext = x1 + 0.5*h*(f1+f2)
dx = xnext - x2
dl = 0.5*h*(dl1+dl2)
dq = 0.5*h*(dq1+dq2)

# OCP functions
# OCP running cost
ocp_runcost = Function('ocp_runcost', [t, x1, u, x2, p, h], [dl])

# OCP boundary cost
ocp_bcscost = Function('ocp_bcscost', [xi, ui, xf, uf, p], [m])

# OCP dynamic
ocp_dyn = Function('ocp_dyn', [t, x1, u, x2, p, h], [dx])

# OCP path
ocp_path = Function('ocp_path', [t, x, u, p], [c])

# OCP boundary conditions
ocp_bcs = Function('ocp_bcs', [xi, ui, xf, uf, p], [b])

# OCP integral constraints
ocp_int = Function('ocp_int', [t, x1, u, x2, p, h], [dq])

# Create gradients
ocp_runcost_grad = Function('ocp_runcost_grad', [t, x1, u, x2, p, h], 
                            [gradient(ocp_runcost(t, x1, u, x2, p, h), vertcat(x1, u, x2, p))])
ocp_bcscost_grad = Function('ocp_bcscost_grad', [xi, ui, xf, uf, p], 
                            [gradient(ocp_bcscost(xi, ui, xf, uf, p), vertcat(xi, ui, xf, uf, p))])

# Create jacobians
ocp_dyn_jac = Function('ocp_dyn_jac', [t, x1, u, x2, p, h], 
                        [jacobian(ocp_dyn(t, x1, u, x2, p, h), vertcat(x1, u, x2, p))])
ocp_path_jac = Function('ocp_path_jac', [t, x, u, p], 
                        [jacobian(ocp_path(t, x, u, p), vertcat(x, u, p))])
ocp_bcs_jac = Function('ocp_bcs_jac', [xi, ui, xf, uf, p], 
                        [jacobian(ocp_bcs(xi, ui, xf, uf, p), vertcat(xi, ui, xf, uf, p))])
ocp_int_jac = Function('ocp_int_jac', [t, x1, u, x2, p, h], 
                        [jacobian(ocp_int(t, x1, u, x2, p, h), vertcat(x1, u, x2, p))])

# Collect all function
ocp_funcs = [ocp_dyn, ocp_path, ocp_bcs, ocp_int, ocp_runcost, ocp_bcscost, 
             ocp_dyn_jac, ocp_path_jac, ocp_bcs_jac, ocp_int_jac, ocp_runcost_grad, ocp_bcscost_grad]

# Create hessians if not skip_hessian
if not skip_hessian:
    # num of path and bcs
    nc = ocp_path.size1_out(0)
    nb = ocp_bcs.size1_out(0)
    nq = ocp_int.size1_out(0)
    # create CASADI variables for multipliers
    sigma = SX.sym('sigma') # cost multiplier
    lamf = SX.sym('lamf',nx) # dynamic multiplier
    lamc = SX.sym('lamc',nc) # path multiplier
    lamb = SX.sym('lamb',nb) # bcs multiplier
    lamq = SX.sym('lamq',nq) # int multiplier
    # build lagragians
    lagb = sigma*ocp_bcscost(xi, ui, xf, uf, p) # boundary lagragian
    lagi = sigma*ocp_runcost(t, x1, u, x2, p, h) # internal lagragian
    if nb > 0: # add bcs lagragian
        lagb += dot(lamb, ocp_bcs(xi, ui, xf, uf, p))
    if nx > 0: # add dynamic lagragian
        lagi += dot(lamf, ocp_dyn(t, x1, u, x2, p, h))
    if nc > 0: # add path lagragian
        lagb += dot(lamc, ocp_path(t, xf, uf, p))
        lagi += dot(lamc, ocp_path(t, x1, u, p))
    if nq > 0: # add integral laggragian
        lagi += dot(lamq, ocp_int(t, x1, u, x2, p, h))
    # generate hessians of lagragians
    hessb, _ = hessian(lagb, vertcat(xi, ui, xf, uf, p))
    hessi, _ = hessian(lagi, vertcat(x1, u, x2, p))
    # create CASADI functions
    ocp_hessb = Function('ocp_hessb', [t, xi, ui, xf, uf, p, sigma, lamc, lamb], [tril(hessb)])
    ocp_hessi = Function('ocp_hessi', [t, x1, u, x2, p, h, sigma, lamc, lamf, lamq], [tril(hessi)])
    # append to ocp_funcs
    ocp_funcs.append(ocp_hessb)
    ocp_funcs.append(ocp_hessi)

# Call to code generator
if __name__ == '__main__':
    print('Generating %s in %s' % (cfilename, outdir))
    cg = CodeGenerator(cfilename, 
                       {'casadi_int': 'int'})
    for item in ocp_funcs:
        cg.add(item)
    cg.generate(outdir)
    print('Done.')