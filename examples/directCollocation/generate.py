# Generate directCollocation problem

from casadi import * # import CASADI
import argparse

# Parse argument
parser = argparse.ArgumentParser()
parser.add_argument('cfilename', type = str, nargs='?', default='directCollocation.c')
parser.add_argument('outdir', type = str, nargs='?', default='./')
parser.add_argument('--skip_hessian', dest='skip_hessian', action='store_true')
args = parser.parse_args()

# Args
cfilename =  args.cfilename
outdir = args.outdir
skip_hessian = args.skip_hessian

# Create CASADI variables
nx = 2
nu = 1
np = 0
t = SX.sym('t') # independent variable
x = SX.sym('x', nx) # state
u = SX.sym('u', nu) # control
p = SX.sym('p', np) # parameter
xi = SX.sym('xi', nx) # initial state
ui = SX.sym('ui', nu) # initial control
xf = SX.sym('xf', nx) # final state
uf = SX.sym('uf', nu) # final control

# Model equations
xdot = vertcat((1-x[1]**2)*x[0]-x[1]+u[0], x[0]) # state derivative
c = [] # path constraint
q = [] # integral constraint
l = x[0]**2+x[1]**2+u[0]**2 # running cost
m = 0 # boundary cost

# Boundary conditions
b = xi

# Collocation points
d = 3 # degree of interpolating polynomial
tau_root = collocation_points(d, "legendre")
C, D, B = collocation_coeff(tau_root)

# Create augmented state including collocation points
# Augmented state is w = [x, x_1, x_2, ..., x_d], where
# x: state at the start of the interval
# x_1,x_2,...,x_d: states at the collocation points
w = [x] # augmented state
wi = [xi] # augmented initial state
wf = [xf] # augmented final state
# Loop over collocation points
for j in range(0,d):
    xj = SX.sym('x_' + str(j+1), nx)
    xij = SX.sym('xi_' + str(j+1), nx)
    xfj = SX.sym('xf_' + str(j+1), nx)
    w.append(xj)
    wi.append(xij)
    wf.append(xfj)

# Integrator using collocation points
h = SX.sym('h')
f = Function('f', [t, x, u, p], [xdot, l, q])
x1 = SX.sym('x1', nx) # initial state of interval 
x2 = SX.sym('x2', nx) # final state of interval 
# Create auygmented states for x1, x2
w1 = [x1]
w2 = [x2]
for j in range(0,d):
    x1j = SX.sym('x1_' + str(j+1), nx)
    x2j = SX.sym('x2_' + str(j+1), nx)
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
w = vertcat(*w)
wi = vertcat(*wi)
wf = vertcat(*wf)
w1 = vertcat(*w1)
w2 = vertcat(*w2)
dw = vertcat(*dw)
nw = w.size1()

# OCP functions
# OCP running cost
ocp_runcost = Function('ocp_runcost', [t, w1, u, w2, p, h], [dl])

# OCP boundary cost
ocp_bcscost = Function('ocp_bcscost', [wi, ui, wf, uf, p], [m])

# OCP dynamic
ocp_dyn = Function('ocp_dyn', [t, w1, u, w2, p, h], [dw])

# OCP path
ocp_path = Function('ocp_path', [t, w, u, p], [c])

# OCP boundary conditions
ocp_bcs = Function('ocp_bcs', [wi, ui, wf, uf, p], [b])

# OCP integral constraints
ocp_int = Function('ocp_int', [t, w1, u, w2, p, h], [dq])

# Create gradients
ocp_runcost_grad = Function('ocp_runcost_grad', [t, w1, u, w2, p, h], 
                            [gradient(ocp_runcost(t, w1, u, w2, p, h), vertcat(w1, u, w2, p))])
ocp_bcscost_grad = Function('ocp_bcscost_grad', [wi, ui, wf, uf, p], 
                            [gradient(ocp_bcscost(wi, ui, wf, uf, p), vertcat(wi, ui, wf, uf, p))])

# Create jacobians
ocp_dyn_jac = Function('ocp_dyn_jac', [t, w1, u, w2, p, h], 
                        [jacobian(ocp_dyn(t, w1, u, w2, p, h), vertcat(w1, u, w2, p))])
ocp_path_jac = Function('ocp_path_jac', [t, w, u, p], 
                        [jacobian(ocp_path(t, w, u, p), vertcat(w, u, p))])
ocp_bcs_jac = Function('ocp_bcs_jac', [wi, ui, wf, uf, p], 
                        [jacobian(ocp_bcs(wi, ui, wf, uf, p), vertcat(wi, ui, wf, uf, p))])
ocp_int_jac = Function('ocp_int_jac', [t, w1, u, w2, p, h], 
                        [jacobian(ocp_int(t, w1, u, w2, p, h), vertcat(w1, u, w2, p))])

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
    lamf = SX.sym('lamf',nw) # dynamic multiplier
    lamc = SX.sym('lamc',nc) # path multiplier
    lamb = SX.sym('lamb',nb) # bcs multiplier
    lamq = SX.sym('lamq',nq) # int multiplier
    # build lagragians
    lagb = sigma*ocp_bcscost(wi, ui, wf, uf, p) # boundary lagragian
    lagi = sigma*ocp_runcost(t, w1, u, w2, p, h) # internal lagragian
    if nb > 0: # add bcs lagragian
        lagb += dot(lamb, ocp_bcs(wi, ui, wf, uf, p))
    if nw > 0: # add dynamic lagragian
        lagi += dot(lamf, ocp_dyn(t, w1, u, w2, p, h))
    if nc > 0: # add path lagragian
        lagb += dot(lamc, ocp_path(t, wf, uf, p))
        lagi += dot(lamc, ocp_path(t, w1, u, p))
    if nq > 0: # add integral laggragian
        lagi += dot(lamq, ocp_int(t, w1, u, w2, p, h))
    # generate hessians of lagragians
    hessb, _ = hessian(lagb, vertcat(wi, ui, wf, uf, p))
    hessi, _ = hessian(lagi, vertcat(w1, u, w2, p))
    # create CASADI functions
    ocp_hessb = Function('ocp_hessb', [t, wi, ui, wf, uf, p, sigma, lamc, lamb], [tril(hessb)])
    ocp_hessi = Function('ocp_hessi', [t, w1, u, w2, p, h, sigma, lamc, lamf, lamq], [tril(hessi)])
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