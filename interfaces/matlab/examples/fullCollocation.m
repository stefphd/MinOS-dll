%------------------------ Full Collocation Example -------------------------%
% This example can be found in the following reference:                     %
% https://github.com/casadi/casadi/blob/main/docs/examples/matlab/direct_   %
% collocation.m                                                             %
% This employs full (i.e. state & control) collocation instead of state     %
% collocation.                                                              %
% --------------------------------------------------------------------------%

clc, clear

if ispc
    addpath('./../../bin');
else
    addpath('./../../lib');
end

% import casadi
import casadi.*

% Solving options
name = 'fullCollocation'; % name of the problem
skip_hessian = false; % use approximated hessian
N = 50; % number of mesh points
d = 3; % degree of interpolating polynomial
ti = 0; % initial time
tf = 10; % final time

% Boundary conditions
xi = [0; 1]; % initial state x1,x2

% Bounds
lbx = [-0.25; -1e9]; % lower bound for state
ubx = [+1e9; +1e9]; % upper bound for state
lbu = -1; % lower bound for control
ubu = +1; % upper bound for control
lbp = []; % lower bound for parameter
ubp = []; % upper bound for parameter
lbc = []; % lower bound for path constraint 
ubc = []; % upper bound for path constraint
lbb = xi; % lower bound for boundary conditions
ubb = xi; % upper bound for boundary conditions
lbq = []; % lower bound for integral constraints
ubq = []; % upper bound for integral constraints

% Create CASADI functions
nx = 2; nu = 1; np = 0;
t = SX.sym('t'); % independent variable
x = SX.sym('x', nx); % state
u = SX.sym('u', nu); % control
p = SX.sym('p', np); % parameter
x0 = SX.sym('x0', nx); % initial state
u0 = SX.sym('u0', nu); % initial control
xn = SX.sym('xn', nx); % final state
un = SX.sym('un', nu); % final control

% Model equations
xdot = [(1-x(2)^2)*x(1)-x(2)+u(1);
        x(1)]; % state derivative
c = []; % path constraint
q = []; % integral constraint
l = x(1)^2+x(2)^2+u(1)^2; % running (Lagrange) cost
m = 0; % boundary (Mayer) cost

% Boundary condtions
b = x0;

% Get collocations
tau_root = collocation_points(d, 'legendre');
[C, D, B] = collocation_coeff(tau_root);

% Create augmented state and control including collocation points
% Augmented state is w = [x, x_1, x_2, ..., x_d], where
% x: state at the start of the interval
% x_1,x_2,...,x_d: states at the collocation points
% Augmented control is z = [u_1, u_2, ..., u_d], where 
% u_1,u_2,...,u_d: controls at the collocation points
w = {x}; % augmented state
w0 = {x0}; % augmented inital state
wn = {xn}; % augmented final state
lbw = lbx; % lower bound for augmented state
ubw = ubx; % upper bound for augmented state
z = {}; % augmented control
z0 = {}; % augmented inital control
zn = {}; % augmented final control
lbz = []; % lower bound for augmented control
ubz = []; % upper bound for augmented control
% Loop over collocation points
for j = 1:d
    xj = SX.sym(['x_' num2str(j)], nx);
    x0j = SX.sym(['x0_' num2str(j)], nx);
    xnj = SX.sym(['xn_' num2str(j)], nx);
    w = [w; {xj}];
    w0 = [w0; {x0j}];
    wn = [wn; {xnj}];
    lbw = [lbw; lbx];
    ubw = [ubw; ubx];
    uj = SX.sym(['u_' num2str(j)], nu);
    u0j = SX.sym(['u0_' num2str(j)], nu);
    unj = SX.sym(['un_' num2str(j)], nu);
    z = [z; {uj}];
    z0 = [z0; {u0j}];
    zn = [zn; {unj}];
    lbz = [lbz; lbu];
    ubz = [ubz; ubu];
end

% Integrator using collocation points
h = SX.sym('h'); % time step
f = Function('f', {t, x, u, p}, {xdot, l, q}); % ocp dynamics [xdot,l,q] = f(t,x,u,p)
x1 = SX.sym('x1', nx); % initial state of interval 
x2 = SX.sym('x2', nx); % final state of interval 
% Create augmented states for x1, x2
w1 = {x1}; % augmented state at the start of interval
w2 = {x2}; % augmented state at the end of interval
for j = 1:d
    x1j = SX.sym(['x1_' num2str(j)], nx);
    x2j = SX.sym(['x2_' num2str(j)], nx);
    w1 = [w1; {x1j}];
    w2 = [w2; {x2j}];
end
% Loop over collocation points
dw = []; % collocation equations
dl = 0; % cost over one mesh interval
dq = 0; % int constraint over one mesh interval
xend = D(1)*x1; % contribution to the end state due to the start state
for j = 1:d
    xp = C(1,j)*w1{1};
    for r=1:d
        % Expression for the state derivative at the collocation point
        xp = xp + C(r+1,j)*w1{r+1};
    end
    % Append collocation equations
    [fj, lj, qj] = f(t+tau_root(j)*h, w1{j+1}, z{j}, p);
    dw = [dw; h*fj-xp];
    % Add contribution to the end state
    xend = xend + D(j+1)*w1{j+1};
    % Add contribution to quadrature function
    dl = dl + B(j)*lj*h;
    dq = dq + B(j)*qj*h;
end
% Continuity constraint of the end state
dw = [dw; x2-xend];

% CAT augmented state & control to obtain vectors instead of cells
w = vertcat(w{:});
w1 = vertcat(w1{:});
w2 = vertcat(w2{:});
w0 = vertcat(w0{:});
wn = vertcat(wn{:});
z = vertcat(z{:});
z0 = vertcat(z0{:});
zn = vertcat(zn{:});

% OCP functions
% use augmented state w instead of x
ocp_runcost = Function('ocp_runcost', {t, w1, z, w2, p, h}, {dl});
ocp_bcscost = Function('ocp_bcscost', {w0, z0, wn, zn, p}, {m});
ocp_dyn = Function('ocp_dyn', {t, w1, z, w2, p, h}, {dw});
ocp_path = Function('ocp_path', {t, w, z, p}, {c});
ocp_bcs = Function('ocp_bcs', {w0, z0, wn, zn, p}, {b});
ocp_int = Function('ocp_int', {t, w1, z, w2, p, h}, {dq});

% Build OCP
buildOCP(name, ...
        ocp_runcost, ocp_bcscost, ocp_dyn, ocp_path, ocp_bcs, ocp_int, ...
        skip_hessian);

% Create guess
x0 = [linspace(xi(1),xi(1),N)
      linspace(xi(2),xi(2),N)];
u0 = linspace(0,0,N);
w0 = x0;
z0 = [];
for j = 1 : d
    w0 = [w0; x0];
    z0 = [z0; u0];
end

% Create problem structure
problem.name = name;
problem.N = N;
problem.ti = ti;
problem.tf = tf;
problem.guess.x = w0;
problem.guess.u = z0;
problem.guess.p = [];
problem.bounds.lbx = lbw; problem.bounds.ubx = ubw;
problem.bounds.lbu = lbz; problem.bounds.ubu = ubz;
problem.bounds.lbp = lbp; problem.bounds.ubp = ubp;
problem.bounds.lbc = lbc; problem.bounds.ubc = ubc;
problem.bounds.lbb = lbb; problem.bounds.ubb = ubb;
problem.bounds.lbq = lbq; problem.bounds.ubq = ubq;

% Call ocp solver
sol = minosMex(problem);
clear minosMex

% Plot the solution
figure
plot(sol.t, sol.x(1,:),sol.t, sol.x(2,:));
xlabel('t'), ylabel('state');
legend('x1','x2')

figure
plot(sol.t, sol.u(1,:));
xlabel('t'), ylabel('control');
legend('u')

figure
subplot(2, 1, 1)
plot(0:sol.stats.num_iter, sol.stats.obj_history)
xlabel('Iteration')
ylabel('Objective')
subplot(2, 1, 2)
semilogy(0:sol.stats.num_iter, sol.stats.infpr_history, ...
         0:sol.stats.num_iter, sol.stats.infdu_history)
xlabel('Iteration')
ylabel('Max norm')
legend('Feasibility', 'Optimality')
xlabel('Iteration')
ylabel('Max norm')