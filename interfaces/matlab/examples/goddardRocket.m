%---------------------- Launch Vehicle Ascent Example ----------------------%
% This example can be found in the following reference:                     %
% Rao, A. V., Benson, D. A., Darby, C. L., Patterson, M. A., Francolin,     % 
% C., Sanders, I., and Huntington, G. T., Algorithm 902: GPOPS, A MATLAB    %
% Software for Solving Multiple-Phase Optimal Control Problems Using The    %
% Gauss Pseudospectral Method, ACM Transactions on Mathematical Software,   %
% Vol. 37, No. 2, April-June 2010, Article No. 22, pp. 1-39.                %
% --------------------------------------------------------------------------%

clc, clear

addpath('./../../bin');

% import casadi
import casadi.*

% Solving options
name = 'goddardRocket'; % name of the problem
scheme = 'trapz'; %integration scheme: trapz, rk4, euler
skip_hessian = false; % use approximated hessian
N = 500; % number of mesh points
ti = 0; % initial time
tf = 1; % final time

% Boundary conditions
hi = 0; % initial height
vi = 0; % initial velocity
mi = 3; % initial mass
mf = 1; % final mass

% Data
g0 = 32.174;
rho0 = 0.002378;
H = 23800;
csqrd = 3.264*g0*H;
cc = sqrt(csqrd);
Tmax = 2*mi*g0;
dragk = 0.7110*Tmax/csqrd;

% Bounds
hmin = 0; % min height
hmax = 30000; % max height
vmin = -1500; % min velocity
vmax = 1500; % max velocity
mmin = 0.2*mi; %min mass
mmax = mi; % max mass
Tfmin = 0; % min final time
Tfmax = 500; % max final time
lbx = [hmin; vmin; mmin]; % lower bound for state
ubx = [hmax; vmax; mmax]; % upper bound for state
lbu = 0; % lower bound for control
ubu = Tmax; % upper bound for control
lbp = Tfmin; % lower bound for parameter
ubp = Tfmax; % upper bound for parameter
lbc = []; % lower bound for path constraint
ubc = []; % upper bound for path constraint
lbb = [hi vi mi 0]; % lower bound for boundary conditions
ubb = [hi vi mi 0]; % upper bound for boundary conditions
lbq = -(mi-mf); % lower bound for integral constraints
ubq = -(mi-mf); % upper bound for integral constraints

% Create CASADI variables
nx = 3; nu = 1; np = 1; % state, control, and parameter dimensions
t = SX.sym('t'); % independent variable
x = SX.sym('x', nx); % state
u = SX.sym('u', nu); % control
p = SX.sym('p', np); % parameter
x0 = SX.sym('x0', nx); % initial state
u0 = SX.sym('u0', nu); % initial control
xn = SX.sym('xn', nx); % final state
un = SX.sym('un', nu); % final control

% Model equations ( using time scaled by Tf=p(1) )
h = x(1); v = x(2); m = x(3); T = u(1); Tf = p(1);
D = dragk.*(v.^2).*exp(-h/H);
xdot = Tf*[v;
           (T-D)/m-g0;
           -T/cc]; % state derivative
c = []; % path constraint
q = -Tf*T/cc;
l = 0; % running (Lagrange) cost
m = -xn(1); % boundary (Mayer) cost: maximize final altitude
% equivalently, we could also use:
% l = -xdot(1);
% m = 0;

% Boundary conditions
b = [x0(1); % initial altitude
     x0(2); % initial velocity
     x0(3); % initial mass
     un(1); % just to avoid undetermined final control
     ];

% Integrator
h = SX.sym('h'); % time step
f = Function('f', {t, x, u, p}, {xdot, l, q}); % ocp dynamics [xdot,l,q] = f(t,x,u,p)
switch scheme
    case 'euler' % Euler integration
        x1 = SX.sym('x1', nx); % initial state of interval 
        x2 = SX.sym('x2', nx); % final state of interval 
        [f1, dl1, dq1] = f(t, x1, u, p);
        xnext = x1 + h*f1; % euler 
        dl = h*dl1; % euler 
        dx = xnext - x2;
        dq = h*dq1;
    case 'trapz'
        x1 = SX.sym('x1', nx); % initial state of interval 
        x2 = SX.sym('x2', nx); % final state of interval 
        [f1, dl1, dq1] = f(t,   x1, u, p); 
        [f2, dl2, dq2] = f(t+h, x2, u, p); 
        xnext = x1 + 1/2 * h * (f1+f2);
        dx = xnext - x2;
        dl = 1/2 * h * (dl1+dl2);
        dq = 1/2 * h * (dq1+dq2);
    case 'rk4' % Fixed step Runge-Kutta 4 integrator
        M = 4; % RK4 steps per interval
        x1 = SX.sym('x1', nx); % initial state of interval 
        x2 = SX.sym('x2', nx); % final state of interval 
        xnext = x1; % init
        dl = 0; % init
        dq = 0; % init
        for j=1:M
           [k1, k1_l, k1_q] = f(t, xnext, u, p);
           [k2, k2_l, k2_q] = f(t + h/M/2, xnext + h/M/2 * k1, u, p);
           [k3, k3_l, k3_q] = f(t + h/M/2, xnext + h/M/2 * k2, u, p);
           [k4, k4_l, k4_q] = f(t + h/M, xnext + h/M * k3, u, p);
           xnext = xnext + h/M/6*(k1 +2*k2 +2*k3 +k4);
           dl = dl + h/M/6*(k1_l + 2*k2_l + 2*k3_l + k4_l);
           dq = dq + h/M/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
        end
        dx = xnext - x2;
end

% OCP functions
% OCP runnig cost function dl = ocp_runcost(t, x1, u, x2, p, h)
% t: time at step k
% x1: state at step k
% x2: state at step k+1
% u: control over the interval [k,k+1]
% p: scalar parameter
% h: time step
% dl: integral (i.e. Lagrange) cost over the interval [k,k+1]
ocp_runcost = Function('ocp_runcost', {t, x1, u, x2, p, h}, {dl});

% OCP boundary cost function m = ocp_bcscost(x0, u0, xn, un, p)
% x0: initial state
% u0: initial control
% xn: final state
% un: final control
% p: scalar parameter
% m: boundary (i.e. Mayer) cost
ocp_bcscost = Function('ocp_bcscost', {x0, u0, xn, un, p}, {m});

% OCP dynamic function dx = ocp_dyn(t, x1, u, x2, p, h)
% t: time at step k
% x1: state at step k
% x2: state at step k+1
% u: control over the interval [k,k+1]
% p: scalar parameter
% h: time step
% dx: rhs of integrator over the interval [k,k+1] (which must be dx=0 for
% continuity)
ocp_dyn = Function('ocp_dyn', {t, x1, u, x2, p, h}, {dx});

% OCP path function c = ocp_path(t, x, u, p)
% t: time at step k
% x: state at step k
% u: control over the interval [k,k+1]
% p: scalar parameter
% c: path constraint at step k
ocp_path = Function('ocp_path', {t, x, u, p}, {c});

% OCP boundary condition function b = ocp_bcs(x0, u0, xn, un, p)
% x0: initial state
% u0: initial control
% xn: final state
% un: final control
% p: scalar parameter
% b: boundary condtions
ocp_bcs = Function('ocp_bcs', {x0, u0, xn, un, p}, {b});

% OCP integral constraint function dq = ocp_int(t, x1, u, x2, p, h)
% t: time at step k
% x1: state at step k
% x2: state at step k+1
% u: control over the interval [k,k+1]
% p: scalar parameter
% h: time step
% dq: integral constraint over the interval [k,k+1]
ocp_int = Function('ocp_int', {t, x1, u, x2, p, h}, {dq});

% Build OCP
buildOCP(name, ...
        ocp_runcost, ocp_bcscost, ocp_dyn, ocp_path, ocp_bcs, ocp_int, ...
        skip_hessian);

% Create guess
x0 = [linspace(hi,hmax,N)
      linspace(vi,vi,N)
      linspace(mi,mf,N)];
u0 = linspace(0,Tmax,N);
p0 = Tfmax;

% Create problem structure
problem.name = name;
problem.N = N;
problem.ti = ti;
problem.tf = tf;
problem.guess.x = x0;
problem.guess.u = u0;
problem.guess.p = p0;
problem.bounds.lbx = lbx; problem.bounds.ubx = ubx;
problem.bounds.lbp = lbp; problem.bounds.ubp = ubp;
problem.bounds.lbu = lbu; problem.bounds.ubu = ubu;
problem.bounds.lbc = lbc; problem.bounds.ubc = ubc;
problem.bounds.lbb = lbb; problem.bounds.ubb = ubb;
problem.bounds.lbq = lbq; problem.bounds.ubq = ubq;

% Call ocp solver
sol = minosMex(problem);
clear minosMex

% Plot the solution
figure
plot(sol.t, sol.x);
xlabel('t'), ylabel('state');
legend('h','v', 'm')

figure
plot(sol.t, sol.u);
xlabel('t'), ylabel('control');
legend('T')

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