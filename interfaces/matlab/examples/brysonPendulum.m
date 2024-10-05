%------------------------ Bryson Pendulum problem ------------------------%
% This problem is taken from Section 9.3.13 of the following reference:   %
% Bryson, A.E. Dynamic Optimization; Addison-Wesley: Menlo Park, CA, USA, %
% 1999.                                                                   %
%-------------------------------------------------------------------------%

clc, clear

if ispc
    addpath('./../../bin');
else
    addpath('./../../lib');
end

% import casadi
import casadi.*

% Solving options
name = 'brysonPendlumum'; % name of the problem
scheme = 'trapz'; %integration scheme: trapz, rk4, euler
skip_hessian = false; % use approximated hessian
N = 500; % number of mesh points
ti = 0; % initial time
tf = 1; % final time

% Data
epsilon = 0.5; % m/(M+m)
wu = 0; % u penalty
umax = 1; % max u
Tfmax = 10; % max final time

% Boundary conditions
xi = [0 0 0 0]; % initial y, theta, yd, thetad
xf = [0 pi 0 0]; % final y, theta, yd, thetad

% Bounds
lbx = [-5 -2*pi -100 -100]; % lower bound for state
ubx = [+5 +2*pi +100 +100]; % upper bound for state
lbu = -umax; % lower bound for control
ubu = +umax; % upper bound for control
lbp = 0; % lower bound for parameter
ubp = Tfmax; % upper bound for parameter
lbc = []; % lower bound for path constraint
ubc = []; % upper bound for path constraint
lbb = [xi xf 0]; % lower bound for boundary conditions
ubb = [xi xf 0]; % upper bound for boundary conditions
lbq = []; % lower bound for integral constraints
ubq = []; % upper bound for integral constraints

% Create CASADI variables
nx = 4; nu = 1; np = 1; % state, control, and parameter dimensions
t = SX.sym('t'); % independent variable
x = SX.sym('x', nx); % state
u = SX.sym('u', nu); % control
p = SX.sym('p', np); % parameter
x0 = SX.sym('x0', nx); % initial state
u0 = SX.sym('u0', nu); % initial control
xn = SX.sym('xn', nx); % final state
un = SX.sym('un', nu); % final control

% Model equations ( using time scaled by Tf=p(1) )
theta = x(2);
ydot = x(3);
thetadot = x(4);
Tf = p(1);
xdot = Tf * [ydot;
             thetadot;
             (epsilon*sin(theta)*thetadot^2+epsilon*cos(theta)*sin(theta)+u)/(1-epsilon*cos(theta)^2);
             -(epsilon*cos(theta)*sin(theta)*thetadot^2+sin(theta)+u*cos(theta))/(1-epsilon*cos(theta)^2);
            ]; % state derivative
c = []; % path constraint
q = []; % integral constraint
l = Tf*wu*u^2; % running (Lagrange) cost
m = Tf; % boundary (Mayer) cost
% equivalently, we could also use:
% l = p(1) + wu*u^2;
% m = 0;

% Boundary conditions
b = [x0;
     xn;
     un];

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
x0 = [linspace(xi(1),xf(1),N)
      linspace(xi(2),xf(2),N)
      linspace(xi(3),xf(3),N)
      linspace(xi(4),xf(4),N)];
u0 = linspace(0,0,N);
p0 = 5;

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

% Call OCP solver
sol = minosMex(problem);
clear minosMex

% Plot the solution
figure
plot(sol.t, [sol.x(1:2,:); sol.u]);
xlabel('t'), ylabel('state and control');
legend('y','\theta','u')

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