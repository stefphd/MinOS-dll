%-------------- Classical Brachistochrone Problem with Auxdata -----------%
% This example is identitical to brachistochrone.m, but                   %
% the OCP includes axmax,aymax as auxdata, in order to allow to run       %
% several simulations with different sets of data without the need to     %
% generate the MEX function again.                                        %
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
name = 'brachistochroneAuxdata'; % name of the problem
scheme = 'trapz'; %integration scheme: trapz, rk4, euler
skip_hessian = false; % use approximated hessian
N = 500; % number of mesh points
ti = 0; % initial time
tf = 1; % final time

% Data
% g = 10; % now this is an auxdata

% Boundary conditions
xi = [0; 0; 0]; % initial x,y,v
xf = [2; 2]; % final x,y

% Bounds
lbx = [0; 0; -50]; % lower bound for state
ubx = [10; 10; 50]; % upper bound for state
lbu = -pi/2; % lower bound for control
ubu = pi/2; % upper bound for control
lbp = 0; % lower bound for parameter
ubp = 2; % upper bound for parameter
lbc = []; % lower bound for path constraint
ubc = []; % upper bound for path constraint
lbb = [xi; xf]; % lower bound for boundary conditions
ubb = [xi; xf]; % upper bound for boundary conditions
lbq = []; % lower bound for integral constraints
ubq = []; % upper bound for integral constraints

% Create CASADI variables
nx = 3; nu = 1; np = 1; na = 1; % state, control, parameter, and auxdata dimensions
t = SX.sym('t'); % independent variable
x = SX.sym('x', nx); % state
u = SX.sym('u', nu); % control
p = SX.sym('p', np); % parameter
x0 = SX.sym('x0', nx); % initial state
u0 = SX.sym('u0', nu); % initial control
xn = SX.sym('xn', nx); % final state
un = SX.sym('un', nu); % final control
auxdata = SX.sym('auxdata', na); % auxdata

% Model equations ( using time scaled by Tf=p(1) )
g = auxdata(1); % g is an auxdata
xdot = p(1)*[x(3)*sin(u(1));
             x(3)*cos(u(1));
             g*cos(u(1))]; % state derivative
c = []; % path constraint
q = []; % integral constraint
l = 0; % running (Lagrange) cost
m = p(1); % boundary (Mayer) cost
% equivalently, we could also use:
% l = p(1);
% m = 0;

% Boundary conditions
b = [x0;
     xn(1:2)];

% Integrator
h = SX.sym('h'); % time step
f = Function('f', {t, x, u, p, auxdata}, {xdot, l, q}); % ocp dynamics [xdot,l,q] = f(t,x,u,p,auxdata)
switch scheme
    case 'euler' % Euler integration
        x1 = SX.sym('x1', nx); % initial state of interval 
        x2 = SX.sym('x2', nx); % final state of interval 
        [f1, dl1, dq1] = f(t, x1, u, p, auxdata);
        xnext = x1 + h*f1; % euler 
        dl = h*dl1; % euler 
        dx = xnext - x2;
        dq = h*dq1;
    case 'trapz'
        x1 = SX.sym('x1', nx); % initial state of interval 
        x2 = SX.sym('x2', nx); % final state of interval 
        [f1, dl1, dq1] = f(t,   x1, u, p, auxdata); 
        [f2, dl2, dq2] = f(t+h, x2, u, p, auxdata); 
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
           [k1, k1_l, k1_q] = f(t, xnext, u, p, auxdata);
           [k2, k2_l, k2_q] = f(t + h/M/2, xnext + h/M/2 * k1, u, p, auxdata);
           [k3, k3_l, k3_q] = f(t + h/M/2, xnext + h/M/2 * k2, u, p, auxdata);
           [k4, k4_l, k4_q] = f(t + h/M, xnext + h/M * k3, u, p, auxdata);
           xnext = xnext + h/M/6*(k1 +2*k2 +2*k3 +k4);
           dl = dl + h/M/6*(k1_l + 2*k2_l + 2*k3_l + k4_l);
           dq = dq + h/M/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
        end
        dx = xnext - x2;
end

% OCP functions
% In this case the OCP functions have also 'auxdata' as the last argument.
% OCP runnig cost function dl = ocp_runcost(t, x1, u, x2, p, h, auxdata)
% t: time at step k
% x1: state at step k
% x2: state at step k+1
% u: control over the interval [k,k+1]
% p: scalar parameter
% h: time step
% auxdata: auxdata
% dl: integral (i.e. Lagrange) cost over the interval [k,k+1]
ocp_runcost = Function('ocp_runcost', {t, x1, u, x2, p, h, auxdata}, {dl});

% OCP boundary cost function m = ocp_bcscost(x0, u0, xn, un, p, auxdata)
% x0: initial state
% u0: initial control
% xn: final state
% un: final control
% p: scalar parameter
% auxdata: auxdata
% m: boundary (i.e. Mayer) cost
ocp_bcscost = Function('ocp_bcscost', {x0, u0, xn, un, p, auxdata}, {m});

% OCP dynamic function dx = ocp_dyn(t, x1, u, x2, p, h, auxdata)
% t: time at step k
% x1: state at step k
% x2: state at step k+1
% u: control over the interval [k,k+1]
% p: scalar parameter
% h: time step
% auxdata: auxdata
% dx: rhs of integrator over the interval [k,k+1] (which must be dx=0 for
% continuity)
ocp_dyn = Function('ocp_dyn', {t, x1, u, x2, p, h, auxdata}, {dx});

% OCP path function c = ocp_path(t, x, u, p, auxdata)
% t: time at step k
% x: state at step k
% u: control over the interval [k,k+1]
% p: scalar parameter
% auxdata: auxdata
% c: path constraint at step k
ocp_path = Function('ocp_path', {t, x, u, p, auxdata}, {c});

% OCP boundary condition function b = ocp_bcs(x0, u0, xn, un, p, auxdata)
% x0: initial state
% u0: initial control
% xn: final state
% un: final control
% p: scalar parameter
% auxdata: auxdata
% b: boundary condtions
ocp_bcs = Function('ocp_bcs', {x0, u0, xn, un, p, auxdata}, {b});

% OCP integral constraint function dq = ocp_int(t, x1, u, x2, p, h)
% t: time at step k
% x1: state at step k
% x2: state at step k+1
% u: control over the interval [k,k+1]
% p: scalar parameter
% h: time step
% auxdata: auxdata
% dq: integral constraint over the interval [k,k+1]
ocp_int = Function('ocp_int', {t, x1, u, x2, p, h, auxdata}, {dq});

% Build OCP
buildOCP(name, ...
        ocp_runcost, ocp_bcscost, ocp_dyn, ocp_path, ocp_bcs, ocp_int, ...
        skip_hessian);

% Create guess
x0 = [linspace(xi(1),xf(1),N)
      linspace(xi(2),xf(2),N)
      linspace(0,0,N)];
u0 = linspace(0,1,N);
p0 = 1;

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
problem.auxdata = 10; % set also auxdata

% Solve with auxdata = 10
sol1 = minosMex(problem);

% Change auxdata, and solve again starting from previous solution
problem = sol1.next_problem;
problem.auxdata = 5;
sol2 = minosMex(problem);

clear minosMex

% Plot the solution
figure
plot(sol1.t, sol1.x, '-', sol2.t, sol2.x, '--');
xlabel('t'), ylabel('state');
legend('x1','y1', 'v1','x2','y2', 'v2')

figure
plot(sol1.t, sol1.u, '-' ,sol2.t, sol2.u, '--');
xlabel('t'), ylabel('control');
legend('u1','u2')

figure
plot(sol1.x(1,:), sol1.x(2,:), '-', sol2.x(1,:), sol2.x(2,:), '--');
xlabel('x'), ylabel('y');
legend('1','2')
axis equal
set(gca,'YDir','reverse');