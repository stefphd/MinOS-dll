% ----------- Space Station Attitude Control Example -----------%
% This example is taken verbatim from the following reference:  %
% Betts, J. T., Practical Methods for Optimal Control Using     %
% Nonlinear Programming, SIAM Press, Philadelphia, 2009.        %
% --------------------------------------------------------------%

clc, clear

if ispc
    addpath('./../../bin');
else
    addpath('./../../lib');
end

% import casadi
import casadi.*

% Solving options
name = 'spaceStation'; % name of the problem
scheme = 'trapz'; %integration scheme: trapz, rk4, euler
skip_hessian = false; % use approximated hessian
N = 500; % number of mesh points
ti = 0; % initial time
tf = 1800; % final time

% Boundary conditions
omegai = [-9.5380685844896e-6; -1.1363312657036e-3; 5.3472801108427e-6];
ri = [2.9963689649816e-3; 1.5334477761054e-1; 3.8359805613992e-3];
hi = [5000; 5000; 5000];
hf = [0;0;0];

% Data
J = [2.80701911616e7 4.822509936e5 -1.71675094448e7
     4.822509936e5 9.5144639344e7 6.02604448e4
     -1.71675094448e7 6.02604448e4 7.6594401336e7];
omegaorb = 0.06511*pi/180;
angmommax = 10000;

% Bounds
omegamin = [-2e-3; -2e-3; -2e-3];
omegamax = [+2e-3; +2e-3; +2e-3];
rmin = [-1; -1; -1];
rmax = [+1; +1; +1];
hmin = [-15000; -15000; -15000];
hmax = [+15000; +15000; +15000];
umin = -150;
umax = +150;

% Problem scaling
% This problem has varibles which may have large values,
% scaling is thus quite important here
% Here, variables are scaled according to their bounds
tscale = tf;
omegascale = omegamax;
rscale = rmax;
hscale = hmax;
xscale = [omegascale; rscale; hscale];
uscale = umax;
cscale = angmommax^2;

% Scaled time
ti = ti/tscale;
tf = tf/tscale;

% Scaled bounds
lbx = [omegamin; rmin; hmin]./xscale;
ubx = [omegamax; rmax; hmax]./xscale;
lbu = [umin; umin; umin]/uscale;
ubu = [umax; umax; umax]/uscale;
lbp = [];
ubp = [];
lbc = 0/cscale;
ubc = angmommax^2/cscale;
lbb = [omegai./omegascale; ri./rscale; hi./hscale; hf./hscale; zeros(9,1)];
ubb = [omegai./omegascale; ri./rscale; hi./hscale; hf./hscale; zeros(9,1)];
lbq = [];
ubq = [];

% Create CASADI variables
nx = 3*3; nu = 3; np = 0; % state, control, and parameter dimensions
t = SX.sym('t'); % independent variable
x = SX.sym('x', nx); % state
u = SX.sym('u', nu); % control
p = SX.sym('p', np); % parameter
x0 = SX.sym('x0', nx); % initial state
u0 = SX.sym('u0', nu); % initial control
xn = SX.sym('xn', nx); % final state
un = SX.sym('un', nu); % final control

% Get variables
xx = x.*xscale;
uu = u.*uscale;
omega = xx(1:3); 
r = xx(4:6); 
h = xx(7:9);
% Model equations
Sr = skew(r);
C = 2/(1+dot(r,r))*(Sr*(Sr - eye(3))) + eye(3);
taugg = 3*omegaorb^2*cross(C(:,3), J*C(:,3));
omegar = -omegaorb*C(:,2);
omegadot = inv(J)*(taugg - cross(omega, J*omega+h) - uu);
rdot = 1/2*(r*r' + eye(3) + Sr)*(omega - omegar);
hdot = uu;
xdot = [omegadot; rdot; hdot]./xscale*tscale;
c = dot(h,h)/cscale;
q = []; % integral constraint
l = 1e-6*dot(uu,uu)/tscale;
m = 0;
f = Function('f', {t, x, u, p}, {xdot, l}); % ocp dynamics [xdot,l] = f(t,x,u,p)

% Boundary conditions
[xdotf, ~] = f(tf, xn, un, []);
b = [x0;
     xn(7:9);
     xdotf(1:3);
     xdotf(4:6);
     un;
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

buildOCP(name, ...
        ocp_runcost, ocp_bcscost, ocp_dyn, ocp_path, ocp_bcs, ocp_int, ...
        skip_hessian);

% Create guess
x0 = [repmat(omegai,[1 N])./omegascale
      repmat(ri,[1 N])./rscale  
      repmat(hi,[1 N])./hscale];
u0 = zeros(nu, N);
p0 = [];

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

% Clean
% delete(['spaceStationMex.' mexext])

figure
plot(sol.t*tscale,sol.x(1:3,:).*omegascale);
legend('omega1','omega2','omega3');
xlabel('t'), ylabel('state');

figure
plot(sol.t*tscale,sol.x(4:6,:).*rscale);
legend('r1','r2','r3');
xlabel('t'), ylabel('state');

figure
plot(sol.t*tscale,sol.x(7:9,:).*hscale);
legend('h1','h2','h3');
xlabel('t'), ylabel('state');

figure
plot(sol.t(1:end-1)*tscale,sol.u(:,1:end-1).*uscale);
legend('u1','u2','u3');
xlabel('t'), ylabel('control');

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