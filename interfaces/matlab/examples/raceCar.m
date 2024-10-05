%---------------------------- RaceCar Problem ----------------------------%
% Free-trajectory minimum-lap-time of a race car with a basic model       %
% Independent variable is s (travelled distance)                          %
% Problem states x are V, n, chi                                          %
% Problem controls u are ax, ay                                           %
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
name = 'raceCar'; % name of the problem
scheme = 'trapz'; %integration scheme: trapz, rk4, euler
skip_hessian = false; % use approximated hessian
N = 500; % number of mesh points

% Data
g = 9.806; % gravity (m/s2)
axmax = 1; % max ax (g)
aymax = 1.5; % max ay (g)
Vmax = 100; % max V (m/s)
s0 = [0 200 250 350 400 800 850 950 1000 1200]; % track s
kappa0 = [0 0 pi/150 pi/150 0 0 pi/150 pi/150 0 0]; % track curvature
rw = 12; % track width

% Bounds
lbx = [0 -rw/2 -pi/2];
ubx = [Vmax +rw/2 +pi/2];
lbu = [-axmax -aymax];
ubu = [+axmax +aymax];
lbp = [];
ubp = [];
lbc = -1;
ubc = 1;
lbb = zeros(5,1);
ubb = zeros(5,1);
lbq = [];
ubq = [];

% Create CASADI variables (use MX to allow using casadi.interpolant)
nx = 3; nu = 2; np = 0; % state, control, and parameter dimensions
s = MX.sym('s'); % independent variable
x = MX.sym('x', nx); % state
u = MX.sym('u', nu); % control
p = MX.sym('p', np); % parameter
x0 = MX.sym('x0', nx); % initial state
u0 = MX.sym('u0', nu); % initial control
xn = MX.sym('xn', nx); % final state
un = MX.sym('un', nu); % final control

% Get variables
V = x(1);
n = x(2);
chi = x(3);
ax = u(1);
ay = u(2);
kappainterp = casadi.interpolant('kappa','linear',{s0},kappa0);
kappa = kappainterp(s);
sdot = V*cos(chi)/(1-n*kappa);

% Model equations
xdot = [ax*g/sdot
        V*sin(chi)/sdot
        ay*g/V/sdot-kappa]; % state derivative
c = (ax/axmax)^2 + (ay/aymax)^2; % path constraint
q = []; % integral constraint
l = 1/sdot; % running (Lagrange) cost
m = 0; % boundary (Mayer) cost

% Boundary condtions
b = [x0-xn;
     u0-un
     ];

% Integrator
h = MX.sym('h'); % time step
f = Function('f', {s, x, u, p}, {xdot, l, q}); % ocp dynamics [xdot,l,q] = f(t,x,u,p)
switch scheme
    case 'euler' % Euler integration
        x1 = MX.sym('x1', nx); % initial state of interval 
        x2 = MX.sym('x2', nx); % final state of interval 
        [f1, dl1, dq1] = f(s, x1, u, p);
        xnext = x1 + h*f1; % euler 
        dl = h*dl1; % euler 
        dx = xnext - x2;
        dq = h*dq1;
    case 'trapz'
        x1 = MX.sym('x1', nx); % initial state of interval 
        x2 = MX.sym('x2', nx); % final state of interval 
        [f1, dl1, dq1] = f(s,   x1, u, p); 
        [f2, dl2, dq2] = f(s+h, x2, u, p); 
        xnext = x1 + 1/2 * h * (f1+f2);
        dx = xnext - x2;
        dl = 1/2 * h * (dl1+dl2);
        dq = 1/2 * h * (dq1+dq2);
    case 'rk4' % Fixed step Runge-Kutta 4 integrator
        M = 4; % RK4 steps per interval
        x1 = MX.sym('x1', nx); % initial state of interval 
        x2 = MX.sym('x2', nx); % final state of interval 
        xnext = x1; % init
        dl = 0; % init
        dq = 0; % init
        for j=1:M
           [k1, k1_l, k1_q] = f(s, xnext, u, p);
           [k2, k2_l, k2_q] = f(s + h/M/2, xnext + h/M/2 * k1, u, p);
           [k3, k3_l, k3_q] = f(s + h/M/2, xnext + h/M/2 * k2, u, p);
           [k4, k4_l, k4_q] = f(s + h/M, xnext + h/M * k3, u, p);
           xnext = xnext + h/M/6*(k1 +2*k2 +2*k3 +k4);
           dl = dl + h/M/6*(k1_l + 2*k2_l + 2*k3_l + k4_l);
           dq = dq + h/M/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
        end
        dx = xnext - x2;
end

% OCP functions
% OCP runnig cost function dl = ocp_runcost(s, x1, u, x2, p, h)
% s: space at step k
% x1: state at step k
% x2: state at step k+1
% u: control over the interval [k,k+1]
% p: scalar parameter
% h: time step
% dl: integral (i.e. Lagrange) cost over the interval [k,k+1]
ocp_runcost = Function('ocp_runcost', {s, x1, u, x2, p, h}, {dl});

% OCP boundary cost function m = ocp_bcscost(x0, u0, xn, un, p)
% x0: initial state
% u0: initial control
% xn: final state
% un: final control
% p: scalar parameter
% m: boundary (i.e. Mayer) cost
ocp_bcscost = Function('ocp_bcscost', {x0, u0, xn, un, p}, {m});

% OCP dynamic function dx = ocp_dyn(s, x1, u, x2, p, h)
% s: space at step k
% x1: state at step k
% x2: state at step k+1
% u: control over the interval [k,k+1]
% p: scalar parameter
% h: time step
% dx: rhs of integrator over the interval [k,k+1] (which must be dx=0 for
% continuity)
ocp_dyn = Function('ocp_dyn', {s, x1, u, x2, p, h}, {dx});

% OCP path function c = ocp_path(s, x, u, p)
% s: space at step k
% x: state at step k
% u: control over the interval [k,k+1]
% p: scalar parameter
% c: path constraint at step k
ocp_path = Function('ocp_path', {s, x, u, p}, {c});

% OCP boundary condition function b = ocp_bcs(x0, u0, xn, un, p)
% x0: initial state
% u0: initial control
% xn: final state
% un: final control
% p: scalar parameter
% b: boundary condtions
ocp_bcs = Function('ocp_bcs', {x0, u0, xn, un, p}, {b});

% OCP integral constraint function dq = ocp_int(s, x1, u, x2, p, h)
% s: space at step k
% x1: state at step k
% x2: state at step k+1
% u: control over the interval [k,k+1]
% p: scalar parameter
% h: time step
% dq: integral constraint over the interval [k,k+1]
ocp_int = Function('ocp_int', {s, x1, u, x2, p, h}, {dq});

% Build OCP
buildOCP(name, ...
        ocp_runcost, ocp_bcscost, ocp_dyn, ocp_path, ocp_bcs, ocp_int, ...
        skip_hessian);

% Guess
x0 = [Vmax/2*ones(1,N)
      linspace(0,0,N)
      linspace(0,0,N)];
u0 = [linspace(0,0,N)
      linspace(0,0,N)];

% Create problem structure
problem.name = name;
problem.N = N;
problem.ti = s0(1);
problem.tf = s0(end);
problem.guess.x = x0;
problem.guess.u = u0;
problem.guess.p = [];
problem.bounds.lbx = lbx; problem.bounds.ubx = ubx;
problem.bounds.lbu = lbu; problem.bounds.ubu = ubu;
problem.bounds.lbp = lbp; problem.bounds.ubp = ubp;
problem.bounds.lbc = lbc; problem.bounds.ubc = ubc;
problem.bounds.lbb = lbb; problem.bounds.ubb = ubb;
problem.bounds.lbq = lbq; problem.bounds.ubq = ubq;

% Initializa using approximated Hessian with few iterations
problem.options.flag_hessian = true;
problem.options.max_iter = 80;
init_sol = minosMex(problem);

% Finalize the problem using exact hessian
problem = init_sol.next_problem;
problem.options.flag_hessian = false;
problem.options.max_iter = 3000;
sol = minosMex(problem);
clear minosMex

% Clean
% delete(['raceCarMex.' mexext])

% Plot the solution
figure
subplot(311)
plot(sol.t, sol.x(1,:));
ylabel('V (m/s)');
subplot(312)
plot(sol.t, sol.x(2,:));
yline(rw/2); yline(-rw/2);
ylabel('n (m)');
subplot(313)
plot(sol.t, sol.x(3,:));
ylabel('chi (rad)');

figure
subplot(211)
plot(sol.t, sol.u(1,:));
ylabel('ax (g)');
subplot(212)
plot(sol.t, sol.u(2,:));
ylabel('ay (g)');

% Track
kappa = interp1(s0,kappa0,sol.t);
theta = cumtrapz(sol.t,kappa);
x0 = cumtrapz(sol.t,cos(theta));
y0 = cumtrapz(sol.t,sin(theta));
xr = x0 - rw/2*sin(theta);
yr = y0 + rw/2*cos(theta);
xl = x0 + rw/2*sin(theta);
yl = y0 - rw/2*cos(theta);
x = x0 - sol.x(2,:).*sin(theta);
y = y0 + sol.x(2,:).*cos(theta);

figure
plot(x0,y0,'k-.',xl,yl,'k-',xr,yr,'k-',x,y,'r-')
hold off
set(gca,'YDir','reverse');
axis equal

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