%------------- Brachistochrone Dual Problem - Indirect Method ------------%
% The problem solved here is given as follows:                            %
%   Maximize x_f                                                          %
% subject to the dynamic constraints                                      %
%    dx/dt = v*sin(u)                                                     %
%    dv/dt = g*cos(u)                                                     %
% and the boundary conditions                                             %
%    x(0) = 0, v(0) = 0                                                   %
%    x(t_f) = FREE, v(t_f) = FREE                                         %
% where t_f is FIXED.                                                     %   
% The problem is solved using an indirect approach, with trapezoidal rule %
% used to integrate the dynamics. Thus, the solution is found solving a   %
% nonlinear system of equations. This is solved using as a NLP, where the %
% the constraints are the system of equations, while the cost is 0.       %
%-------------------------------------------------------------------------%

clc, clear

if ispc
     addpath('./../../bin');
 else
     addpath('./../../lib');
 end

% import casadi
import casadi.*

% Data
name = 'brachistochroneDualIndirect'; % name of the problem
skip_hessian = false; % use approximated hessian
g = 10; % gravity
xi = [0; 0]; % initial x,v
ti = 0; % initial time
tf = 1; % final time
N = 100; % number of mesh points

% Variables
nx = 2; nu = 1;
t = SX.sym('t'); % independent variable
x = SX.sym('x', nx); % state
lam = SX.sym('lam',nx); % costate
u = SX.sym('u', nu); % control
p = SX.sym('p', 0); % parameter
x0 = SX.sym('x0', nx); % initial state
lam0 = SX.sym('lam0', nx); % initial costate
u0 = SX.sym('u0', nu); % initial control
xn = SX.sym('xn', nx); % final state 
lamn = SX.sym('lamn', nx); % final costate 
un = SX.sym('un', nu); % final control

% Model equations
f = [x(2)*sin(u(1));
     g*cos(u(1))]; % dynamics
q = []; % integral constraint
l = 0; % lagrange
m = -xn(1); % mayer
H = l + lam'*f; % hamiltonian

% ODE
xdot = gradient(H, lam); % dynamics xdot=DH/dlam
lamdot = -gradient(H, x); % coequations lamdot=DH/dx

% Optimality
Hu = jacobian(H, u); % dH/du
% sol of Hu=0 is
% uopt = atan(lam(1)*x(2) / (g*lam(2)));
c = Hu; % path constraint
% l = H;

% Boundary conditions
b = [x0-xi; % fixed x at initial
     lamn-gradient(m,xn) % free x at final
     ];

% Integrator
h = SX.sym('h'); % time step
f = Function('f', {t, x, lam, u, p}, {xdot, lamdot, l, q}); % ocp dynamics [xdot,lamdot,l,q] = f(t,x,lam,u,p)
x1 = SX.sym('x1', nx); % initial state of interval 
x2 = SX.sym('x2', nx); % final state of interval
lam1 = SX.sym('lam1', nx); % initial costate of interval 
lam2 = SX.sym('lam2', nx); % final costate of interval 
[xd1, ld1, dl1, dq1] = f(t,   x1, lam1, u, p); 
[xd2, ld2, dl2, dq2] = f(t+h, x2, lam2, u, p); 
xnext = x1   + 1/2 * h * (xd1+xd2);
lnext = lam1 + 1/2 * h * (ld1+ld2);
dx = xnext - x2;
dlam = lnext - lam2;
dl = 1/2 * h * (dl1+dl2);
dq = 1/2 * h * (dq1+dq2);

% Collect x,lam into augmented state w
w = [x; lam];
w0 = [x0; lam0];
wn = [xn; lamn];
w1 = [x1; lam1];
w2 = [x2; lam2];
dw = [dx; dlam];

% OCP functions
% use augmented state w instead of x
ocp_runcost = Function('ocp_runcost', {t, w1, u, w2, p, h}, {dl});
ocp_bcscost = Function('ocp_bcscost', {w0, u0, wn, un, p}, {m});
ocp_dyn = Function('ocp_dyn', {t, w1, u, w2, p, h}, {dw});
ocp_path = Function('ocp_path', {t, w, u, p}, {c});
ocp_bcs = Function('ocp_bcs', {w0, u0, wn, un, p}, {b});
ocp_int = Function('ocp_int', {t, w1, u, w2, p, h}, {dq});

% Build OCP
buildOCP(name, ...
        ocp_runcost, ocp_bcscost, ocp_dyn, ocp_path, ocp_bcs, ocp_int, ...
        skip_hessian);

% Create guess
x0 = [linspace(xi(1),1,N)
      linspace(0,0,N)
      linspace(-1, -1, N);
      linspace(-1, 0, N)];
u0 = linspace(0,1,N);
p0 = [];

% Create problem structure
problem.name = name;
problem.N = N;
problem.ti = ti;
problem.tf = tf;
problem.guess.x = x0;
problem.guess.u = u0;
problem.guess.p = p0;
problem.bounds.lbx = -100*ones(size(w)); problem.bounds.ubx = 100*ones(size(w));
problem.bounds.lbp = []; problem.bounds.ubp = [];
problem.bounds.lbu = -pi/2; problem.bounds.ubu = +pi/2;
problem.bounds.lbc = zeros(size(c)); problem.bounds.ubc = zeros(size(c));
problem.bounds.lbb = zeros(size(b)); problem.bounds.ubb = zeros(size(b));
problem.bounds.lbq = []; problem.bounds.ubq = [];

% Call OCP solver
sol = minosMex(problem);
clear minosMex

% Compute y using trapz integration
ydot1 = sol.x(2,1:end-1).*cos( sol.u(1:end-1) );
ydot2 = sol.x(2,2:end).*cos( sol.u(1:end-1) );
ydot = 0.5*(ydot1+ydot2);
h = (tf-ti)/(N-1);
y = [0 cumsum(ydot*h)];

% Plot the solution
figure
plot(sol.t, sol.x);
xlabel('t'), ylabel('state & costate');
legend('x','v','\lambda_x','\lambda_v')

figure
plot(sol.t, sol.u);
xlabel('t'), ylabel('control');
legend('u')

figure
plot(sol.x(1,:), y);
xlabel('x'), ylabel('y');
axis equal
set(gca,'YDir','reverse');

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