% BUILDOCP - Generate the C code using CASADI and build the library file with name 'name'.
% 
% Syntax
%
% BUILDOCP(name, ocp_runcost, ocp_bcscost, ocp_dyn, ocp_path, ocp_bcs, ocp_int)
%
% BUILDOCP(name, ocp_runcost, ocp_bcscost, ocp_dyn, ocp_path, ocp_bcs, ocp_int, skip_hessian)
%
% BUILDOCP(name, ocp_runcost, ocp_bcscost, ocp_dyn, ocp_path, ocp_bcs, ocp_int, skip_hessian, outdir)
%
% Input Arguments
%   name     - Name of the problem and the generated library file
%   ocp_runcost - CASADI running cost function with signature dl = ocp_runcost(t, x1, u, x2, p, h, auxdata). 
%              The function computes the integral of the running (i.e. lagrange) cost
%              over one mesh interval, with initial time t, initial state x1, 
%              final state x2, control u (constant during the interval), parameter p,
%              time step h, and (optionally) auxdata.
%   ocp_bcscost - CASADI boundary cost function with signature m = ocp_bcscost(x0, u0, xn, un, p, auxdata). 
%              The function computes the boundary (i.e. Mayer) cost at the
%              initial state x0 and control u0, and final state xn,
%              control un, parameter p, and (optionally) auxdata.
%   ocp_dyn  - CASADI dynamic function with signature f = ocp_dyn(t, x1, u, x2, p, h, auxdata). 
%              The function integrates the OCP differential equations f(x,xdot,u,t)=0 over
%              one mesh interval, with initial time t, initial state x1, final state
%              x2, control u (constant during the interval), parameter p, and time step h.
%   ocp_path - CASADI path constraint with signature c = ocp_path(t, x, u, p, auxdata). 
%              The function calculates the OCP path constraint at one
%              mesh point, at time t, state x, control u, parameter p, and (optionally) auxdata.
%   ocp_bcs  - CASADI boundary conditions with signature b = ocp_bcs(x0,u0,xn,un, auxdata).
%              The function enforces the boundary conditions at the
%              initial state x0 and control u0, and final state xn,
%              control un, parameter p, and (optionally) auxdata.
%   ocp_int  - CASADI integral constraints with signature q = ocp_int(t, x1, u, x2, p, h, auxdata). 
%              The function computes the integral constraint over one mesh interval,
%              with initial time t, initial state x1, final state x2, control  u (constant 
%              during the interval), parameter p, time step h, and (optionally) auxdata.
%
% Optional Input Arguments
%   skip_hessian - FALSE (defualt) or TRUE. Set to TRUE to skip the exact Hessian calculation. 
%              In such case approximated Hessian is used.
%   outdir   - Output directory (default is current working directory, i.e. './'). 
%              Must end with '/'.
%
% Copyright (C) 2024 <a href="mailto:stefano.lovato@unipd.it">Stefano Lovato</a>