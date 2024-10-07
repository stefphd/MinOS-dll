function buildOCP(name, ocp_runcost, ocp_bcscost, ocp_dyn, ocp_path, ocp_bcs, ocp_int, skip_hessian, outdir)
    % Checks args
    arguments 
        name (1,:) char
        ocp_runcost casadi.Function {check_fun(ocp_runcost, 'runcost', ocp_runcost)}
        ocp_bcscost casadi.Function {check_fun(ocp_bcscost, 'bcscost', ocp_runcost)}
        ocp_dyn casadi.Function {check_fun(ocp_dyn, 'dyn', ocp_runcost)}
        ocp_path casadi.Function {check_fun(ocp_path, 'path', ocp_runcost)}
        ocp_bcs casadi.Function {check_fun(ocp_bcs, 'bcs', ocp_runcost)}
        ocp_int casadi.Function {check_fun(ocp_int, 'int', ocp_runcost)}
        skip_hessian (1,1) {islogical} = false
        outdir (1,:) char = './'
    end

    % Codegen directories
    if ~endsWith(outdir, '/')
        outdir = strcat(outdir, '/');
    end
    basedir = fileparts(mfilename('fullpath')); % folder of this file
    basedir = [cd(cd(basedir)) '/']; % resolve basedir path

    % Generate ocp functions
    import casadi.* % import casadi
    nx = ocp_runcost.size1_in(1);
    nu = ocp_runcost.size1_in(2);
    np = ocp_runcost.size1_in(4);
    na = -1;
    if (ocp_runcost.n_in() > 6)
        na = ocp_runcost.size1_in(6);
    end
    % variables
    x = casadi.MX.sym('x', nx);
    x1 = casadi.MX.sym('x1', nx); 
    x2 = casadi.MX.sym('x2', nx);
    u = casadi.MX.sym('u', nu); 
    xi = casadi.MX.sym('xi', nx); 
    ui = casadi.MX.sym('ui', nu); 
    xf = casadi.MX.sym('xf', nx);
    uf = casadi.MX.sym('uf', nu);
    p = casadi.MX.sym('p', np);
    t = casadi.MX.sym('t'); 
    h = casadi.MX.sym('h');
    % arguments
    args_runcost = {t, x1, u, x2, p, h};
    args_bcscost = {xi, ui, xf, uf, p};
    args_dyn = {t, x1, u, x2, p, h};
    args_path = {t, x, u, p};
    args_bcs = {xi, ui, xf, uf, p};
    args_int = {t, x1, u, x2, p, h};
    % add auxdata
    if na>=0
        auxdata = casadi.MX.sym('auxdata', na);
        args_runcost{end+1} = auxdata;
        args_bcscost{end+1} = auxdata;
        args_dyn{end+1} = auxdata;
        args_path{end+1} = auxdata;
        args_bcs{end+1} = auxdata;
        args_int{end+1} = auxdata;
    end

    % gradients
    ocp_runcost_grad = casadi.Function('ocp_runcost_grad', args_runcost, ...
                                    {gradient(ocp_runcost(args_runcost{:}), [x1; u; x2; p])});
    ocp_bcscost_grad = casadi.Function('ocp_bcscost_grad', args_bcscost, ...
                                    {gradient(ocp_bcscost(args_bcscost{:}), [xi; ui; xf; uf; p])});
    % jacobians
    ocp_dyn_jac = casadi.Function('ocp_dyn_jac', args_dyn, ...
                                    {jacobian(ocp_dyn(args_dyn{:}), [x1; u; x2; p])});
    ocp_path_jac = casadi.Function('ocp_path_jac', args_path, ...
                                    {jacobian(ocp_path(args_path{:}), [x; u; p])});
    ocp_bcs_jac = casadi.Function('ocp_bcs_jac', args_bcs, ...
                                    {jacobian(ocp_bcs(args_bcs{:}), [xi; ui; xf; uf; p])});
    ocp_int_jac = casadi.Function('ocp_int_jac', args_int, ...
                                    {jacobian(ocp_int(args_int{:}), [x1; u; x2; p])});

    % collect all function
    ocp_funcs = {ocp_dyn, ocp_path, ocp_bcs, ocp_int, ocp_runcost, ocp_bcscost, ...
                 ocp_dyn_jac, ocp_path_jac, ocp_bcs_jac, ocp_int_jac, ocp_runcost_grad, ocp_bcscost_grad}; 

    % hessians if not skip_hessian
    if ~skip_hessian
        % num of path and bcs
        nc = ocp_path.size1_out(0);
        nb = ocp_bcs.size1_out(0);
        nq = ocp_int.size1_out(0);
        % create CASADI variables for multipliers
        sigma = casadi.MX.sym('sigma'); % cost multiplier
        lamf = casadi.MX.sym('lamf', nx); % dynamic multiplier
        lamc = casadi.MX.sym('lamc', nc); % path multiplier
        lamb = casadi.MX.sym('lamb', nb); % bcs multiplier
        lamq = casadi.MX.sym('lamq', nq); % int multiplier
        % build lagragians
        lagb = sigma*ocp_bcscost(args_bcscost{:}); % boundary lagragian
        lagi = sigma*ocp_runcost(args_runcost{:}); % internal lagragian
        if nb > 0 % add bcs lagragian
            lagb = lagb + lamb'*ocp_bcs(args_bcs{:});
        end
        if nx > 0 % add dynamic lagragian
            lagi = lagi + lamf'*ocp_dyn(args_dyn{:});
        end
        if nc > 0 % add path lagragian
            lagb = lagb + lamc'*ocp_path(args_path{1}, xf, uf, args_path{4:end});
            lagi = lagi + lamc'*ocp_path(args_path{1}, x1, u, args_path{4:end});
        end
        if nq > 0 % add integral lagragian
            lagi = lagi + lamq'*ocp_int(args_int{:});
        end
        % generate hessians of lagragians
        hessb = tril(hessian(lagb,[xi; ui; xf; uf; p]));
        hessi = tril(hessian(lagi,[x1; u; x2; p]));
        % create CASADI functions
        args_hessb = {t, xi, ui, xf, uf, p, sigma, lamc, lamb};
        args_hessi = {t, x1, u, x2, p, h, sigma, lamc, lamf, lamq};
        if na>=0
            args_hessb{end+1} = auxdata;
            args_hessi{end+1} = auxdata;
        end
        ocp_hessb = casadi.Function('ocp_hessb', args_hessb, { hessb });
        ocp_hessi = casadi.Function('ocp_hessi', args_hessi, { hessi });
        % append to ocp_funcs
        ocp_funcs{end+1} = ocp_hessb;
        ocp_funcs{end+1} = ocp_hessi;
    end

    % Generate code
    fprintf("Generating C code...\n")
    cfilename = [name '.c'];
    codegen_options.casadi_int = 'int'; % use casadi_int = int for consistency
    cg = CodeGenerator(cfilename, codegen_options); 
    % Append ocp functions to cg
    for k = 1 : numel(ocp_funcs)
        cg.add(ocp_funcs{k});
    end
    cg.generate(outdir); % generate c and h file in out dir
    pause(0) % just to print out all fprintf

    % For windows add <basedir> to PATH
    if ispc
        add2path(basedir)
    end

    % Compile library
    % Select the C compiler
    if ispc % use local GCC in Windows
        cc = ['"' basedir 'gcc/bin/gcc' '"'];
    else % GCC otherwise
        cc = 'gcc';
    end
    % Test the C compiler
    [exit, ~] = system([cc ' -v']); % ~ is to suppress output messages
    if exit ~= 0
        error('minosMex:buildFailed','Unable to find C compiler ''%s''.', cc);
    end
    % Define library output name
    if ispc % dll for windows
        libext = 'dll'; 
    else % so otherwise
        libext = 'so';
    end
    libname = [name '.' libext];
    % Build command
    cc_args = ['-shared -O1 -fPIC ' outdir cfilename ' -o ' outdir libname];
    cc_cmd = [cc ' ' cc_args];
    % Run the build process
    tic;
    exit = system(cc_cmd);
    compileTime = toc;
    if (exit ~= 0)
        error('minosMex:buildFailed','Unable to build library ''%s'' from file ''%s''', libname, cfilename);
    end
    % Print end message
    fprintf("Library %s built in %.2fs\n", libname, compileTime);
    % Clean
    delete([outdir cfilename])
    % Done
    pause(0) % just to print out all fprintf

end

function add2path(path)
    PATH = getenv('PATH');
    % Add path to the PATH
    if ~contains(PATH, path)
        setenv('PATH',[path ';' PATH]);
    end
end

function check_name(fun, name)
    if ~strcmp(fun.name, name)
        error('buildMex:invalidInputs', '%s must have name ''%s''.', name, name);
    end
end

function check_nargs(fun, nin, nout)
    if (fun.n_in ~= nin) || (fun.n_out ~= nout)
        error('buildMex:invalidInputs','%s must have %d inputs and %d output (found %d and %d).', fun.name, nin, nout, fun.n_in, fun.n_out);
    end
end

function check_argin(fun, i, sz)
    sz0 = fun.size_in(i);
    if ~isequal(sz0, sz)
        error('buildMex:invalidInputs','Input[%d] of %s must have size %d-by-%d (found %d-by-%d).', i, fun.name, sz(1), sz(2), sz0(1), sz0(2));
    end
end

function check_argout(fun, i, sz)
    sz0 = fun.size_out(i);
    if ~isequal(sz0(~isinf(sz)), sz(~isinf(sz)))
        error('buildMex:invalidInputs','Input[%d] of %s must have size %d-by-%d (found %d-by-%d).', i, fun.name, sz(1), sz(2), sz0(1), sz0(2));
    end
end

function check_fun(fun, flag, reffun)
    nx = reffun.size1_in(1);
    nu = reffun.size1_in(2);
    np = reffun.size1_in(4);
    na = -1;
    if (reffun.n_in() > 6)
        na = reffun.size1_in(6);
    end
    switch (flag)
        case 'runcost'
            check_name(fun, 'ocp_runcost');
            check_nargs(fun, 6+(na>=0), 1);
            check_argin(fun, 0, [1 1]);
            check_argin(fun, 1, [nx 1]);
            check_argin(fun, 2, [nu 1]);
            check_argin(fun, 3, [nx 1]);
            check_argin(fun, 4, [np 1]);
            check_argin(fun, 5, [1 1]);
            if na>=0
                check_argin(fun, 6, [na 1]);
            end
            check_argout(fun, 0, [1 1]);
        case 'bcscost'
            check_name(fun, 'ocp_bcscost');
            check_nargs(fun, 5+(na>=0), 1);
            check_argin(fun, 0, [nx 1]);
            check_argin(fun, 1, [nu 1]);
            check_argin(fun, 2, [nx 1]);
            check_argin(fun, 3, [nu 1]);
            check_argin(fun, 4, [np 1]);
            if na>=0
                check_argin(fun, 5, [na 1]);
            end
            check_argout(fun, 0, [1 1]);
        case 'dyn' 
            check_name(fun, 'ocp_dyn');
            check_nargs(fun, 6+(na>=0), 1);
            check_argin(fun, 0, [1 1]);
            check_argin(fun, 1, [nx 1]);
            check_argin(fun, 2, [nu 1]);
            check_argin(fun, 3, [nx 1]);
            check_argin(fun, 4, [np 1]);
            check_argin(fun, 5, [1 1]);
            if na>=0
                check_argin(fun, 6, [na 1]);
            end
            check_argout(fun, 0, [nx 1]);
        case 'path'
            check_name(fun, 'ocp_path');
            check_nargs(fun, 4+(na>=0), 1);
            check_argin(fun, 0, [1 1]);
            check_argin(fun, 1, [nx 1]);
            check_argin(fun, 2, [nu 1]);
            check_argin(fun, 3, [np 1]);
            if na>=0
                check_argin(fun, 4, [na 1]);
            end
            sz = fun.size_out(0);
            if sz(1)>0, check_argout(fun, 0, [inf 1]); end
        case 'bcs'
            check_name(fun, 'ocp_bcs');
            check_nargs(fun, 5+(na>=0), 1);
            check_argin(fun, 0, [nx 1]);
            check_argin(fun, 1, [nu 1]);
            check_argin(fun, 2, [nx 1]);
            check_argin(fun, 3, [nu 1]);
            check_argin(fun, 4, [np 1]);
            if na>=0
                check_argin(fun, 5, [na 1]);
            end
            sz = fun.size_out(0);
            if sz(1)>0, check_argout(fun, 0, [inf 1]); end
        case 'int' 
            check_name(fun, 'ocp_int');
            check_nargs(fun, 6+(na>=0), 1);
            check_argin(fun, 0, [1 1]);
            check_argin(fun, 1, [nx 1]);
            check_argin(fun, 2, [nu 1]);
            check_argin(fun, 3, [nx 1]);
            check_argin(fun, 4, [np 1]);
            check_argin(fun, 5, [1 1]);
            if na>=0
                check_argin(fun, 6, [na 1]);
            end
            sz = fun.size_out(0);
            if sz(1)>0, check_argout(fun, 0, [inf 1]); end
    end
end