-- LUA script file for minos-cli
-- Usage
--   minos-cli -r script.lua

-- Function to create a linear space
function lspace(start, stop, num)
    local step = (stop - start) / (num - 1)
    local arr = {}
    for i = 0, num - 1 do
        arr[i + 1] = start + i * step
    end
    return arr
end

-- Function to create a vector filled with a specific value
function vector(value, length)
    local vec = {}
    for i = 1, length do
        vec[i] = value
    end
    return vec
end

-- Build problem
cfile = 'brachistochrone.c'
outname = 'myproblem'
build(cfile, outname)

-- problem definition
problem = {
    name = outname,                -- name of the problem
    N = 500,                       -- number of mesh points
    ti = 0,                        -- initial time
    tf = 1,                        -- final time
}
-- problem guess
problem.guess = {
	x = {
            lspace(0, 2, problem.N),
            lspace(0, 2, problem.N),
            lspace(0, 0, problem.N),
        },                         -- state guess 
	u = {
            lspace(0, 1, problem.N)
        },                          -- control guess
	p = { 1, },                     -- parameter guess
    -- guess on multipliers - OPTIONAL
    -- lam_x = { },                   -- state multiplier guess
    -- lam_u = { },                   -- control multiplier guess
    -- lam_p = { },                   -- parameter multiplier guess
    -- lam_f = { },                   -- dynamic multiplier guess
    -- lam_c = { },                   -- path multiplier guess
    -- lam_b = { },                   -- boundary multiplier guess
    -- lam_q = { },                   -- integral multiplier guess
}
-- problem bounds
problem.bounds = {
	lbx = {0, 0, -50},         -- lower bound for state
	ubx = {10, 10, 50},        -- upper bound for state
	lbu = { 0 },               -- lower bound for parameter
	ubu = { 2 },               -- upper bound for parameter
	lbp = { -math.pi / 2 },    -- lower bound for control
	ubp = {  math.pi / 2 },    -- upper bound for control
	lbc = { },                 -- lower bound for path constraint
	ubc = { },                 -- upper bound for path constraint
	lbb = {0, 0, 0, 2, 2},     -- lower bound for boundary conditions
	ubb = {0, 0, 0, 2, 2},     -- upper bound for boundary conditions
	lbq = { },                 -- lower bound for integral constraints
	ubq = { },                 -- upper bound for integral constraints
}
-- problem auxdata (OPTIONAL)
-- problem.auxdata = { }          -- auxiliary data
-- problem mesh (OPTIONAL)
-- problem.mesh = vector(1 / (problem.N - 1), problem.N - 1) -- mesh fractions
-- problem options (OPTIONAL)
-- problem.options = {
	-- nlpsolver = 'ipopt',      -- NLP solver to use
	-- flag_hessian = false,     -- use approx Hessian
	-- max_iter = 3000,          -- max number of iterations
	-- mu_init = 1,              -- initial barrier parameter for interior-point NLP solvers (e.g. IPOPT)
	-- logfile = 'nlp.log',      -- NLP log filename
    -- outfile = 'output.txt',   -- Name of the output file (if any)
	-- print_itersol = 0,        -- iteration interval to print output file
    -- display = false,          -- display to monitor
-- }

-- Solve problem
solution = solve(problem)