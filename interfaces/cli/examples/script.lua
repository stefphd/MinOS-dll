-- LUA script file for minos-cli
-- Usage
--   minos-cli -r script.lua

-- Build problem
cfile = 'brachistochrone.c'
outname = "luaproblem"
build(cfile, outname)

-- Load problem table from problem.lua
require("problem")
problem["name"] = outname -- change problem name to outname

-- Solve problem
solution = solve(problem)

-- Solve again from previous solution
next_problem = solution.next_problem
next_solution = solve(next_problem)