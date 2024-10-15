/* File minos_cli.cpp
    Command line interface for MINOS
    Copyright (C) 2024 Stefano Lovato
*/

#include "minos.h"
#include "minos_header.h"

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <chrono>

#ifdef _WIN32
#include <windows.h>
#define popen _popen
#define pclose _pclose
#else
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif

/* LUA headers */
extern "C" {
    #include "lua.h"
    #include "lualib.h"
    #include "lauxlib.h"
}

#define UNVALID_USAGE -1
#define BUILD_USAGE 0
#define SOLVE_USAGE 1
#define RUN_USAGE 2

const std::string helpstr = 
R"(Command-line interface usage

  minos-cli -b <csource>
  minos-cli -b <csource> -o <outname>
  minos-cli -b <csource> -o <outname> -d <outdir>

Build the OCP library from the C file.

  minos-cli -s <luafile>

Solve the LUA problem.

  minos-cli -r <luafile>

Run the LUA script.

Options:
  -h [--help]            = Print this help
  -v [--version]         = Print the header with version
  -b [--build] <csource> = Specify the name of the C file to buildy
  -o [--out]   <outname> = Specify the name of the output library
  -d [--dir]   <outdir>  = Specify the output directory of the library
  -s [--solve] <luafile> = Specify the LUA file with problem to solve
  -r [--run]   <luafile> = Specify the LUA file to run as a script)";

void print_help() {
    std::cout << projectHeader << std::endl;
    std::cout << helpstr << std::endl;
    exit(0);
}

void print_version() {
    std::cout << projectHeader << std::endl;
    exit(0);
}

void print_errusage(std::vector<std::string> args) {
    std::string strerr = "Unknown usage 'minos-cli";
    for (int i = 0; i < args.size(); ++i) strerr += " " + args[i];
    strerr += "'";
    strerr += "\n                 run 'minos-cli -h' for all supported options";
    throw std::runtime_error(strerr);
}

std::string rem_ext(std::string filename) {
    size_t dot_pos = filename.rfind('.');
    if (dot_pos == std::string::npos || dot_pos == 0) {
        return filename;
    }
    return filename.substr(0, dot_pos);
}

int parse_args(int argc, char* argv[], std::string& infile, std::string& outfile, std::string& outdir) {
    // convert char* argv into std::vector<std::string>
    std::vector<std::string> args(argc);
    for (int i = 0; i < argc; ++i) {
        args[i] = std::string(argv[i]); // Use emplace_back to add each argument
    }
    // too many args
    if (argc > 6) {
        print_errusage(args); // print error and throw
        return UNVALID_USAGE;
    }
    // 1 arg: -h or -v expected
    if (argc == 1) {
        if ((args[0] == "-h") || (args[0] == "--help")) { // print help and exit
            print_help(); 
        }
        else if ((args[0] == "-v") || (args[0] == "--version")) { // print version and exit
            print_version();
        }
    }
    // for now a pair of args is expected
    if ((argc % 2) != 0) { // odd num of args
        print_errusage(args); // print error and throw
        return UNVALID_USAGE;
    }
    // solve usage: -s <luafile.lua>
    if ((args[0] == "-s") || (args[0] == "--solve")) {
        infile = args[1];
        if (argc > 2) {    
            print_errusage(args); // print error and throw
            return UNVALID_USAGE;
        }
        return SOLVE_USAGE;
    }
    // run usage: -s <luafile.lua>
    if ((args[0] == "-r") || (args[0] == "--run")) {
        infile = args[1];
        if (argc > 2) {    
            print_errusage(args); // print error and throw
            return UNVALID_USAGE;
        }
        return RUN_USAGE;
    }
    // build usage: -b <ocp.h> [-o <outname> -d <outdir>]
    if ((args[0] == "-b") || (args[0] == "--build")) {
        infile = args[1];
    } else {
        print_errusage(args); // print error and throw
        return UNVALID_USAGE;
    }
    // init optional flags
    outdir = "."; // current dir
    outfile = rem_ext(infile); 
    // check other build args
    int flago = 1, flagd = 1; // 0 when flag consumed
    for (int i = 2; i < 6; i+=2) {
        if (argc > i) {
            if (flago && ((args[i] == "-o") || (args[i] == "--out"))) {
                flago = 0; // consume flag
                outfile = args[i+1];
            } else if (flagd && ((args[i] == "-d") || (args[i] == "--dir"))) {
                flagd = 0; // consume flag
                outdir = args[i+1];
            } else {
                print_errusage(args); // print error and throw
                return UNVALID_USAGE;
            }
        }
    }
    // fix end char of outdir
    if (outdir.back() == '/') {
        outdir.pop_back();
    }
    return BUILD_USAGE;
}


bool exist_file(const std::string &filename) {
    std::ifstream infile(filename);
    return infile.good();
}

int create_dir(const std::string& dir) {
    if (dir == ".") {
        return 0;
    }
#ifdef _WIN32
    // Check if the directory exists
    DWORD ftyp = GetFileAttributesA(dir.c_str());
    if (ftyp == INVALID_FILE_ATTRIBUTES) {
        // Directory does not exist, create it
        if (CreateDirectoryA(dir.c_str(), NULL) || GetLastError() == ERROR_ALREADY_EXISTS) {
            return 0; // Directory created or already exists
        }
    }
#else
    // Check if the directory exists
    struct stat sb;
    if (stat(dir.c_str(), &sb) == -1) {
        // Directory does not exist, create it
        if (mkdir(dir.c_str(), 0755) == 0) {
            return 0; // Directory created
        }
    }
#endif
    return 1; // Failed to create or directory already exists
}

int exec_cmd(const std::string &cmd, bool display_out) {
    std::string full_cmd = cmd;
    if (display_out)  {
        std::cout << cmd << std::endl;
    } else {
#ifdef _WIN32
        full_cmd = cmd + " > NUL 2>&1";
#else
        full_cmd = cmd + " > /dev/null 2>&1";
#endif
    }
    FILE* pipe = popen(full_cmd.c_str(), "r");
    return pclose(pipe);
}

std::string get_env(const std::string &name) {
#ifdef _WIN32
    char* buf = nullptr;
    size_t sz = 0;
    if (_dupenv_s(&buf, &sz, name.c_str()) == 0 && buf != nullptr) {
        std::string result(buf);
        free(buf);
        return result;
    }
    return "";
#else
    const char* val = getenv(name.c_str());
    return val == nullptr ? "" : std::string(val);
#endif
}

void set_env(const std::string &name, const std::string &value) {
#ifdef _WIN32
    _putenv_s(name.c_str(), value.c_str());
#else
    setenv(name.c_str(), value.c_str(), 1);
#endif
}

std::vector<std::string> split_path(const std::string& path, const char delimiter) {
    std::vector<std::string> paths;
    size_t start = 0;
    size_t end = path.find(delimiter);
    while (end != std::string::npos) {
        paths.push_back(path.substr(start, end - start));
        start = end + 1;
        end = path.find(delimiter, start);
    }
    paths.push_back(path.substr(start)); // Add the last path
    return paths;
}

int find_executable(const std::string& exec_name, std::string& exec_path) {
    std::string path_env =  get_env("PATH");
#ifdef _WIN32
    const char delimiter = ';';
#else
    const char delimiter = ':';
#endif
    std::vector<std::string> paths = split_path(path_env, delimiter);
    for (const std::string& dir : paths) {
        std::string full_path = dir + "/" + exec_name;
        if (exist_file(full_path)) {
            exec_path = full_path;
            return 0;
        }
    }
    return 1;
}

void reset_path(const std::string &oldpath) {
    set_env("PATH", oldpath);
}

std::vector<double> mat_to_vec(const std::vector<std::vector<double>>& mat) {
    // Check if matrix is non-empty
    if (mat.empty()) {
        return {};
    }
    // Result vector to hold the rolled-out matrix
    std::vector<double> vec;
    // Get number of rows and columns
    size_t rows = mat.size();
    size_t cols = mat[0].size();
    // Loop through each column
    for (size_t col = 0; col < cols; ++col) {
        // Loop through each row in the current column
        for (size_t row = 0; row < rows; ++row) {
            vec.push_back(mat[row][col]);
        }
    }
    return vec;
}

std::vector<std::vector<double>> vec_to_mat(const std::vector<double>& vec, size_t n, size_t m) {
    // Result matrix to hold the reshaped data
    std::vector<std::vector<double>> mat(n, std::vector<double>(m));
    // Fill the matrix by iterating over the vector
    for (size_t col = 0; col < m; ++col) {
        for (size_t row = 0; row < n; ++row) {
            mat[row][col] = vec[col * n + row];
        }
    }
    return mat;
}


std::vector<std::vector<double>> read_luamatrix(lua_State* L, const char* key) {
    std::vector<std::vector<double>> mat;
    // Push the key (as a string) to access the table
    lua_pushstring(L, key);
    lua_gettable(L, -2); // Get the table associated with the key
    // Check if the retrieved value is a table
    if (!lua_istable(L, -1)) {
        std::ostringstream strerr;
        strerr << "Expecting table for '" << key << "'";
        lua_pop(L, 1); // Pop the invalid value
        throw std::runtime_error(strerr.str());
        return mat;
    }
    // Loop through the outer table (the rows)
    for (int i = 1; true; ++i) { // Lua is 1-indexed
        lua_pushinteger(L, i); // Push the row index
        lua_gettable(L, -2); // Get the row table (mat[i])
        // Check if the row is nil, indicating no more rows
        if (lua_isnil(L, -1)) {
            lua_pop(L, 1); // Pop nil
            break; // No more rows
        }
        // Ensure the row is a table
        if (!lua_istable(L, -1)) {
            std::ostringstream strerr;
            strerr << "Expecting table for '" << key << "[" << i << "]'";
            lua_pop(L, 1); // Pop the invalid row
            throw std::runtime_error(strerr.str());
            return mat;
        }
        std::vector<double> row;
        for (int j = 1; true; ++j) { // Loop through the columns
            lua_pushinteger(L, j); // Push the column index
            lua_gettable(L, -2); // Get mat[i][j]
            // Check if the column value is nil, indicating no more columns
            if (lua_isnil(L, -1)) {
                lua_pop(L, 1); // Pop nil
                break; // No more columns
            }
            // Check if the value is a number
            if (lua_isnumber(L, -1)) {
                row.push_back(lua_tonumber(L, -1)); // Add value to the row
            } else {
                std::ostringstream strerr;
                strerr << "Expecting number for '" << key << "[" << i << "][" << j << "]'";
                lua_pop(L, 1); // Pop invalid value
                throw std::runtime_error(strerr.str());
                return mat;
            }
            lua_pop(L, 1); // Pop the column value
        }
        mat.push_back(row); // Add the row to the matrix
        lua_pop(L, 1); // Pop the row table
    }
    lua_pop(L, 1); // Pop the table associated with the key
    return mat; // Return the constructed matrix
}

std::vector<double> read_luavector(lua_State* L, const char* key) {
    std::vector<double> vec;
    // Push the key (as a string) to access the table
    lua_pushstring(L, key);
    lua_gettable(L, -2); // Get the table associated with the key
    // Check if the retrieved value is a table
    if (!lua_istable(L, -1)) {
        std::ostringstream strerr;
        strerr << "Expecting table for '" << key << "'";
        lua_pop(L, 1); // Pop the invalid value
        throw std::runtime_error(strerr.str());
        return vec; 
    }
    // Loop through the table (1D vector)
    for (int i = 1; true; ++i) { // Lua is 1-indexed
        lua_pushinteger(L, i); // Push the index
        lua_gettable(L, -2); // Get the element vec[i]
        // Check if the value is nil, indicating no more elements
        if (lua_isnil(L, -1)) {
            lua_pop(L, 1); // Pop nil
            break; // No more elements
        }
        // Check if the value is a number
        if (lua_isnumber(L, -1)) {
            vec.push_back(lua_tonumber(L, -1)); // Add the value to the vector
        } else {
            std::ostringstream strerr;
            strerr << "Expecting number for '" << key << "[" << i << "]'" << std::endl;
            lua_pop(L, 1); // Pop the invalid value
            throw std::runtime_error(strerr.str());
        return vec; 
        }
        lua_pop(L, 1); // Pop the element value
    }
    lua_pop(L, 1); // Pop the table associated with the key
    return vec; // Return the constructed vector
}

std::string read_luastring(lua_State* L, int index, const char* key) {
    std::string str;
    lua_getfield(L, index, key); // Get the string associated with the key
    if (!lua_isstring(L, -1)) {
        std::ostringstream strerr;
        strerr << "Expecting string for '" << key << "'";
        lua_pop(L, 1); // Pop the invalid value
        throw std::runtime_error(strerr.str());
        return str;
    }
    str = lua_tostring(L, -1);
    lua_pop(L, 1);
    return str;
}

double read_luanumber(lua_State* L, int index, const char* key) {
    double val;
    lua_getfield(L, index, key); // Get the number associated with the key
    if (!lua_isnumber(L, -1)) {
        std::ostringstream strerr;
        strerr << "Expecting number for '" << key << "'";
        lua_pop(L, 1); // Pop the invalid value
        throw std::runtime_error(strerr.str());
        return val;
    }
    val = lua_tonumber(L, -1);
    lua_pop(L, 1);
    return val;
}

int read_luainteger(lua_State* L, int index, const char* key) {
    int val;
    lua_getfield(L, index, key); // Get the number associated with the key
    if (!lua_isinteger(L, -1)) {
        std::ostringstream strerr;
        strerr << "Expecting integer for '" << key << "'";
        lua_pop(L, 1); // Pop the invalid value
        throw std::runtime_error(strerr.str());
        return val;
    }
    val = lua_tointeger(L, -1);
    lua_pop(L, 1);
    return val;
}

bool read_luabool(lua_State* L, int index, const char* key) {
    bool val;
    lua_getfield(L, index, key); // Get the value associated with the key
    if (!lua_isboolean(L, -1)) { // Check if it's a boolean
        std::ostringstream strerr;
        strerr << "Expecting boolean for '" << key << "'";
        lua_pop(L, 1); // Pop the invalid value
        throw std::runtime_error(strerr.str());
        return val;
    }
    val = lua_toboolean(L, -1); // Retrieve the boolean value
    lua_pop(L, 1); // Pop the value from the Lua stack
    return val;
}

bool exist_luakey(lua_State* L, const char* key) {
    lua_pushstring(L, key); // Push the key onto the stack
    lua_gettable(L, -2); // Get the value from the table using the key
    bool exists = !lua_isnil(L, -1); // If the value is not nil, the key exists
    lua_pop(L, 1); // Pop the value (whether nil or not) off the stack
    return exists;
}

bool check_vector_size(std::vector<double> vec, size_t n, const char* key) {
    if (vec.size() != n) {
        std::ostringstream strerr;
        strerr << "Expecting " << n << " elements for '" << key << "' (found " << vec.size() << ")";
        throw std::runtime_error(strerr.str());
        return false;
    }
    return true;
}

bool check_matrix_size(std::vector<std::vector<double>> mat, size_t n, size_t m, const char* key) {
    if (mat.size() != n) {
        std::ostringstream strerr;
        strerr << "Expecting " << n << " elements for '" << key << "' (found " << mat.size() << ")";
        throw std::runtime_error(strerr.str());
        return false;
    }
    for (int i = 0; i < n; ++i) {
        if (mat[i].size() != m) {
            std::ostringstream strerr;
            strerr << "Expecting " << m << " elements for '" << key << "[" << (i + 1) << "]' (found " << mat[i].size() << ")"; // one-based
            throw std::runtime_error(strerr.str());
            return false;
        }
    }
    return true;
}

int read_luaproblem(lua_State *L, int index, OCPInterface **ocp, std::string &outfile) {
    // Declarations
    std::string name;
    int N, nx, nu, np, nc, nb, nq, na;
    double ti, tf;
    std::vector<std::vector<double>> x, u, lamx, lamu, lamf, lamc;
    std::vector<double> p, lamp, lamb, lamq, lbx, ubx, lbu, ubu, lbp, ubp,
                        lbc, ubc, lbb, ubb, lbq, ubq, mesh, auxdata,
                        xv, uv, lamxv, lamuv, lamfv, lamcv;
    std::string nlpsolver, logfile;
    int flag_hessian = -1, max_iter = -1, print_itersol = -1, display = -1;
    double mu_init = -1;
    // Get problem definitions
    name = read_luastring(L, index, "name");
    N = read_luainteger(L, index, "N");
    ti = read_luanumber(L, index, "ti");
    tf = read_luanumber(L, index, "tf");
    // Create OCP object
    try {
        (*ocp) = new OCPInterface(name, N, ti, tf);
    } catch (const std::exception &e) {
        throw e;
        return 1;
    }
    // Get OCP dimensions
    (*ocp)->get_dims(&nx, &nu, &np, &nc, &nb, &nq, NULL, NULL, NULL, NULL, &na);
    // Get guess
    if (!exist_luakey(L, "guess")) {
        throw std::runtime_error("Expecting table for 'guess'");
        return 1;
    }
    lua_getfield(L, index, "guess");
    x = read_luamatrix(L, "x");
    u = read_luamatrix(L, "u");
    p = read_luavector(L, "p");
    // Retrieve optional guess
    if (exist_luakey(L, "lam_x")) lamx = read_luamatrix(L, "lam_x");
    if (exist_luakey(L, "lam_u")) lamu = read_luamatrix(L, "lam_u");
    if (exist_luakey(L, "lam_p")) lamp = read_luavector(L, "lam_p");
    if (exist_luakey(L, "lam_f")) lamf = read_luamatrix(L, "lam_f");
    if (exist_luakey(L, "lam_c")) lamc = read_luamatrix(L, "lam_c");
    if (exist_luakey(L, "lam_b")) lamb = read_luavector(L, "lam_b");
    if (exist_luakey(L, "lam_q")) lamq = read_luavector(L, "lam_q");
    lua_pop(L, 1); // pop guess
    // Check sizes
    if (!check_matrix_size(x, nx, N, "x")) { return 1; }
    if (!check_matrix_size(u, nu, N, "u")) { return 1; }
    if (!check_vector_size(p, np, "p")) { return 1; }
    if ((lamx.size()) && !check_matrix_size(lamx, nx, N, "lam_x")) { return 1; }
    if ((lamu.size()) && !check_matrix_size(lamu, nu, N, "lam_u")) { return 1; }
    if ((lamp.size()) && !check_vector_size(lamp, np, "lam_p")) { return 1; }
    if ((lamf.size()) && !check_matrix_size(lamf, nx, N-1, "lam_f")) { return 1; }
    if ((lamc.size()) && !check_matrix_size(lamc, nc, N, "lam_c")) { return 1; }
    if ((lamb.size()) && !check_vector_size(lamb, nb, "lam_b")) { return 1; }
    if ((lamq.size()) && !check_vector_size(lamq, nq, "lam_q")) { return 1; }
    // Transfor matrices into a unique vector by rolling out by columns
    xv = mat_to_vec(x);
    uv = mat_to_vec(u);
    lamxv = mat_to_vec(lamx);
    lamuv = mat_to_vec(lamu);
    lamfv = mat_to_vec(lamf);
    lamcv = mat_to_vec(lamc);
    // Get bounds
    if (!exist_luakey(L, "bounds")) {
        throw std::runtime_error("Expecting table for 'bounds'");
        return 1;
    }
    lua_getfield(L, index, "bounds");
    lbx = read_luavector(L, "lbx"); ubx = read_luavector(L, "ubx");
    lbu = read_luavector(L, "lbu"); ubu = read_luavector(L, "ubu");
    lbp = read_luavector(L, "lbp"); ubp = read_luavector(L, "ubp");
    lbc = read_luavector(L, "lbc"); ubc = read_luavector(L, "ubc");
    lbb = read_luavector(L, "lbb"); ubb = read_luavector(L, "ubb");
    lbq = read_luavector(L, "lbq"); ubq = read_luavector(L, "ubq");
    lua_pop(L, 1); // pop bounds
    // Check sizes
    if (!check_vector_size(lbx, nx, "lbx") || !check_vector_size(ubx, nx, "ubx")) { return 1; }
    if (!check_vector_size(lbu, nu, "lbu") || !check_vector_size(ubu, nu, "ubu")) { return 1; }
    if (!check_vector_size(lbp, np, "lbp") || !check_vector_size(ubp, np, "ubp")) { return 1; }
    if (!check_vector_size(lbc, nc, "lbc") || !check_vector_size(ubc, nc, "ubc")) { return 1; }
    if (!check_vector_size(lbb, nb, "lbb") || !check_vector_size(ubb, nb, "ubb")) { return 1; }
    if (!check_vector_size(lbq, nq, "lbq") || !check_vector_size(ubq, nq, "ubq")) { return 1; }
    // Get auxdata
    if (exist_luakey(L, "auxdata")) {
        auxdata = read_luavector(L, "auxdata");
        if (!check_vector_size(auxdata, na, "auxdata")) {return 1; }
    }
    // Get mesh
    if (exist_luakey(L, "mesh")) {
        mesh = read_luavector(L, "mesh");
        if (!check_vector_size(mesh, N-1, "mesh")) {return 1; }
        // check sum to unity
        double sum = 0;
        double eps = 1.0e-9;
        for (int i = 0; i < (N-1); ++i) sum += mesh[i];
        if (((sum-1) > eps) || ((sum-1) < -eps)) {
            std::ostringstream strerr;
            char sumstr[256];
            sprintf(sumstr, "%.9f", sum);
            strerr << "Expecting sum to 1 for 'mesh' (sum to " << sumstr << ")";
            throw std::runtime_error(strerr.str());
            return 1;
        }
    }
    // Get options
    if (exist_luakey(L, "options")) {
        lua_getfield(L, index, "options");
        if (exist_luakey(L, "nlpsolver"))
            nlpsolver = read_luastring(L, -1, "nlpsolver");
        if (exist_luakey(L, "flag_hessian"))
            flag_hessian = (int) read_luabool(L, -1, "flag_hessian");
        if (exist_luakey(L, "max_iter"))
            max_iter = read_luainteger(L, -1, "max_iter");
        if (exist_luakey(L, "mu_init"))
            mu_init = read_luanumber(L, -1, "mu_init");
        if (exist_luakey(L, "logfile"))
            logfile = read_luastring(L, -1, "logfile");
        if (exist_luakey(L, "outfile"))
            outfile = read_luastring(L, -1, "outfile");
        if (exist_luakey(L, "print_itersol"))
            print_itersol = read_luainteger(L, -1, "print_itersol");
        if (exist_luakey(L, "display"))
            display = (int) read_luabool(L, -1, "display");
        lua_pop(L, 1); // pop options   
    }
    // Set OCP
    // Set bounds
    double* lbxp = (lbx.size()) ? lbx.data() : NULL;
    double* ubxp = (ubx.size()) ? ubx.data() : NULL;
    double* lbup = (lbu.size()) ? lbu.data() : NULL;
    double* ubup = (ubu.size()) ? ubu.data() : NULL;
    double* lbpp = (lbp.size()) ? lbp.data() : NULL;
    double* ubpp = (ubp.size()) ? ubp.data() : NULL;
    double* lbcp = (lbc.size()) ? lbc.data() : NULL;
    double* ubcp = (ubc.size()) ? ubc.data() : NULL;
    double* lbbp = (lbb.size()) ? lbb.data() : NULL;
    double* ubbp = (ubb.size()) ? ubb.data() : NULL;
    double* lbqp = (lbq.size()) ? lbq.data() : NULL;
    double* ubqp = (ubq.size()) ? ubq.data() : NULL;
    (*ocp)->set_bounds(lbxp, ubxp, lbup, ubup, lbpp, ubpp, lbcp, ubcp, lbbp, ubbp, lbqp, ubqp);
    // Set guess
    double* xp = (x.size()) ? xv.data() : NULL;
    double* up = (u.size()) ? uv.data() : NULL;
    double* pp = (p.size()) ? p.data() : NULL;
    double* lampp = (lamp.size()) ? lamp.data() : NULL;
    double* lambp = (lamb.size()) ? lamb.data() : NULL;
    double* lamqp = (lamq.size()) ? lamq.data() : NULL;
    double* lamxp = (lamx.size()) ? lamxv.data() : NULL;
    double* lamup = (lamu.size()) ? lamuv.data() : NULL;
    double* lamfp = (lamf.size()) ? lamfv.data() : NULL;
    double* lamcp = (lamc.size()) ? lamcv.data() : NULL;
    (*ocp)->set_guess(xp, up, pp, lamxp, lamup, lampp, lamfp, lamcp, lambp, lamqp);
    // Set auxdata
    if (auxdata.size()) (*ocp)->set_auxdata(auxdata.data());
    // Set mesh
    if (mesh.size()) (*ocp)->set_mesh(mesh.data());
    // Set options
    if (nlpsolver.size())   (*ocp)->set_option(OCPInterface::NLPSOLVER,     nlpsolver);
    if (logfile.size())     (*ocp)->set_option(OCPInterface::LOGFILE,       logfile);
    if (flag_hessian >=0 )  (*ocp)->set_option(OCPInterface::FLAG_HESSIAN,  flag_hessian);
    if (max_iter > 0)       (*ocp)->set_option(OCPInterface::MAX_ITER,      max_iter);
    if (mu_init > 0)        (*ocp)->set_option(OCPInterface::MU_INIT,       mu_init);
    if (print_itersol >= 0) (*ocp)->set_option(OCPInterface::PRINT_ITERSOL, print_itersol);
    if (display >= 0)       (*ocp)->set_option(OCPInterface::DISPLAY,       display);
    // Return
    return 0;
}

void write_luamatrix(lua_State *L, const char *key, std::vector<std::vector<double>> mat) {
    lua_pushstring(L, key);
    lua_newtable(L);
    for (size_t i = 0; i < mat.size(); ++i) {
        lua_pushnumber(L, i + 1); // one-based
        lua_newtable(L);
        for (size_t j = 0; j < mat[i].size(); ++j) {
            lua_pushnumber(L, j + 1); // one-based
            lua_pushnumber(L, mat[i][j]);
            lua_settable(L, -3); // set row[j+1]=mat[i][j]
        }
        lua_settable(L, -3); // set table[i+1] = row
    }
    lua_settable(L, -3); // set main[key] = table
}

void write_luavector(lua_State *L, const char *key, std::vector<double> vec) {
    lua_pushstring(L, key);
    lua_newtable(L);
    for (size_t i = 0; i < vec.size(); ++i) {
        lua_pushnumber(L, i + 1); // one-based
        lua_pushnumber(L, vec[i]);
        lua_settable(L, -3); // set table[i+1]=vec[i]
    }
    lua_settable(L, -3); // set main[key] = table
}

void write_luanumber(lua_State *L, const char *key, double val) {
    lua_pushstring(L, key);
    lua_pushnumber(L, val);
    lua_settable(L, -3); // set main[key] = val
}

void write_luastring(lua_State *L, const char *key, std::string str) {
    lua_pushstring(L, key);
    lua_pushstring(L, str.c_str());
    lua_settable(L, -3); // set main[key] = str
}

void write_luabool(lua_State *L, const char *key, bool val) {
    lua_pushstring(L, key);
    lua_pushboolean(L, val);
    lua_settable(L, -3); // set main[key] = val
}

void write_luainteger(lua_State *L, const char *key, int val) {
    lua_pushstring(L, key);
    lua_pushinteger(L, val);
    lua_settable(L, -3); // set main[key] = val
}

void write_luasolution(lua_State *L, OCPInterface* ocp, std::string outfile) {
    // Init variables
    int nx, nu, np, nc, nb, nq, nz, ng, nnzj, nnzh, na, N, num_iter;
    std::string name, nlpsolver, logfile;
    double ti, tf, objval, m, ttot, talg, teval, mu_curr, flag_hessian, print_itersol, max_iter, display;
    // Get
    name = ocp->get_name();
    N = ocp->get_N();
    ocp->get_t(&ti, &tf);
    ocp->get_dims(&nx, &nu, &np, &nc, &nb, &nq, &nz, &ng, &nnzj, &nnzh, &na);
    num_iter = ocp->get_num_iter();
    // Init variables
    std::vector<double> t(N), x(nx*N), u(nu*N), p(np),
                        lamx(nx*N), lamu(nu*N), lamp(np),
                        lamf(nx*(N-1)), lamc(nc*N), lamb(nb), lamq(nq), 
                        f(nx*(N-1)), c(nc*N), b(nb), q(nq), l(N-1);
    std::vector<double> obj_history(num_iter+1), infpr_history(num_iter+1), infdu_history(num_iter+1);
    std::vector<double> lbx(nx), ubx(nx), lbu(nu), ubu(nu), lbp(np), ubp(np),
                        lbc(nc), ubc(nc), lbb(nb), ubb(nb), lbq(nq), ubq(nq);
    std::vector<double> auxdata(na), mesh(N-1);
    // Get solution
    ocp->get_sol(&objval, t.data(),
            x.data(), u.data(), p.data(),
            lamx.data(), lamu.data(), lamp.data(),
            lamf.data(), lamc.data(), lamb.data(), lamq.data(),
            f.data(), c.data(), b.data(), q.data(), 
            l.data(), &m);
    ocp->get_cpu_time(ttot, talg, teval);
    ocp->get_history(obj_history.data(), infpr_history.data(), infdu_history.data());
    mu_curr = ocp->get_mu_curr();
    ocp->get_bounds(lbx.data(), ubx.data(), lbu.data(), ubu.data(), lbp.data(), ubp.data(), 
                    lbc.data(), ubc.data(), lbb.data(), ubb.data(), lbq.data(), ubq.data());
    ocp->get_auxdata(auxdata.data());
    ocp->get_mesh(mesh.data());
    ocp->get_option(OCPInterface::NLPSOLVER, nlpsolver);
    ocp->get_option(OCPInterface::LOGFILE, logfile);
    ocp->get_option(OCPInterface::FLAG_HESSIAN, &flag_hessian);
    ocp->get_option(OCPInterface::PRINT_ITERSOL, &print_itersol);
    ocp->get_option(OCPInterface::MAX_ITER, &max_iter);
    ocp->get_option(OCPInterface::DISPLAY, &display);
    // Instantiate new empty table
    lua_newtable(L);
    // Write solution
    write_luanumber(L, "objval", objval);
    write_luavector(L, "t", t);
    write_luamatrix(L, "x", vec_to_mat(x, nx, N));
    write_luamatrix(L, "u", vec_to_mat(u, nu, N));
    write_luavector(L, "p", p);
    write_luamatrix(L, "lam_x", vec_to_mat(lamx, nx, N));
    write_luamatrix(L, "lam_u", vec_to_mat(lamu, nu, N));
    write_luavector(L, "lam_p", lamp);
    write_luamatrix(L, "lam_f", vec_to_mat(lamf, nx, N-1));
    write_luamatrix(L, "lam_c", vec_to_mat(lamc, nc, N));
    write_luavector(L, "lam_b", lamb);
    write_luavector(L, "lam_q", lamq);
    write_luamatrix(L, "f", vec_to_mat(f, nx, N-1));
    write_luamatrix(L, "c", vec_to_mat(c, nc, N));
    write_luavector(L, "b", b);
    write_luavector(L, "q", q);
    write_luavector(L, "l", l);
    write_luanumber(L, "m", m);
    // Write stats
    lua_pushstring(L, "stats");
    lua_newtable(L);
    write_luanumber(L, "num_iter", num_iter);
    write_luavector(L, "obj_history", obj_history);
    write_luavector(L, "infpr_history", infpr_history);
    write_luavector(L, "infdu_history", infdu_history);
    write_luanumber(L, "mu_curr", mu_curr);
    write_luanumber(L, "nz", nz);
    write_luanumber(L, "ng", ng);
    write_luanumber(L, "nnzj", nnzj);
    write_luanumber(L, "nnzh", nnzh);
    write_luanumber(L, "ttot", ttot);
    write_luanumber(L, "talg", talg);
    write_luanumber(L, "teval", teval);
    lua_settable(L, -3); // set main["stats"] = stats
    // Write next_problem
    lua_pushstring(L, "next_problem");
    lua_newtable(L);
    write_luastring(L, "name", name);
    write_luainteger(L, "N", N);
    write_luanumber(L, "ti", ti);
    write_luanumber(L, "tf", tf);
    // guess table
    lua_pushstring(L, "guess");
    lua_newtable(L);
    write_luamatrix(L, "x", vec_to_mat(x, nx, N));
    write_luamatrix(L, "u", vec_to_mat(u, nu, N));
    write_luavector(L, "p", p);
    write_luamatrix(L, "lam_x", vec_to_mat(lamx, nx, N));
    write_luamatrix(L, "lam_u", vec_to_mat(lamu, nu, N));
    write_luavector(L, "lam_p", lamp);
    write_luamatrix(L, "lam_f", vec_to_mat(lamf, nx, N-1));
    write_luamatrix(L, "lam_c", vec_to_mat(lamc, nc, N));
    write_luavector(L, "lam_b", lamb);
    write_luavector(L, "lam_q", lamq);
    lua_settable(L, -3); // set next_problem["guess"] = guess
    // bound table
    lua_pushstring(L, "bounds");
    lua_newtable(L);
    write_luavector(L, "lbx", lbx);
    write_luavector(L, "ubx", ubx);
    write_luavector(L, "lbu", lbu);
    write_luavector(L, "ubu", ubu);
    write_luavector(L, "lbp", lbp);
    write_luavector(L, "ubp", ubp);
    write_luavector(L, "lbc", lbc);
    write_luavector(L, "ubc", ubc);
    write_luavector(L, "lbb", lbb);
    write_luavector(L, "ubb", ubb);
    write_luavector(L, "lbq", lbq);
    write_luavector(L, "ubq", ubq);
    lua_settable(L, -3); // set next_problem["bounds"] = bounds
    // auxdata table
    write_luavector(L, "auxdata", auxdata);
    // mesh table
    write_luavector(L, "mesh", mesh);
    // options table
    lua_pushstring(L, "options");
    lua_newtable(L);
    write_luastring(L, "nlpsolver", nlpsolver);
    write_luastring(L, "logfile", logfile);
    write_luanumber(L, "mu_init", mu_curr);
    write_luainteger(L, "max_iter", max_iter);
    write_luainteger(L, "print_itersol", print_itersol);
    write_luabool(L, "flag_hessian", (bool) flag_hessian);
    write_luabool(L, "display", (bool) display);
    if (outfile.size()) {
        write_luastring(L, "outfile", outfile);
    }
    lua_settable(L, -3); // set next_problem["options"] = options
    lua_settable(L, -3); // set main["next_problem"] = next_problem
}

int build(std::string cc, std::string csource, std::string outfile, 
    std::string outdir) {
    // Check if C file exists
    if (!exist_file(csource)) {
        throw std::runtime_error("Unable to find source C file: " + csource);
        return 1;
    }
    // Get current PATH
    std::string oldpath = get_env("PATH");
    // Test the C compiler
    std::string output;
    int exitcode = exec_cmd(cc + " --version", false);
#ifdef _WIN32
    // If on Windows and global GCC found
    if (exitcode == 0) {
        std::string ccpath;
        exitcode = find_executable(cc + ".exe", ccpath);
        if (exitcode == 0) {
            std::cout << "Using user C compiler at '" << ccpath << "'\n";
        }
    }
    // If global GCC not found, try local GCC distribution
    if (exitcode != 0) {
        char buffer[MAX_PATH];
        // Get the path of the executable
        GetModuleFileNameA(NULL, buffer, MAX_PATH);
        std::string basedir(buffer);
        size_t pos = basedir.find_last_of("\\/");
        if (pos != std::string::npos) {
            basedir = basedir.substr(0, pos) + "/";
        }
        set_env("PATH", basedir + "gcc/bin;" + oldpath);
        exitcode = exec_cmd(cc + " --version", false);
    }
#endif
    // Throw error if GCC is still not found
    if (exitcode != 0) {
        reset_path(oldpath);
        throw std::runtime_error("Unable to find C compiler: " + cc);
        return 1;
    }
    // Creatr outdir if not exists
    exitcode = create_dir(outdir);
    if (exitcode != 0) {
        reset_path(oldpath);
        throw std::runtime_error("Failed to create directory: " + outdir);
        return 1;
    }
    // Set the library extension based on the platform
#ifdef _WIN32
    std::string libext = "dll";
#else
    std::string libext = "so";
#endif
    std::string libname = outdir + "/" + outfile + "." + libext;
    // Build command
    std::string cc_cmd = cc + " -shared -O1 -fPIC \"" + csource + "\" -o \"" + libname + "\"";
    // Run the build process
    auto start = std::chrono::high_resolution_clock::now();
    std::string build_out;
    exitcode = exec_cmd(cc_cmd, true);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> compileTime = end - start;
    // Check exit status and reset PATH if needed
    if (exitcode != 0) {
        reset_path(oldpath);
        throw std::runtime_error("Unable to build library '" + libname + "' from file '" + csource + "'\n" + build_out);
        return 1;
    }
    // Print end message
    std::cout << "Library " << libname << " built in " << compileTime.count() << "s\n";
    return 0;
}

int solve(const std::string &luafile) {
    // Start LUA
    lua_State* L = luaL_newstate();    // Create a new Lua state
    luaL_openlibs(L);                  // Load Lua libraries
    // Load LUA file
    if (luaL_dofile(L, luafile.c_str())) {
        std::string strerr = lua_tostring(L,-1);
        lua_close(L);
        throw std::runtime_error("Failed to load LUA file: " + strerr);
        return 1;
    }
    // Read problem and allocate OCP
    OCPInterface *ocp = NULL;
    std::string outfile;
    //Check problem type
    lua_getglobal(L, "problem"); // push problem on the top of the stack
    if (!lua_istable(L, -1)) {
        throw std::runtime_error("Expecting table for 'problem'");
        return 1;
    }
    // Parse solve input
    try {
        read_luaproblem(L, -1, &ocp, outfile); // read problem
    } catch (const std::exception &e) {
        if (ocp) delete ocp;
        throw e;
        return 1;
    }
    // Solve OCP
    try {
        ocp->solve();
    } catch (const std::exception &e) {
        delete ocp;
        throw e;
        return 1;
    }
    // Write outfile if required
    if (outfile.size()) {
        std::ofstream outstream;
        outstream.open(outfile);
        if (!outstream.is_open()) {
            delete ocp;
            throw std::runtime_error("Failed to open outfile file: " + outfile);
        return 1;
        }
        outstream << *ocp;
        outstream.close();
    }
    // Free mem and return
    lua_close(L);
    delete ocp;
    return 0;
}

int run_build(lua_State *L) {
    // Check first argument: infile
    std::string infile = luaL_checkstring(L, 1);
    // Check second argument: outfile
    std::string outfile = luaL_optstring(L, 2, rem_ext(infile).c_str());
    // Check third argument: outdir
    std::string outdir = luaL_optstring(L, 3, ".");
    // Call to build
    try {
        build("gcc", infile, outfile, outdir);
    } catch (const std::exception &e) {
        return luaL_error(L, "%s", e.what()); // pass error to LUA
    }
    return 0;
}

int run_solve(lua_State *L) {
    OCPInterface *ocp = NULL;
    std::string outfile;
    // Check input argument
    luaL_checktype(L, 1, LUA_TTABLE);
    // Parse solve input
    try {
        read_luaproblem(L, 1, &ocp, outfile);
    } catch (const std::exception &e) {
        if (ocp) delete ocp;
        return luaL_error(L, "%s", e.what()); // pass error to LUA
    }
    // Call to solve
    try {
        ocp->solve();
    } catch (const std::exception &e) {
        delete ocp;
        return luaL_error(L, "%s", e.what()); // pass error to LUA
    }
    // Write solution
    write_luasolution(L, ocp, outfile);
    // Write outfile if required
    if (outfile.size()) {
        std::ofstream outstream;
        outstream.open(outfile);
        if (!outstream.is_open()) {
            delete ocp;
            return luaL_error(L, "Failed to open outfile file: %s", outfile.c_str()); // pass error to LUA
        }
        outstream << *ocp;
        outstream.close();
    }
    // Free mem
    delete ocp;
    // Return num of outputs
    return 1;
}

int run(const std::string &luafile) {
    // Start LUA
    lua_State* L = luaL_newstate();    // Create a new Lua state
    luaL_openlibs(L);                  // Load Lua libraries
    // Register the C function
    lua_register(L, "build", run_build);
    lua_register(L, "solve", run_solve);
    // Run LUA file
    if (luaL_dofile(L, luafile.c_str())) {
        std::string strerr = lua_tostring(L,-1);
        lua_close(L);
        throw std::runtime_error("Failed to run LUA file: " + strerr);
        return 1;
    }
    // Return 
    lua_close(L);
    return 0;
}

int main(int argc, char* argv[]) {
    argv++; argc--; // Consume first arg that is the filename
    // No args: print help and exit
    if (argc < 1) print_help();
    // Parse args
    std::string infile, outfile, outdir;
    int usage;
    try {
        usage = parse_args(argc, argv, infile, outfile, outdir);
    } catch (const std::exception &e) {
        std::cerr << "minos-cli error: " << e.what() << std::endl;
        return 1;
    }
    // Call function
    int exit;
    try {
        switch (usage) {
            case BUILD_USAGE:
                build("gcc", infile, outfile, outdir);
                break;
            case SOLVE_USAGE:
                solve(infile);
                break;
            case RUN_USAGE:
                run(infile);
                break;
        } 
    } catch (const std::exception &e) {
        std::cerr << "minos-cli error: " << e.what() << std::endl;
        return 1;
    }
    // Return success
    return 0;
}