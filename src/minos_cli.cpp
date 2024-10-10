/* File minos_cli.cpp
    Command line interface for MINOS
    Copyright (C) 2024 Stefano Lovato
*/

#include "minos.h"
#include "minos_header.h"
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#define BUILD_USAGE 0
#define RUN_USAGE 1

const std::string helpstr = 
R"(Command-line interface usage

  minos-cli -b <ocp.c>
  minos-cli -b <ocp.c> -o <outname>
  minos-cli -b <ocp.c> -o <outname> -d <outdir>

Build the OCP library from the C file.

  minos-cli -r <luafile.lua>

Run the solver from the LUA file.

Options:
  -h                   = Print this help
  -v                   = Print the header with version
  -b <ocp.c>           = Specify the name of the C file to build
  -o <outname>         = Specify the name of the output library
  -d <outdit>          = Specify the output directory of the library
  -r <luafile.lua>     = Specify the LUA file to run the solver)";

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
    std::string strerr = "minos-cli error: unknown usage 'minos-cli";
    for (int i = 0; i < args.size(); ++i) strerr += " " + args[i];
    strerr += "'";
    std::cerr << strerr << std::endl;
    std::cerr << "minos-cli error: run 'minos-cli -h' for all supported options" << std::endl;
    exit(1);
}

int parse_args(int argc, char* argv[], std::string& luafile, std::string& cfile, std::string& outfile, std::string& outdir) {
    // convert char* argv into std::vector<std::string>
    std::vector<std::string> args(argc);
    for (int i = 0; i < argc; ++i) {
        args[i] = std::string(argv[i]); // Use emplace_back to add each argument
    }
    // 1 arg: -h or -v expected
    if (argc == 1) {
        if (args[0] == "-h") { // print help and exit
            print_help(); 
        } else if (args[0] == "-v") { // print version and exit
            print_version();
        }
    }
    // for now a pair of args is expected
    if ((argc % 2) != 0) { // odd num of args
        print_errusage(args); // print error and eixt
    }
    // run usage: -r <luafile.lua>
    if (args[0] == "-r") {
        if (argc > 2) {
            print_errusage(args); // print error and eixt
        }
        luafile = args[1];
        return RUN_USAGE;
    }
    // build usage: -b <ocp.h> [-o <outname> -d <outdir>]
    if (args[0] == "-b") {
        cfile = args[1];
    } else {
        print_errusage(args); // print error and eixt
    }
    // init
    outdir = "./"; // current dir
    outfile = cfile; // TODO remove extension
    // check other build args
    int flago = 1, flagd = 1;
    for (int i = 2; i < 6; i+=2) {
        if (argc > i) {
            if (flago && (args[i] == "-o")) {
                flago = 0; // consume flag
                outfile = args[i+1];
            } else if (flagd && (args[i] == "-d")) {
                flagd = 0; // consume flag
                outdir = args[i+1];
            } else {
                print_errusage(args); // print error and eixt
            }
        }
    }
    return BUILD_USAGE;
}

int build(std::string cfile, std::string outfile, std::string outdir) {
    // look at python to understand what to do (various check for gcc)
    // some calls to system(cmd);
    std::cout << "cfile: " << cfile << std::endl;
    std::cout << "outfile: " << outfile << std::endl;
    std::cout << "outdir: " << outdir << std::endl;
    return 0;
}

int run(std::string luafile) {
    // here parse lua or whatever scripting language selected
    // at the end create OCPInterface pointer, set various stuf, and call solver.
    // here always print outsol to file, with filename specified in lua
    std::cout << "luafile: " << luafile << std::endl;
    return 0;
}

int main(int argc, char* argv[]) {
    argv++; argc--; // consume first arg that is the filename
    // no args: print help and exit
    if (argc < 1) print_help();
    // parse args
    std::string luafile, cfile, outfile, outdir;
    int usage = parse_args(argc, argv, luafile, cfile, outfile, outdir);
    // call function
    int exit;
    switch (usage) {
        case BUILD_USAGE:
            exit = build(cfile, outfile, outdir);
            break;
        case RUN_USAGE:
            exit = run(luafile);
            break;
    }
    return exit;
}