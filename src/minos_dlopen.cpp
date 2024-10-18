/*  File minos_dlopen.cpp 
    OCPInterface class implementation with load library
    Copyright (C) 2024 Stefano Lovato
*/

#include "minos.h" 
#include "macros.h" // for nlp macros
#include <stdexcept>  // For std::runtime_error
#if _WIN32
#include <Windows.h>
#else
#include <dlfcn.h> //dlopen
#include <link.h>
#endif

// Macro to load library and function
#if _WIN32 // Windows
    #define LIBEXT ".dll"
    #define LOADLIBRARY(name, libhandle) name += LIBEXT;                           \
        HMODULE libhandle = LoadLibraryA(name.c_str());                            \
        if (!libhandle) { /* try with name "lib<name>" */                          \
            std::string libname = "lib" + name;                                    \
            libhandle = LoadLibraryA(libname.c_str());                             \
        }                                                                          \
        if (!libhandle) {                                                          \
            DWORD errorCode = GetLastError();                                      \
            char errorMessage[256];                                                \
            FormatMessageA(                                                        \
                FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,        \
                NULL,                                                              \
                errorCode,                                                         \
                0,                                                                 \
                errorMessage,                                                      \
                sizeof(errorMessage),                                              \
                NULL);                                                             \
            throw std::runtime_error("Failed to load DLL " + name + ": " +         \
                                    std::string(errorMessage));                    \
            return NULL;                                                           \
        }                                           
    #define LOADFUNCTION(libhandle, FuncType, functionName, doCheck)                   \
        reinterpret_cast<FuncType>(GetProcAddress((HMODULE) libhandle, functionName)); \
        if (doCheck && !GetProcAddress((HMODULE) libhandle, functionName)) {           \
            CHAR dllname[MAX_PATH];                                                    \
            GetModuleFileNameA((HMODULE) libhandle, dllname, sizeof(dllname));         \
            FreeLibrary((HMODULE) libhandle);                                          \
            throw std::runtime_error("Failed to find function " + std::string(functionName) + " in DLL " + std::string(dllname)); \
            return NULL;                                                               \
        }
    #define FREELIBRARY(libhandle) FreeLibrary((HMODULE) libhandle);
#else // Linux and MacOS
    #include <dlfcn.h>
    #ifdef __APPLE__ // MacOS
        #define LIBEXT ".dylib"
    #else // Linux
        #define LIBEXT ".so"
    #endif
    #define LOADLIBRARY(name, libhandle) name += LIBEXT;                     \
        void* libhandle = dlopen(name.c_str(), RTLD_LAZY);                   \
        if (!libhandle) { /* try with name "lib<name>" */                    \
            std::string libname = "lib" + name;                              \
            libhandle = dlopen(libname.c_str(), RTLD_LAZY);                  \
        }                                                                    \
        if (!libhandle) {                                                    \
            throw std::runtime_error("Failed to load shared library: " + name + " - " + dlerror()); \
            return NULL;                                                     \
        }
    #define LOADFUNCTION(libhandle, FuncType, functionName, doCheck)         \
        reinterpret_cast<FuncType>(dlsym(libhandle, functionName));          \
        if (doCheck && !dlsym(libhandle, functionName)) {                    \
            struct link_map *map;                                            \
            dlinfo(libhandle, RTLD_DI_LINKMAP, &map);                        \
            std::string libname = std::string(map->l_name);                  \
            dlclose(libhandle);                                              \
            throw std::runtime_error("Failed to find function " + std::string(functionName) + " in library " + libname); \
            return NULL;                                                    \
        }
    #define FREELIBRARY(libhandle) dlclose((void*) libhandle);
#endif

// Macro to load OPC functions with a given suffix
#define LOADOCPFUNCS(var, flag) (var).eval = LOADFUNCTION(libhandle, EvalFunc, #var, flag); \
                                (var).alloc = LOADFUNCTION(libhandle, AllocFunc, #var "_alloc_mem", flag); \
                                (var).free = LOADFUNCTION(libhandle, FreeFunc, #var "_free_mem", flag); \
                                (var).nin = LOADFUNCTION(libhandle, NInFunc, #var "_n_in", flag); \
                                (var).spin = LOADFUNCTION(libhandle, SpInFunc, #var "_sparsity_in", flag); \
                                (var).spout = LOADFUNCTION(libhandle, SpOutFunc, #var "_sparsity_out", flag); \
                                (var).work = LOADFUNCTION(libhandle, WorkFunc, #var "_work", flag);

/** Load the OCP library and import the related functions */
void* OCPInterface::load_ocplib(
    std::string name
) {
    // Load library
    LOADLIBRARY(name, libhandle);
    // Load functions
    LOADOCPFUNCS(ocp_dyn, true);
    LOADOCPFUNCS(ocp_path, true);
    LOADOCPFUNCS(ocp_bcs, true);
    LOADOCPFUNCS(ocp_int, true);
    LOADOCPFUNCS(ocp_runcost, true);
    LOADOCPFUNCS(ocp_bcscost, true);
    LOADOCPFUNCS(ocp_dyn_jac, true);
    LOADOCPFUNCS(ocp_path_jac, true);
    LOADOCPFUNCS(ocp_bcs_jac, true);
    LOADOCPFUNCS(ocp_int_jac, true);
    LOADOCPFUNCS(ocp_runcost_grad, true);
    LOADOCPFUNCS(ocp_bcscost_grad, true);
    LOADOCPFUNCS(ocp_hessb, false);
    LOADOCPFUNCS(ocp_hessi, false);
    // Check hessian
    if (!ocp_hessb.eval || !ocp_hessi.eval) {
        ocp_hessb.eval = NULL;
        ocp_hessi.eval = NULL;
    }
    // Ok, return true
    return (void*) libhandle;
}

/** Import the NLP library */
void* OCPInterface::load_nlplib(
    std::string name
) {
    // Load library
    LOADLIBRARY(name, libhandle);
    // Ok, return true
    return (void*) libhandle;
}

/** Import the NLP solver function */
void* OCPInterface::import_nlpsolve(
    void* libhandle
) {
    // Check if arg is NULL
    if (!libhandle) { return NULL; }
    // Load function
    void* callSolve_ptr = LOADFUNCTION(libhandle, void*, "callSolve", true);
    // Return
    return callSolve_ptr;
}

/** Free the specified library */
void OCPInterface::free_library(
    void* libhandle
) {
    if (libhandle) {
        FREELIBRARY(libhandle);
    }
}
