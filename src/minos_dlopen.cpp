/*  File minos_dlopen.cpp 
    OCPInterface class implementation with load library
    Copyright (C) 2024 Stefano Lovato
*/

#define MAKE_MINOS

#include "minos.h" 
#include "macros.h" // for nlp macros
#include <stdexcept>  // For std::runtime_error
#if _WIN32
#include <Windows.h>
#else
#include <dlfcn.h> //dlopen
#endif

// Macro to load library and function
#if _WIN32 // Windows
    #define LIBEXT ".dll"
    #define LOADLIBRARY(name, libhandle) name += LIBEXT;                           \
        HMODULE libhandle = LoadLibraryA(name.c_str());                            \
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
            return false;                                                          \
        }                                           
    #define LOADFUNCTION(libhandle, FuncType, functionName, doCheck)         \
        reinterpret_cast<FuncType>(GetProcAddress(libhandle, functionName)); \
        if (doCheck && !GetProcAddress(libhandle, functionName)) {           \
            CHAR dllname[MAX_PATH] = { 0 };                                  \
            GetModuleFileNameA(libhandle, dllname, sizeof(dllname));         \
            FreeLibrary(libhandle);                                          \
            throw std::runtime_error("Failed to find function " + std::string(functionName) + " in DLL " + std::string(dllname)); \
            return false;                                                    \
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
            return false;                                                    \
        }
    #define LOADFUNCTION(libhandle, FuncType, functionName, doCheck)         \
        reinterpret_cast<FuncType>(dlsym(libhandle, functionName));          \
        if (doCheck && !dlsym(libhandle, functionName)) {                    \
            char* err = dlerror();                                           \
            dlclose(libhandle);                                              \
            throw std::runtime_error("Failed to find function " + std::string(functionName) + " in library - " + (err ? err : "")); \
            return false;                                                    \
        }
    #define FREELIBRARY(libhandle) dlclose((void*) libhandle);
#endif

// Macro to load OPC functions with a given suffix
#define LOADOCPFUNCS(suffix, FuncType)  ocp_dyn  ## suffix = LOADFUNCTION(libhandle, FuncType, "ocp_dyn" #suffix, true);  \
                                        ocp_path ## suffix = LOADFUNCTION(libhandle, FuncType, "ocp_path" #suffix, true); \
                                        ocp_bcs  ## suffix = LOADFUNCTION(libhandle, FuncType, "ocp_bcs" #suffix, true); \
                                        ocp_int  ## suffix = LOADFUNCTION(libhandle, FuncType, "ocp_int" #suffix, true); \
                                        ocp_runcost  ## suffix = LOADFUNCTION(libhandle, FuncType, "ocp_runcost" #suffix, true); \
                                        ocp_bcscost  ## suffix = LOADFUNCTION(libhandle, FuncType, "ocp_bcscost" #suffix, true); \
                                        ocp_dyn_jac  ## suffix = LOADFUNCTION(libhandle, FuncType, "ocp_dyn_jac" #suffix, true); \
                                        ocp_path_jac ## suffix = LOADFUNCTION(libhandle, FuncType, "ocp_path_jac" #suffix, true); \
                                        ocp_bcs_jac  ## suffix = LOADFUNCTION(libhandle, FuncType, "ocp_bcs_jac" #suffix, true); \
                                        ocp_int_jac  ## suffix = LOADFUNCTION(libhandle, FuncType, "ocp_int_jac" #suffix, true); \
                                        ocp_runcost_grad ## suffix = LOADFUNCTION(libhandle, FuncType, "ocp_runcost_grad" #suffix, true); \
                                        ocp_bcscost_grad ## suffix = LOADFUNCTION(libhandle, FuncType, "ocp_bcscost_grad" #suffix, true); \
                                        ocp_hessb ## suffix = LOADFUNCTION(libhandle, FuncType, "ocp_hessb" #suffix, false); \
                                        ocp_hessi ## suffix = LOADFUNCTION(libhandle, FuncType, "ocp_hessi" #suffix, false);

/** Load library - implementation depends on the machine */
bool OCPInterface::load_ocplib(
    std::string name
) {
    // Load library
    LOADLIBRARY(name, libhandle);
    // Load functions
    LOADOCPFUNCS(,EvalFunc);
    LOADOCPFUNCS(_alloc_mem, AllocFunc);
    LOADOCPFUNCS(_free_mem, FreeFunc);
    LOADOCPFUNCS(_n_in, NInFunc);
    LOADOCPFUNCS(_sparsity_in, SpInFunc);
    LOADOCPFUNCS(_sparsity_out, SpOutFunc);
    LOADOCPFUNCS(_work, WorkFunc);
    // Save the library handle
    ocp_lib = (void*) libhandle;
    // Check hessian
    if (!ocp_hessb || !ocp_hessi) {
        ocp_hessb = NULL;
        ocp_hessi = NULL;
    }
    // Ok, return true
    return true;
}

bool OCPInterface::load_nlplib(
    std::string name
) {
    // Load library
    LOADLIBRARY(name, libhandle);
    // Load function
    callSolve = LOADFUNCTION(libhandle, SolveFunc, "callSolve", true);
    // Save the library handle
    nlp_lib = (void*) libhandle;
    // Ok, return true
    return true;
}

void OCPInterface::free_library() {
    if (ocp_lib) FREELIBRARY(ocp_lib);
    if (nlp_lib) FREELIBRARY(nlp_lib);
    ocp_lib = NULL;
    nlp_lib = NULL;
}
