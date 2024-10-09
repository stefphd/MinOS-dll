/*  File minos_print.cpp 
    OCPInterface class implementation with out print
    Copyright (C) 2024 Stefano Lovato
*/

#define MAKE_MINOS

#include "minos.h" 
#include "macros.h" // for nlp macros
#include <stdarg.h> // for va_list

/* Forward declarations */
constexpr char OCPInterface::print_headerFormat[];
constexpr char OCPInterface::print_dataFormat[];
constexpr char OCPInterface::print_strFormat[];
constexpr char OCPInterface::print_intFormat[];
constexpr char OCPInterface::print_doubleFormat[];
constexpr char OCPInterface::print_indexFormat[];
constexpr char OCPInterface::print_rangeFormat[];
constexpr char OCPInterface::print_varnameFormat[];
constexpr char OCPInterface::print_vardataFormat[];

std::ostream& operator<<(
    std::ostream& os,
    OCPInterface& ocpInterface
) {
    // get dims
    int nx, nu, np, nc, nb, nq, nz, ng, nnzj, nnzh, na;
    int N = ocpInterface.get_N();
    ocpInterface.get_dims(&nx, &nu, &np, &nc, &nb, &nq, &nz, &ng, &nnzj, &nnzh, &na);
    // init variables
    double objval, m, inf_pr = NAN, inf_du = NAN;
    std::vector<double> t(N), x(nx*N), u(nu*N), p(np),
                        lam_x(nx*N), lam_u(nu*N), lam_p(np),
                        lam_f(nx*(N-1)), lam_c(nc*N), lam_b(nb), lam_q(nq), 
                        f(nx*(N-1)), c(nc*N), b(nb), q(nq), l(N-1);
    // call to get_sol
    ocpInterface.get_sol(&objval, t.data(),
            x.data(), u.data(), p.data(),
            lam_x.data(), lam_u.data(), lam_p.data(),
            lam_f.data(), lam_c.data(), lam_b.data(), lam_q.data(),
            f.data(), c.data(), b.data(), q.data(), 
            l.data(), &m);
    double ttot, talg, teval;
    ocpInterface.get_cpu_time(ttot, talg, teval);
    int num_iter = ocpInterface.get_num_iter();
    //double mu_curr = ocpInterface.get_mu_curr();
    if (num_iter>=0) {
        inf_pr = ocpInterface.infpr_history[num_iter];
        inf_du = ocpInterface.infdu_history[num_iter];
    }
    // Print to os
    os << OCPInterface::get_header();
    os << "\n";

    os << "# Problem settings\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_strFormat, "name:", ocpInterface.name.c_str()) << "\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_doubleFormat, "ti:", ocpInterface.ti) << "\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_doubleFormat, "tf:", ocpInterface.tf) << "\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_intFormat, "N:", N) << "\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_strFormat, "solver:", ocpInterface.nlpsolver.c_str()) << "\n";
    os << "\n";

    os << "# Problem dimensions\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_intFormat, "nx:", nx) << "\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_intFormat, "nu:", nu) << "\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_intFormat, "np:", np) << "\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_intFormat, "nc:", nc) << "\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_intFormat, "nb:", nb) << "\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_intFormat, "nq:", nq) << "\n";
    os << "\n";

    os << "# Problem bounds\n";
    std::string varname;
    for (int i = 0; i < nx; ++i) {
        varname = "x[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]:";
        os << "\t" << OCPInterface::format_str(OCPInterface::print_rangeFormat, varname.c_str(), ocpInterface.lbx[i], ocpInterface.ubx[i]) << "\n";
    }
    for (int i = 0; i < nu; ++i) {
        varname = "u[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]:";
        os << "\t" << OCPInterface::format_str(OCPInterface::print_rangeFormat, varname.c_str(), ocpInterface.lbu[i], ocpInterface.ubu[i]) << "\n";
    }
    for (int i = 0; i < np; ++i) {
        varname = "p[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]:";
        os << "\t" << OCPInterface::format_str(OCPInterface::print_rangeFormat, varname.c_str(), ocpInterface.lbp[i], ocpInterface.ubp[i]) << "\n";
    }
    for (int i = 0; i < nc; ++i) {
        varname = "c[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]:";
        os << "\t" << OCPInterface::format_str(OCPInterface::print_rangeFormat, varname.c_str(), ocpInterface.lbc[i], ocpInterface.ubc[i]) << "\n";
    }
    for (int i = 0; i < nb; ++i) {
        varname = "b[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]:";
        os << "\t" << OCPInterface::format_str(OCPInterface::print_rangeFormat, varname.c_str(), ocpInterface.lbb[i], ocpInterface.ubb[i]) << "\n";
    }
    for (int i = 0; i < nq; ++i) {
        varname = "q[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]:";
        os << "\t" << OCPInterface::format_str(OCPInterface::print_rangeFormat, varname.c_str(), ocpInterface.lbq[i], ocpInterface.ubq[i]) << "\n";
    }
    os << "\n";

    if (na > 0) {
        os << "# Auxdata\n";
        for (int i = 0; i < na; ++i)  {
            varname = "auxdata[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]:";
            os << "\t" << OCPInterface::format_str(OCPInterface::print_doubleFormat, varname.c_str(), ocpInterface.auxdata[i]) << "\n";
        }
        os << "\n";
    }

    os << "# Statistics\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_intFormat, "NLP nz:", nz) << "\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_intFormat, "NLP ng:", ng) << "\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_intFormat, "NLP jacobian nnz:", nnzj) << "\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_intFormat, "NLP hessian nnz:", nnzh) << "\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_intFormat, "Number of iters:", num_iter) << "\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_doubleFormat, "Final inf pr:", inf_pr)  << "\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_doubleFormat, "Final inf du:", inf_du) << "\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_doubleFormat, "CPU total time (s):", ttot) << "\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_doubleFormat, "CPU algorithm time (s):", talg) << "\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_doubleFormat, "CPU evaluation time (s):", teval) << "\n";
    os << "\n";

    os << "# Solution\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_doubleFormat, "obj:", objval) << "\n";
    os << "\t" << OCPInterface::format_str(OCPInterface::print_doubleFormat, "m:", m) << "\n";
    for (int i = 0; i < np; ++i) {
        varname = "p[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]:";
        os << "\t" << OCPInterface::format_str(OCPInterface::print_doubleFormat, varname.c_str(), p[i]) << "\n";
    }
    for (int i = 0; i < np; ++i) {
        varname = "lam_p[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]:";
        os << "\t" << OCPInterface::format_str(OCPInterface::print_doubleFormat, varname.c_str(), lam_p[i]) << "\n";
    }
    for (int i = 0; i < nb; ++i) {
        varname = "b[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]:";
        os << "\t" << OCPInterface::format_str(OCPInterface::print_doubleFormat, varname.c_str(), b[i]) << "\n";
    }
    for (int i = 0; i < nb; ++i) {
        varname = "lam_b[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]:";
        os << "\t" << OCPInterface::format_str(OCPInterface::print_doubleFormat, varname.c_str(), lam_b[i]) << "\n";
    }
    for (int i = 0; i < nq; ++i) {
        varname = "q[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]:";
        os << "\t" << OCPInterface::format_str(OCPInterface::print_doubleFormat, varname.c_str(), q[i]) << "\n";
    }
    for (int i = 0; i < nq; ++i) {
        varname = "lam_q[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]:";
        os << "\t" << OCPInterface::format_str(OCPInterface::print_doubleFormat, varname.c_str(), lam_q[i]) << "\n";
    }
    os << "\n";
    // create header and data format
    std::ostringstream header;
    std::string dataFormat;
    header << OCPInterface::format_str(OCPInterface::print_varnameFormat, "t") << "\t";
    for (int i = 0; i < nx; ++i) {
        varname = "x[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]";
        header << OCPInterface::format_str(OCPInterface::print_varnameFormat, varname.c_str()) << "\t";
    }
    for (int i = 0; i < nu; ++i) {
        varname = "u[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]";
        header << OCPInterface::format_str(OCPInterface::print_varnameFormat, varname.c_str()) << "\t";
    }
    for (int i = 0; i < nx; ++i) {
        varname = "lam_x[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]";
        header << OCPInterface::format_str(OCPInterface::print_varnameFormat, varname.c_str()) << "\t";
    }
    for (int i = 0; i < nu; ++i) {
        varname = "lam_u[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]";
        header << OCPInterface::format_str(OCPInterface::print_varnameFormat, varname.c_str()) << "\t";
    }
    for (int i = 0; i < nx; ++i) {
        varname = "lam_f[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]";
        header << OCPInterface::format_str(OCPInterface::print_varnameFormat, varname.c_str()) << "\t";
    }
    for (int i = 0; i < nc; ++i) {
        varname = "lam_c[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]";
        header << OCPInterface::format_str(OCPInterface::print_varnameFormat, varname.c_str()) << "\t";
    }
    for (int i = 0; i < nx; ++i) {
        varname = "f[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]";
        header << OCPInterface::format_str(OCPInterface::print_varnameFormat, varname.c_str()) << "\t";
    }
    for (int i = 0; i < nc; ++i) {
        varname = "c[" + OCPInterface::format_str(OCPInterface::print_indexFormat, i) + "]";
        header << OCPInterface::format_str(OCPInterface::print_varnameFormat, varname.c_str()) << "\t";
    }
    header << OCPInterface::format_str(OCPInterface::print_varnameFormat, "l");
    os << header.str() << "\n";
    
    for (int k = 0; k < N; ++k) {
        os << OCPInterface::format_str(OCPInterface::print_vardataFormat, t[k]) << "\t";
        for (int i = 0; i < nx; ++i) os << OCPInterface::format_str(OCPInterface::print_vardataFormat, x[nx*k+i]) << "\t";
        for (int i = 0; i < nu; ++i) os << OCPInterface::format_str(OCPInterface::print_vardataFormat, u[nu*k+i]) << "\t";
        for (int i = 0; i < nx; ++i) os << OCPInterface::format_str(OCPInterface::print_vardataFormat, lam_x[nx*k+i]) << "\t";
        for (int i = 0; i < nu; ++i) os << OCPInterface::format_str(OCPInterface::print_vardataFormat, lam_u[nu*k+i]) << "\t";
        for (int i = 0; i < nx; ++i) {
            if ( k == (N-1) ) os << OCPInterface::format_str(OCPInterface::print_vardataFormat, NAN) << "\t";
            else os << OCPInterface::format_str(OCPInterface::print_vardataFormat, lam_f[nx*k+i]) << "\t";
        }
        for (int i = 0; i < nc; ++i) os << OCPInterface::format_str(OCPInterface::print_vardataFormat, lam_c[nc*k+i]) << "\t";
        for (int i = 0; i < nx; ++i) {
            if ( k == (N-1) ) os << OCPInterface::format_str(OCPInterface::print_vardataFormat, NAN) << "\t";
            else os << OCPInterface::format_str(OCPInterface::print_vardataFormat, f[nx*k+i]) << "\t";
        }
        for (int i = 0; i < nc; ++i) os << c[nc*k+i] << "\t";
        if ( k == (N-1) ) os << OCPInterface::format_str(OCPInterface::print_vardataFormat, NAN) << "\t";
        else os << OCPInterface::format_str(OCPInterface::print_vardataFormat, l[k]) << "\n";
    }
    os << "\n";
    return os;
}


void OCPInterface::print_sol(
    std::ostream &os
) {
    os << *this;
}

std::string OCPInterface::toString() {
    std::ostringstream oss;
    oss << *this;
    return oss.str();
}

std::string OCPInterface::get_header() {
    int numLines, maxLength, lineLegth;
    std::ostringstream str;
    std::vector<std::string> lines;
    std::vector<int> lengths;
    // Set header lines
    lines.push_back(PROJECT_DESCRIPTION);
    lines.push_back(OCPInterface::format_str("Version %s", get_version()));
    lines.push_back(OCPInterface::format_str("Copyright (C) %s %s <%s>", COPYRIGHT_YEAR, COPYRIGHT_NAME, COPYRIGHT_CONTACT));
    //lines.push_back(OCPInterface::format_str("Copyright (C) %s %s", COPYRIGHT_YEAR, COPYRIGHT_NAME));
    //lines.push_back(OCPInterface::format_str("<%s>", COPYRIGHT_CONTACT));
    numLines = lines.size();
    // Resize lengths vector to match lines and populate lengths
    lengths.resize(numLines);
    for (int i = 0; i < numLines; ++i)
        lengths[i] = lines[i].length();
    // Find max length
    maxLength = 0;
    for (int i = 0; i < numLines; ++i)
        maxLength = (lengths[i] > maxLength) ? lengths[i] : maxLength;
    // Add white spaces to strings
    maxLength += 1*2; // add white spaces at the line borders
    for (int i = 0; i < numLines; ++i) {
        int numLeftSpaces  = (maxLength - lines[i].length()) / 2;
        int numRightSpaces = (maxLength - lines[i].length()) - numLeftSpaces; // may differ from numLeftSpaces due to odd length
        // Prepend left spaces
        lines[i].insert(0, numLeftSpaces, ' '); 
        // Append right spaces
        lines[i].append(numRightSpaces, ' ');
        // Prepend and append '*'
        lines[i].insert(0, 1, '*'); 
        lines[i].append(1, '*');
    }
    maxLength += 2; // include '*' at the beggining and the end
    // Add first and end line of '*'
    std::string starline = "";
    starline.append(maxLength, '*');
    lines.insert(lines.begin(), starline);
    lines.push_back(starline);
    // Print to ostringsteam
    numLines = lines.size();
    for (int i = 0; i < numLines; ++i)
        str << lines[i] << std::endl;
    return str.str();
}

const char* OCPInterface::get_version() {
    return PROJECT_VERSION;
}

/** Print statistics */
void OCPInterface::print_stats(
    const int status,
    const int succed,
    const int acceptable,
    const int userstop,
    const int maxiter,
    const int infeasible
) {
    if (status == succed)
        (*print_funptr)("EXIT: Optimal Solution Found.\n");
    else if (status == acceptable)
        (*print_funptr)("EXIT: Solved To Acceptable Level.\n");
    else if (status == userstop) 
        (*print_funptr)("EXIT: Stopping optimization at current point as requested by user.\n");
    else if (status == maxiter) 
        (*print_funptr)("EXIT: Maximum Number of Iterations Exceeded.\n");
    else if (status == infeasible) 
        (*print_funptr)("EXIT: Converged to a point of local infeasibility. Problem may be infeasible.\n");
    else
        (*print_funptr)("EXIT: Unable To Solve (Status: %d).\n", status);

    (*print_funptr)("Elapsed time is %.6f seconds.\n", (tcpu_tot   != NAN) ? tcpu_tot : 0);
    (*print_funptr)("Algorithm time is %.6f seconds.\n", (tcpu_eval  != NAN) ? (tcpu_tot-tcpu_eval) : tcpu_tot);
    (*print_funptr)("Function evaluation time is %.6f seconds.\n", (tcpu_eval != NAN) ? tcpu_eval : 0);
    (*print_funptr)("Number of Iterations is %d\n", num_iter);
}

/** Print iteration */
bool OCPInterface::print_iter(
    const int iter,
    const double obj,
    const double inf_pr,
    const double inf_du,
    const bool force
) {
    // Record history
    this->num_iter = iter;
    this->J_opt = obj;
    this->obj_history[iter] = obj;
    this->infpr_history[iter] = inf_pr;
    this->infdu_history[iter] = inf_du;

    // Check interrupts, call function pointed by int_funptr 
    if ((*int_funptr)() && !force) {
        return false; // stop NLP solver
    }

    // Check for not display
    if (!this->display) { 
        return true; 
    }

    // Print sol to file every print_itersol iterations or if forced
    if ((print_itersol>0) && (((iter % print_itersol) == 0) ^ force)) {
        /* Print solution to file */
        std::ofstream outfile;
        char strbuf[256];
        sprintf(strbuf, "sol-iter%d.txt", iter);
        outfile.open(strbuf);
        outfile << *this;
        outfile.close();
    }
    
    // Print header every HEADERPRINTINT iters, do not print if forced
    if ((( iter % HEADERPRINTINT) == 0) && !force)
        (*print_funptr)(print_headerFormat, "iter", "obj", "inf_pr", "inf_du");
    
    // Print data every ITERPRINTINT or if forced
    if ((( iter % ITERPRINTINT) == 0) ^ force)
        (*print_funptr)(print_dataFormat, iter, obj, inf_pr, inf_du);
    
    return true;
}

/** Implementation of format strings */
std::string OCPInterface::format_str(
    const char *fmt, 
    ...
) {
    // handle variable args
    va_list ap, apcp;
    va_start(ap, fmt);
    va_start(ap, fmt); va_copy(apcp, ap); 
    int len = vsnprintf(NULL, 0, fmt, apcp);
    va_end(apcp);
    if (len == -1) { return ""; }
    char *ret = (char*) malloc((size_t) len + 1);
    if (!ret) { return ""; }
    len = vsprintf(ret, fmt, ap);
    va_end(ap);
    if (len == -1) {
        free(ret);
        return "";
    }
    std::string str(ret);
    free(ret);
    return str;
}