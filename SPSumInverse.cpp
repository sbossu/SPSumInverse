// SPSumInverse.cpp : This file contains the 'main' console function

#include <iostream>
#include <string>
#include <format>
#include <random>

//x=0 z=1 tol=1e-6 -x 2>&1 <$(IntDir)abc.txt

//Declare low-level C algorithm
extern "C" int c_SPSumInverse(int n, double tol, double oa[], double xb[], double oc[],
    double Q[], double Pi[], double zv[], double obeta[], bool verbose = true);

//Declare utility functions
double str2d(std::string s);
int findany(std::string s, const char range[], unsigned int offset = 0);
int readarray(double*& x, unsigned int offset);

#define ABORT(errcode,msg,...){ fprintf(stderr,"Error %d: " #msg "\n",errcode,__VA_ARGS__); return errcode; }

/***** MAIN FUNCTION *****/

/// <summary>Syntax: SPSumInverse [x=...] [z=...] [tol=...] [-q,--quiet]</summary>
/// <returns>>=0 if successul, <0 otherwise</returns>
int main(int argc, char* argv[]) {
    double tol = 1e-9, x = 0, z = 1;
    int i, j, n, result = 0;
    bool verbose = true;
    double* oa = nullptr, * xb = nullptr, * oc = nullptr, * Q, * Pi, * zv, * obeta;

    //Read arguments
    for (i = 1; i < argc; i++)
        switch (argv[i][0]) {
        case 'x':
            if (argv[i][1] = '=')
                if (isfinite(x = str2d(&argv[i][2])))
                    break;
            goto badarg;
        case 'z':
            if (argv[i][1] = '=')
                if (isfinite(z = str2d(&argv[i][2])) && z != 0)
                    break;
            goto badarg;
        case 't':
            if (!strncmp(argv[i], "tol=", 4))
                if (isfinite(tol = str2d(&argv[i][4])) && tol > 0)
                    break;
            goto badarg;
        case '-':
            if (argv[i][1] == 'q' || !strcmp(argv[i], "--quiet")) {
                verbose = false;
                break;
            }
            goto badarg;
        default:
        badarg:
            ABORT(-257, "Invalid argument '%s'. Expected [x=...] [z=...] [tol=...] [-q,--quiet]", argv[i]);
            break;
        }

    //Read data from standard input
    if (verbose)
        std::cout << "Enter on each line the values of a[], b[], c[] separated by a single space, tab, comma or semicolon:\n";
    n = readarray(oa, 1) - 1;
    if (n <= 0) ABORT(-273, "Invalid a[] input");
    if (readarray(xb, 1) != n + 1) ABORT(-274, "Invalid b[] input");
    xb[0] = x;
    if (readarray(oc, 1) != n + 1) ABORT(-275, "Invalid c[] input");

    //Initialize output arrays
    Q = new double[n * (n + 1) / 2];
    Pi = new double[n * (n + 1) / 2 - 1];
    zv = new double[n + 1];
    zv[0] = z;
    obeta = new double[n + 1];

    //Call low-level C algorithm and display to standard output
    result = c_SPSumInverse(n, tol, oa, xb, oc, Q, Pi, zv, obeta, verbose);
    if (result >= 0) {
        if (verbose) std::cout << "Inverse matrix:\n";
        std::cout.precision(std::numeric_limits<double>::max_digits10 - 1);
        double qvalue;
        for (i = 1; i <= n; i++)
            for (j = 1; j <= n; j++) {
                qvalue = (i <= j ? Q[(2 * n - i) * (i - 1) / 2 + j - 1] : Q[(2 * n - j) * (j - 1) / 2 + i - 1]);
                if (!isfinite(qvalue)) ABORT(-289, "Invalid output Q[%d,%d] = %f", i, j, qvalue);
                std::cout << std::scientific << qvalue;
                if (j == n) std::cout << std::endl; else std::cout << ',';
            }
    }
    return result;
}

/***** UTILITY FUNCTIONS *****/

/// <summary>Convert string to double</summary>
/// <returns>Value if successful, nan otherwise</returns>
double str2d(std::string s) {
    if (s.length() == 0)
        return NAN;
    char* end;
    errno = 0;
    double res = std::strtod(s.c_str(), &end);
    if (s.c_str() != end && errno == 0)
        return res;
    else
        return NAN;
}
/// <summary>Find position within string of any separator in specified range</summary>
/// <param name="s">String to search</param>
/// <param name="range">Null-terminated array of separators</param>
/// <param name="offset">Position to start search</param>
/// <returns>Position if successful, -1 otherwise</returns>
int findany(std::string s, const char range[], unsigned int offset) {
    static int i, j;
    static char c;
    for (i = offset; i < s.length(); i++)
        for (c = range[0], j = 0; c; c = range[++j])
            if (s[i] == c) return i;
    return std::string::npos;
}
/// <summary>Read array of doubles from standard input</summary>
/// <param name="x">Output array of doubles</param>
/// <param name="offset">Position of first input in array</param>
/// <returns>Array size >0 if successful, <=0 otherwise</returns>
int readarray(double*& x, unsigned int offset) {
    const char separators[] = ",; \t";
    std::string line;
    std::getline(std::cin, line);
    if (!line.length()) return 0;
    int n = offset, pos = -1;
    do { //Count number of separators
        n++;
        pos = findany(line, separators, pos + 1);
    } while (pos != std::string::npos);

    //Allocate and read output array
    x = new double[n];
    int start = 0, end = 0;
    for (int i = offset; i < n; i++) {
        end = findany(line, separators, start);
        if (end < start) end = line.length();
        x[i] = str2d(line.substr(start, end - start + 1));
        if (!isfinite(x[i])) {
            delete[] x;
            x = nullptr;
            return -i - 1;
        }
        start = end + 1;
    }
    return n;
}