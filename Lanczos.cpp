/* 
 * File:   Lanczos.cpp
 * Author: sigve
 * 
 * Created on July 6, 2012, 2:35 PM
 */

#include "Lanczos.h"
#include "includes/lib.h"

////////////////////////////////////////////////////////////////////////////////

Lanczos::Lanczos() {

    // Generating Hamilton Matrix
    double G = 1;
    double delta = 0;
    H << 2 * delta - G << -G << -G << -G << -G << 0 << endr
            << -G << 4 * delta - G << -G << -G << -G << -G << endr
            << -G << -G << 6 * delta - G << 0 << -G << -G << endr
            << -G << -G << 0 << 6 * delta - G << -G << -G << endr
            << -G << -G << -G << -G << 8 * delta - G << -G << endr
            << 0 << -G << -G << -G << -G << 10 * delta - G << endr;
    int d = 1000;

    H = randu(d, d);
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < i; j++) {
            H(i, j) = H(j, i);
        }
    }
    iterations = 100;
}

////////////////////////////////////////////////////////////////////////////////

Lanczos::Lanczos(mat H, int iterations) : H(H), iterations(iterations) {
    if (iterations > H.n_rows) {
        cout << "Number of iterations cannot be larger then the matrix." << endl;
        cout << "Setting iterations equal to the size of the matrix." << endl;
        this->iterations = H.n_rows;
    }
}

////////////////////////////////////////////////////////////////////////////////

Lanczos::Lanczos(const Lanczos& orig) {
}

////////////////////////////////////////////////////////////////////////////////

Lanczos::~Lanczos() {
}

////////////////////////////////////////////////////////////////////////////////

void Lanczos::runLanczos() {
    cout << "Starting Lanczos algorithm" << endl;

    // Running the Lanczos algorithm
    int n = H.n_rows;
    mat I = zeros(n, n);
    I.eye();

    vec alpha = zeros(n + 1, 1);
    vec beta = zeros(n + 1, 1);

    mat q = zeros(n, n + 1);
    vec r = 2 * randu(n, 1);
    r = 1.0 / norm(r, 2) * r;
    q.col(1) = r;
    int k = 0;
    beta(0) = 1;

    while (k < iterations) {
        q.col(k + 1) = r / beta(k);
        q.col(k + 1) = gramSchmith(q, k);
        k++;
        alpha(k) = dot(q.col(k), H * q.col(k));
        r = (H - alpha(k) * I) * q.col(k) - beta(k - 1) * q.col(k - 1);
        beta(k) = norm(r, 2);
    }

    //    I = zeros(n, n);
    //
    //    for (int i = 0; i < n - 1; i++) {
    //        I(i, i) = alpha(i + 1);
    //        I(i + 1, i) = I(i, i + 1) = beta(i + 1);
    //    }
    //
    //    eigvalLanczos = eig_sym(I);

    // Using TQLI
    cout << "Diagonalize the tridiagonal matrix" << endl;
    double *d = new double[iterations];
    double *b = new double[iterations];
    for (int i = 0; i < iterations; i++) {
        d[i] = alpha(i + 1);
        b[i] = beta(i);
    }

    double **z = (double**) calloc(iterations, sizeof (double*));

    for (int i = 0; i < iterations; i++)
        z[i] = (double*) calloc(iterations, sizeof (double));

    tqli(d, b, iterations, z);

    vec diag = zeros(iterations);

    for (int i = 0; i < iterations; i++)
        diag(i) = d[i];

    eigvalLanczos = sort(diag, 0);

    // Freeing memory
    for (int i = 0; i < iterations; i++)
        delete z[i];
    delete d, b;
}

////////////////////////////////////////////////////////////////////////////////

vec Lanczos::gramSchmith(const mat &q, int i) {
    vec newQ = q.col(i + 1);
    newQ -= q.col(i) * dot(q.col(i), q.col(i + 1));

    for (int p = 0; p < i; p++) {
        newQ -= q.col(p) * dot(q.col(p), q.col(i + 1));
    }
    newQ /= norm(newQ, 2);

    return newQ;
}

////////////////////////////////////////////////////////////////////////////////

void Lanczos::compareSolution(int n) {
    cout << "Starting comparison" << endl;
    // Armadillos function for finding the eigenvalues
    vec eigval = eig_sym(H);
    n--;

    cout << eigvalLanczos.subvec(0, n) << endl;
    cout << eigval.subvec(0, n) << endl;
    cout << "Difference = " << endl << eigval.subvec(0, n) - eigvalLanczos.subvec(0, n) << endl;
    //    n--;
    //
    //    int end;
    //    if (n > iterations)
    //        end = eigval.n_elem;
    //    else
    //        end = eigval.n_elem - iterations + n;
    //
    //    cout << eigvalLanczos.subvec(eigvalLanczos.n_elem - iterations, end) << endl;
    //    cout << eigval.subvec(0, n) << endl;
    //    cout << "Difference = " << endl << eigval.subvec(0, n) - eigvalLanczos.subvec(eigvalLanczos.n_elem - iterations, end) << endl;
}

////////////////////////////////////////////////////////////////////////////////

vec Lanczos::getEigenvalues() {
    return eigvalLanczos;
}

////////////////////////////////////////////////////////////////////////////////

//#if 0 // Bills algorithm

//    vec qNew = zeros(n, 1);
//    vec qOld = zeros(n, 1);
//    vec w;
//    vec v_prev;
//    vec v;
//    vec b = randu(n, 1);
//    v = 1.0 / norm(b, 2) * b;
//    v_prev = zeros(n, 1);
//    beta(0) = 1;
//    for (int k = 1; k <= n; k++) {
//        w = H*v;
//        alpha(k) = dot(v, w);
//        w = w - beta(k - 1) * v_prev - alpha(k) * v;
//        beta(k) = norm(w, 2);
//        v_prev = v;
//        v = 1.0 / beta(k) * w;
//    }
//
//    // Constructing the tridiagonal matrix
//    mat T = zeros(n, n);
//    for (int i = 0; i < n - 1; i++) {
//        T(i, i) = alpha(i + 1);
//        T(i + 1, i) = T(i, i + 1) = beta(i + 1);
//    }
//
//    T(n - 1, n - 1) = alpha(n);
//    vec eigvalLanczos = eig_sym(T);
//    //    cout << "T = " << endl << T << endl;
//    cout << "Lanczos " << endl;
//    cout << eigvalLanczos << endl;
//#else
//    // Morten's algorithm
//    mat q = zeros(n,n);
//    vec r = randu(n, 1);
//    r = 1.0 / norm(r, 2) * r;
//    qOld = zeros(n, 1);
//    int k = 0;
//    beta(0) = 1;
//    while (k < n) {
//        qNew = r / beta(k);
//        k++;
//        alpha(k) = dot(qNew, H * qNew);
//        r = (H - alpha(k) * I) * qNew - beta(k - 1) * qOld;
//        qOld = qNew;
//        beta(k) = norm(r, 2);
//    }
//    
//    for (int i = 0; i < n - 1; i++) {
//        T(i, i) = alpha(i + 1);
//        T(i + 1, i) = T(i, i + 1) = beta(i + 1);
//    }
//    
//    T(n - 1, n - 1) = alpha(n);
//    eigvalLanczos = eig_sym(T);
//    cout << "Mortens Lanczos " << endl;
//    cout << eigvalLanczos << endl;

//#endif