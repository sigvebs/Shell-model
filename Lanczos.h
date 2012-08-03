/* 
 * File:   Lanczos.h
 * Author: sigve
 *
 * Created on July 6, 2012, 2:35 PM
 */

#ifndef LANCZOS_H
#define	LANCZOS_H

#include<iostream>
#include<fstream>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

class Lanczos {
public:
    Lanczos();
    Lanczos(mat, int);
    Lanczos(const Lanczos& orig);
    virtual ~Lanczos();
    void runLanczos();

    vec gramSchmith(const mat &, int);
    void compareSolution(int);
    vec getEigenvalues();
private:
    mat H;
    int iterations;
    vec eigvalLanczos;
};

#endif	/* LANCZOS_H */

