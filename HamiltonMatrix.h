/* 
 * File:   HamiltonMatrix.h
 * Author: sigve
 *
 * Created on July 9, 2012, 12:00 PM
 */

#ifndef HAMILTONMATRIX_H
#define	HAMILTONMATRIX_H

#include<iostream>
#include<fstream>
#include <vector>
#include <armadillo>
using namespace std;
using namespace arma;

#include "defines.h"
#include <bitset>

class HamiltonMatrix {
public:
    HamiltonMatrix();
    HamiltonMatrix(vector<int>, mat, vec);
    HamiltonMatrix(vector<bitset<BITS> >, mat, vec);
    HamiltonMatrix(const HamiltonMatrix& orig);
    virtual ~HamiltonMatrix();

    void computeMatrixElements();
//    int addParticle(int, int);
//    int removeParticle(int, int);
    bitset<BITS> addParticle(int, bitset<BITS> );
    bitset<BITS> removeParticle(int, bitset<BITS> );
    void int2bin(int);
//    int sign(int, int);
    int sign(int, bitset<BITS> state);
    mat getHamiltonian();
    vec getSingleEnergy();
private:
    mat H;
    vector<int> SLBin;
    vector<bitset<BITS> > SLBit;
    mat interactions;
    vec oneBodyEnergy;
};

#endif	/* HAMILTONMATRIX_H */

