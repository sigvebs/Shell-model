/* 
 * File:   main.cpp
 * Author: sigve
 *
 * Created on July 5, 2012, 1:33 PM
 */


using namespace std;

#include <cstdlib>
#include "defines.h"
#include "ImportSlater.h"
#include "Lanczos.h"
#include "HamiltonMatrix.h"
#include "includes/lib.h"
#include "includes/SimpleIni.h"

/*
 * 
 */


int main(int argc, char** argv) {

    // Decelerations
    string fileName;
    int M;
    int Np;
    int iterations;
    int comparisons;

    // Loading INI-file 
    SI_Error rc;
    CSimpleIniA ini;
    ini.SetUnicode();
    ini.LoadFile("param.ini");
    if (rc < 0) {
        cout << "INI-file not found." << endl;
        exit(1);
    }

    // Collecting parameters from the INI-file
    fileName = ini.GetValue("main", "filename", "default");
    M = (int) ini.GetDoubleValue("main", "M");
    Np = (int) ini.GetDoubleValue("main", "Np");
    iterations = (int) ini.GetDoubleValue("main", "iterations");
    comparisons = (int) ini.GetDoubleValue("main", "comparisons");

    cout
            << "Loading " << fileName
            << " with M = " << M * 0.5
            << " and "
            << Np << " particles"
            << " using " << BITS << " bits." << endl;

    // Importing atomic configurations from file
    ImportSlater slaterDet(M, Np, fileName, 0);
    vector<bitset<BITS> > slaterdetBin = slaterDet.getBitStates();
    mat interactions = slaterDet.getInteractionElements();
    vec oneBodyEnergy = slaterDet.getOneBodyEnergy();

    // Creating the Hamilton matrix
    HamiltonMatrix H(slaterdetBin, interactions, oneBodyEnergy);
    H.computeMatrixElements();
    
    // Running Lanczos algorithm
    Lanczos L(H.getHamiltonian(), iterations);
    L.runLanczos();
    vec eigenvalues = L.getEigenvalues();
    cout << eigenvalues.subvec(0,comparisons);
    
    exit(0);
    return 0;
}