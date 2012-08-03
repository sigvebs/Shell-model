/* 
 * File:   ImportSlater.h
 * Author: sigve
 *
 * Created on July 5, 2012, 1:36 PM
 */

#ifndef IMPORTSLATER_H
#define	IMPORTSLATER_H

#include "defines.h"
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

#include <bitset>

class ImportSlater {
public:
    ImportSlater();
    ImportSlater(int, int, string, bool);
    ImportSlater(const ImportSlater& orig);
    virtual ~ImportSlater();

    void readSlaterFromFile();
    void split(const string& str, const string& delimiters, vector<string>& tokens, int);
    vec odometer(const vec &, int, int);
    int stateToBinary(const vec &);

    vector<vec> getStates();
    vector<int> getBinStates();
    vector<bitset<BITS> > getBitStates();
    vec getOneBodyEnergy();
    mat getInteractionElements();
private:
    string fileName;
    int M;
    int Np;
    vector<vec> states;
    vector<int> statesBin;
    vector<bitset<BITS> > binStates;
    bool usePairs;
    vec localEnergies;
    vec localE;
    mat interactionElements;
};

#endif	/* IMPORTSLATER_H */

