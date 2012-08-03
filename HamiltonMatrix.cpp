/* 
 * File:   HamiltonMatrix.cpp
 * Author: sigve
 * 
 * Created on July 9, 2012, 12:00 PM
 */

#include "HamiltonMatrix.h"

////////////////////////////////////////////////////////////////////////////////

HamiltonMatrix::HamiltonMatrix() {
}

////////////////////////////////////////////////////////////////////////////////

HamiltonMatrix::HamiltonMatrix(vector<int> SLBin, mat interactions, vec oneBodyEnergy) : SLBin(SLBin), interactions(interactions), oneBodyEnergy(oneBodyEnergy) {
}

////////////////////////////////////////////////////////////////////////////////

HamiltonMatrix::HamiltonMatrix(vector<bitset<BITS> > SLBit, mat interactions, vec oneBodyEnergy) : SLBit(SLBit), interactions(interactions), oneBodyEnergy(oneBodyEnergy) {
}

////////////////////////////////////////////////////////////////////////////////

HamiltonMatrix::HamiltonMatrix(const HamiltonMatrix& orig) {
}

////////////////////////////////////////////////////////////////////////////////

HamiltonMatrix::~HamiltonMatrix() {
}

////////////////////////////////////////////////////////////////////////////////

void HamiltonMatrix::computeMatrixElements() {
    cout << "Computing the Hamilton matrix " << endl;
    int s;
    int nStates = SLBit.size();

    //    int newState;
    bitset<BITS> newState;
    int i, j, k, l;
    double interactionElement;    
    H = zeros(nStates, nStates); 
    
    for (int st = 0; st < SLBit.size(); st++) {
        for (int b = 0; b < interactions.n_rows; b++) {
            newState = SLBit[st];
            i = interactions(b, 0);
            j = interactions(b, 1);
            k = interactions(b, 2);
            l = interactions(b, 3);
            interactionElement = interactions(b, 4);

            // Using the two body operator. Mathematics: a+(i)a+(j)a-(l)a-(k)|psi>
            s = 1;
            newState = removeParticle(k, newState);
            s *= sign(k, newState);

            newState = removeParticle(l, newState);
            s *= sign(l, newState);

            newState = addParticle(j, newState);
            s *= sign(j, newState);

            newState = addParticle(i, newState);
            s *= sign(i, newState);

            if (newState != 0) {
                // Searching for the new state to compute the matrix elements.
                for (int st2 = 0; st2 < SLBit.size(); st2++) {
                    if (newState == SLBit[st2]) {
                        H(st, st2) += interactionElement * s;
                    }
                }
            }
        }

        // One body energies
        H(st, st) += oneBodyEnergy(st);
    }

    // Symmetrizing the Hamilton matrix
    for (int i = 0; i < H.n_rows; i++) {
        for (int j = 0; j < i; j++) {
            H(j, i) = H(i, j);
        }
    }
    cout << "Completed generating the Hamilton matrix" << endl;
}

////////////////////////////////////////////////////////////////////////////////

int HamiltonMatrix::sign(int n, bitset<BITS> state) {
    int s = 1;
    for (int i = 0; i < n - 1; i++) {
        if (state[i] != 0) {
            s *= -1;
        }
    }
    return s;
}

////////////////////////////////////////////////////////////////////////////////

bitset<BITS> HamiltonMatrix::addParticle(int n, bitset<BITS> state) {
    // Vacuum state
    if (state.any() == false)
        return state;

    bitset<BITS> a;
    a.set(n - 1);
    bitset<BITS> comp = a & state;

    if (comp.count() == false) {
        state.set(n - 1);
    } else {
        state.reset();
    }

    return state;
}

////////////////////////////////////////////////////////////////////////////////

bitset<BITS> HamiltonMatrix::removeParticle(int n, bitset<BITS> state) {
    // Vacuum state
    if (state.any() == false)
        return state;

    bitset<BITS> a;
    a.set(n - 1);
    bitset<BITS> comp = a & state;

    if (comp.count() == true) {
        state.set(n - 1, 0);
    } else {
        state.reset();
    }

    return state;
}
////////////////////////////////////////////////////////////////////////////////

mat HamiltonMatrix::getHamiltonian() {
    return H;
}

////////////////////////////////////////////////////////////////////////////////

void HamiltonMatrix::int2bin(int a) {
    int buf_size = 33;
    char *buffer = new char[buf_size];
    buffer[buf_size - 1] = '\0';
    buffer += (buf_size - 1);

    for (int i = 31; i >= 0; i--) {
        *(--buffer) = (a & 1) + '0';
        a >>= 1;
    }

    cout << buffer << endl;
    delete buffer;
}

////////////////////////////////////////////////////////////////////////////////


/*
////////////////////////////////////////////////////////////////////////////////

void HamiltonMatrix::computeMatrixElements() {
    int n;
    int s;
    int nStates = SLBin.size();

    int newState;
    int i, j, k, l;
    double interactionElement;
    H = zeros(nStates, nStates);
    for (int st = 0; st < SLBin.size(); st++) {
        for (int b = 0; b < interactions.n_rows; b++) {
            newState = SLBin[st];
            i = interactions(b, 0);
            j = interactions(b, 1);
            k = interactions(b, 2);
            l = interactions(b, 3);
            interactionElement = interactions(b, 4);

            // Using the two body operator. Mathematics: a+(i)a+(j)a-(l)a-(k)|psi>
            s = 1;
            newState = removeParticle(k, newState);
            s *= sign(k, newState);

            newState = removeParticle(l, newState);
            s *= sign(l, newState);

            newState = addParticle(j, newState);
            s *= sign(j, newState);

            newState = addParticle(i, newState);
            s *= sign(i, newState);

            if (newState != 0) {
                // Searching for the new state to compute the matrix elements.
                for (int st2 = 0; st2 < SLBin.size(); st2++) {
                    if (newState == SLBin[st2]) {
                        H(st, st2) += interactionElement * s;
                    }
                }
            }
        }

        // One body energies
        H(st, st) += oneBodyEnergy(st);
    }

    // Symmetrizing the Hamilton matrix
    for (int i = 0; i < H.n_rows; i++) {
        for (int j = 0; j < i; j++) {
            H(j, i) = H(i, j);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

int HamiltonMatrix::sign(int n, int state) {
    int s = 1;
    int set;
    for (int i = 0; i < n - 1; i++) {
        set = state & (1 << i);
        if (set != 0) {
            s *= -1;
        }
    }
    return s;
}

////////////////////////////////////////////////////////////////////////////////

int HamiltonMatrix::addParticle(int n, int state) {
    if (state == 0)
        return state;
    int a = pow(2, n - 1);
    int comp = a & state;
    int newState;

    if (comp == 0) {
        newState = a + state;
    } else {
        newState = 0;
    }

    return newState;
}

////////////////////////////////////////////////////////////////////////////////

int HamiltonMatrix::removeParticle(int n, int state) {
    if (state == 0)
        return state;
    int a = pow(2, n - 1);
    int comp = a & state;
    int newState;

    if (comp != 0) {
        newState = state - a;
    } else {
        newState = 0;
    }

    return newState;
}
////////////////////////////////////////////////////////////////////////////////
 */