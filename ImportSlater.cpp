/* 
 * File:   ImportSlater.cpp
 * Author: sigve
 * 
 * Created on July 5, 2012, 1:36 PM
 */

#include "ImportSlater.h"

////////////////////////////////////////////////////////////////////////////////

ImportSlater::ImportSlater() {
    cout << "ImportSlater: Not implemented! " << endl;
    exit(1);
}
////////////////////////////////////////////////////////////////////////////////

ImportSlater::ImportSlater(int M, int Np, string fileName, bool usePairs) : M(M), Np(Np), fileName(fileName), usePairs(usePairs) {
    readSlaterFromFile();
}

////////////////////////////////////////////////////////////////////////////////

ImportSlater::ImportSlater(const ImportSlater& orig) {
}

////////////////////////////////////////////////////////////////////////////////

ImportSlater::~ImportSlater() {
}

////////////////////////////////////////////////////////////////////////////////

void ImportSlater::readSlaterFromFile() {
    int rows; // Equal to number of configurations
    int counter;
    int columns = 5;

    // Example input file
    //    6
    //     1     0     2     5    -5
    //     2     0     2     5    -3
    //     3     0     2     5    -1
    //     4     0     2     5     1
    //     5     0     2     5     3
    //     6     0     2     5     5
    //   0.00000
    //   0.00000
    //   0.00000
    //   0.00000
    //   0.00000
    //   0.00000
    //     6
    //   1   6   1   6  -1.00000
    //   1   6   2   5   1.00000
    //   1   6   3   4  -1.00000
    //   2   5   2   5  -1.00000
    //   2   5   3   4   1.00000
    //   3   4   3   4  -1.00000
    string line;
    ifstream myfile(fileName.c_str());
    vector<string> v;
    string delimeter = " ";

    // Reading configurations from file, line by line and storing the results
    // in a vector.
    int numberOfStates;
    counter = 0;
    if (myfile.is_open()) {
        while (myfile.good()) {
            getline(myfile, line);

            if (counter == 0) {
                rows = atoi(line.c_str());
                localEnergies = zeros(rows, 1);
            } else if (counter <= rows) {
                // Reading configurations
                split(line, delimeter, v, columns);
            } else if (counter <= 2 * rows) {
                // Reading the local energies
                localEnergies(counter - rows - 1) = atof(line.c_str());
            } else {
                split(line, delimeter, v, columns);
            }
            counter++;
        }


        myfile.close();
    } else {
        cout << "Unable to open file";
    }

    mat configurations = zeros(rows, columns);

    counter = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            configurations(i, j) = atof(v[counter].c_str());
            counter++;
        }
    }

    // Constructing interaction matrix -----------------------------------------
    numberOfStates = atof(v[counter].c_str());
    cout << "Total number of states = " << numberOfStates << endl;

    counter++; // Skipping the number of conf.

    interactionElements = zeros(numberOfStates, columns);
    for (int i = 0; i < numberOfStates; i++) {
        for (int j = 0; j < columns; j++) {
            interactionElements(i, j) = atof(v[counter].c_str());
            counter++;
        }
    }

    // Creating all possible states --------------------------------------------
    int Nsp = rows + 1;
    vec state = zeros(Np, 1);
    int numberM = 0;
    double localM;

    bitset < BITS > binState;

    // Creating the initial state
    for (int i = 0; i < Np; i++) {
        state(i) = i + 1;
        binState.set(i);
    }

    counter = 0;
    localM = 0;
    for (int i = 0; i < Np; i++) {
        localM += configurations(state(i) - 1, 4);
    }

    if (localM == M) {
        bool pair = true;
        // Checks whether a pair is moved or not
        if (usePairs) {
            for (int i = 0; i < Np; i += 2) {
                if ((int) state(i) + 1 != (int) state(i + 1) || (int) state(i + 1) % 2 != 0) {
                    pair = false;
                    break;
                }
            }
        }
        if (pair) {
            statesBin.push_back(stateToBinary(state));
            binStates.push_back(binState);
            states.push_back(state);
            numberM++;
        }
    }

    while (true) {
        state = odometer(state, Nsp, Np);

        if (state(0) == 0)
            break;

        // Computing the total M of a state
        localM = 0;
        bool correctM = false;
        for (int i = 0; i < Np; i++) {
            localM += configurations(state(i) - 1, 4);
        }
        if (localM == M) {
            localM = 0;
            for (int i = 0; i < Np; i++) {
                localM += configurations(state(i) - 1, 4);
            }
            correctM = true;
        }

        // Checking pairs
        bool pair = true;
        // Checks whether a pair is moved or not
        if (usePairs) {
            for (int i = 0; i < Np; i += 2) {
                if ((int) state(i) + 1 != (int) state(i + 1) || (int) state(i + 1) % 2 != 0) {
                    pair = false;
                    break;
                }
            }
        }

        // If everything is OK the states are stored
        if (pair && correctM) {
            // Constructing the binary representation
            binState.reset();
            for (int i = 0; i < state.size(); i++) {
                try {
                    binState.set(state(i) - 1);
                } catch (exception& e) {
                    cout << "Exception: " << e.what() << endl;
                    cout
                            << "To run with this basis the number of bits must be increased. Current number of bits are " << BITS
                            << " , increase BITS > " << state(i) - 1 << ". Change BITS in defines.h and recompile." << endl;
                    exit(1);
                }
            }

            // Storing states.
            statesBin.push_back(stateToBinary(state));
            binStates.push_back(binState);
            states.push_back(state);
            numberM++;
        }
    }

    // Computing the local energy of a state
    double sum;
    vec st;
    localE = zeros(states.size());
    for (int i = 0; i < states.size(); i++) {
        sum = 0;
        st = states[i];
        for (int j = 0; j < st.size(); j++) {
            sum += localEnergies(st(j) - 1);
        }
        localE(i) = sum;
    }

    cout << "Number of M = " << 0.5 * M << " states is " << numberM << endl;
    cout << "Completed generating Slater Determinants." << endl;
}
////////////////////////////////////////////////////////////////////////////////

vector<vec> ImportSlater::getStates() {
    return states;
}

////////////////////////////////////////////////////////////////////////////////

vec ImportSlater::getOneBodyEnergy() {
    //    return localEnergies;
    return localE;
}

////////////////////////////////////////////////////////////////////////////////

mat ImportSlater::getInteractionElements() {
    return interactionElements;
    ;
}

////////////////////////////////////////////////////////////////////////////////

vector<int> ImportSlater::getBinStates() {
    return statesBin;
}

////////////////////////////////////////////////////////////////////////////////

vector<bitset<BITS> > ImportSlater::getBitStates() {
    return binStates;
}

////////////////////////////////////////////////////////////////////////////////

vec ImportSlater::odometer(const vec &oldState, int Nsp, int N) {
    vec newState = oldState;
    double l;

    for (int j = N - 1; j >= 0; j--) {
        if (newState(j) < Nsp - N + j) {
            l = newState(j);
            for (int k = j; k < N; k++) {
                newState(k) = l + 1 + k - j;
            }
            return newState;
        }
    }

    newState = zeros(Np, 1);
    return newState;
}

////////////////////////////////////////////////////////////////////////////////

void ImportSlater::split(const string& str, const string& delimiters, vector<string>& tokens, int columns) {
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of(delimiters, lastPos);
    int i = 0;
    while (string::npos != pos || string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
        i++;
        if (i >= columns)
            break;
    }
}

////////////////////////////////////////////////////////////////////////////////

int ImportSlater::stateToBinary(const vec & state) {
    int binRepresentation = 0;
    int nParticles = state.size();

    for (int i = 0; i < nParticles; i++) {
        binRepresentation += pow(2, state(i) - 1);
    }

    return binRepresentation;
}

////////////////////////////////////////////////////////////////////////////////