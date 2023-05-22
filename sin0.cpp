#include "solution.h"

int main() {
    vector<double> sin_inits0(SPACE_SCALE_N_0); // sin initial conditions
    vector<double> sin_inits1(SPACE_SCALE_N_0);
    vector<double> sin_inits2(SPACE_SCALE_N_0);

    cout << check_stability(3, TIME_SPLIT_T_0, SPACE_SPLIT_H_0) << endl;
    cout << check_stability(3, TIME_SPLIT_T_1, SPACE_SPLIT_H_1) << endl;
    cout << check_stability(3, TIME_SPLIT_T_2, SPACE_SPLIT_H_2) << endl;

    for (int i = 0; i < SPACE_SCALE_N_0; i++) { // init sin with h = SPACE_SPLIT_0
        sin_inits0.at(i) = sin(M_PI * SPACE_SPLIT_H_0 * i);
    }

    for (int i = 0; i < SPACE_SCALE_N_0; i++) { // init sin with h = SPACE_SPLIT_1
        sin_inits1.at(i) = sin(M_PI * SPACE_SPLIT_H_1 * i);
    }

    for (int i = 0; i < SPACE_SCALE_N_0; i++) { // init sin with h = SPACE_SPLIT_2
        sin_inits2.at(i) = sin(M_PI * SPACE_SPLIT_H_2 * i);
    }

    vector<double> rows0(SPACE_SCALE_N_0, 0);

    vector<vector<double>> sin_0(TIME_SCALE_M_0, rows0);
    vector<vector<double>> sin_1(TIME_SCALE_M_0, rows0);
    vector<vector<double>> sin_2(TIME_SCALE_M_0, rows0);

    sin_0[0] = sin_inits0; // initial conditions for sin solution
    sin_1[0] = sin_inits1;
    sin_2[0] = sin_inits2;

    for (int n = 0; n < SPACE_SCALE_N_0; n++) { //////////// SIN 0, AVERAGE ///////////// STEP1
        sin_0[1][n] = sin_0[0][n] - TIME_SPLIT_T_0/SPACE_SPLIT_H_0 *(roll(sin_0[0], n + 1) + sin_0[0][n] + 
            roll(sin_0[0], n - 1))*(roll(sin_0[0], n + 1) - roll(sin_0[0], n - 1)) - TIME_SPLIT_T_0/(2*pow(SPACE_SPLIT_H_0, 3)) * 
            (roll(sin_0[0], n + 2) - 2*roll(sin_0[0], n + 1) + 2*roll(sin_0[0], n - 1) - roll(sin_0[0], n - 2));
    }

    for (int m = 1; m < TIME_SCALE_M_0 - 1; m++) { //////////// SIN 0, AVERAGE /////////////
        for (int n = 0; n < SPACE_SCALE_N_0; n++) {
            sin_0[m + 1][n] = sin_0[m - 1][n] - 2 * TIME_SPLIT_T_0/SPACE_SPLIT_H_0 * (roll(sin_0[m], n + 1) + sin_0[m][n] + 
                roll(sin_0[m], n - 1))*(roll(sin_0[m], n + 1) - roll(sin_0[m], n - 1)) - TIME_SPLIT_T_0/pow(SPACE_SPLIT_H_0, 3) *
                (roll(sin_0[m], n + 2) - 2*roll(sin_0[m], n + 1) + 2*roll(sin_0[m], n - 1) - roll(sin_0[m], n - 2));
        }
    }

    for (int n = 0; n < SPACE_SCALE_N_0; n++) { //////////// SIN 1, AVERAGE ///////////// STEP1
        sin_1[1][n] = sin_1[0][n] - TIME_SPLIT_T_1/SPACE_SPLIT_H_1 *(roll(sin_1[0], n + 1) + sin_1[0][n] + 
            roll(sin_1[0], n - 1))*(roll(sin_1[0], n + 1) - roll(sin_1[0], n - 1)) - TIME_SPLIT_T_1/(2*pow(SPACE_SPLIT_H_1, 3)) * 
            (roll(sin_1[0], n + 2) - 2*roll(sin_1[0], n + 1) + 2*roll(sin_1[0], n - 1) - roll(sin_1[0], n - 2));
    }

    for (int m = 1; m < TIME_SCALE_M_0 - 1; m++) { //////////// SIN 1, AVERAGE /////////////
        for (int n = 0; n < SPACE_SCALE_N_0; n++) {
            sin_1[m + 1][n] = sin_1[m - 1][n] - 2 * TIME_SPLIT_T_1/SPACE_SPLIT_H_1 * (roll(sin_1[m], n + 1) + sin_1[m][n] + 
                roll(sin_1[m], n - 1))*(roll(sin_1[m], n + 1) - roll(sin_1[m], n - 1)) - TIME_SPLIT_T_1/pow(SPACE_SPLIT_H_1, 3) *
                (roll(sin_1[m], n + 2) - 2*roll(sin_1[m], n + 1) + 2*roll(sin_1[m], n - 1) - roll(sin_1[m], n - 2));
        }
    }

    for (int n = 0; n < SPACE_SCALE_N_0; n++) { //////////// SIN 2, AVERAGE ///////////// STEP1
        sin_2[1][n] = sin_2[0][n] - TIME_SPLIT_T_2/SPACE_SPLIT_H_2 *(roll(sin_2[0], n + 1) + sin_2[0][n] + 
            roll(sin_2[0], n - 1))*(roll(sin_2[0], n + 1) - roll(sin_2[0], n - 1)) - TIME_SPLIT_T_2/(2*pow(SPACE_SPLIT_H_2, 3)) * 
            (roll(sin_2[0], n + 2) - 2*roll(sin_2[0], n + 1) + 2*roll(sin_2[0], n - 1) - roll(sin_2[0], n - 2));
    }

    for (int m = 1; m < TIME_SCALE_M_0 - 1; m++) { //////////// SIN 2, AVERAGE /////////////
        for (int n = 0; n < SPACE_SCALE_N_0; n++) {
            sin_2[m + 1][n] = sin_2[m - 1][n] - 2 * TIME_SPLIT_T_2/SPACE_SPLIT_H_2 * (roll(sin_2[m], n + 1) + sin_2[m][n] + 
                roll(sin_2[m], n - 1))*(roll(sin_2[m], n + 1) - roll(sin_2[m], n - 1)) - TIME_SPLIT_T_2/pow(SPACE_SPLIT_H_2, 3) *
                (roll(sin_2[m], n + 2) - 2*roll(sin_2[m], n + 1) + 2*roll(sin_2[m], n - 1) - roll(sin_2[m], n - 2));
        }
    }

    string line;
    ifstream in("input.txt");
    if (in.is_open()) {
        while (getline(in, line))
            cout << line << endl;
    }
    in.close();

    ofstream outsin0; // for recording
    outsin0.open("sin0 1000*1000 SPLIT 0.00001*0.12.txt"); // record sin-0 1000*1000 solution
    if (outsin0.is_open()) {
        outsin0 << TIME_SPLIT_T_0 << " " << SPACE_SPLIT_H_0 << " " << TIME_SCALE_M_0*TIME_SPLIT_T_0 << " " << SPACE_SCALE_N_0*SPACE_SPLIT_H_0 << endl;
        outsin0 << "{";
        for (int m = 0; m < TIME_SCALE_M_0 - 1; m++) {
            outsin0 << "{";
            for (int n = 0; n < SPACE_SCALE_N_0 - 1; n++) {
                outsin0 << sin_0[m][n] << ", ";
            }
            outsin0 << sin_0[m][SPACE_SCALE_N_0 - 1] << "}, ";
        }
        outsin0 << "{ ";
        for (int n = 0; n < SPACE_SCALE_N_0 - 1; n++) {
            outsin0 << sin_0[TIME_SCALE_M_0 - 1][n] << ", ";
        }
        outsin0 << sin_0[TIME_SCALE_M_0 - 1][SPACE_SCALE_N_0 - 1] << "}}";
    }
    outsin0.close();

    ofstream outsin1; // for recording
    outsin1.open("sin0 1000*1000 SPLIT 0.00005*0.15.txt"); // record sin-1 1000*1000 solution
    if (outsin1.is_open()) {
        outsin1 << TIME_SPLIT_T_1 << " " << SPACE_SPLIT_H_1 << " " << TIME_SCALE_M_0*TIME_SPLIT_T_1 << " " << SPACE_SCALE_N_0*SPACE_SPLIT_H_1 << endl;
        outsin1 << "{";
        for (int m = 0; m < TIME_SCALE_M_0 - 1; m++) {
            outsin1 << "{";
            for (int n = 0; n < SPACE_SCALE_N_0 - 1; n++) {
                outsin1 << sin_1[m][n] << ", ";
            }
            outsin1 << sin_1[m][SPACE_SCALE_N_0 - 1] << "}, ";
        }
        outsin1 << "{ ";
        for (int n = 0; n < SPACE_SCALE_N_0 - 1; n++) {
            outsin1 << sin_1[TIME_SCALE_M_0 - 1][n] << ", ";
        }
        outsin1 << sin_1[TIME_SCALE_M_0 - 1][SPACE_SCALE_N_0 - 1] << "}}";
    }
    outsin1.close();

    ofstream outsin2; // for recording
    outsin2.open("sin0 1000*1000 SPLIT 0.0001*0.3.txt"); // record sin-2 1000*1000 solution
    if (outsin2.is_open()) {
        outsin2 << TIME_SPLIT_T_2 << " " << SPACE_SPLIT_H_2 << " " << TIME_SCALE_M_0*TIME_SPLIT_T_2 << " " << SPACE_SCALE_N_0*SPACE_SPLIT_H_2 << endl;
        outsin2 << "{";
        for (int m = 0; m < TIME_SCALE_M_0 - 1; m++) {
            outsin2 << "{";
            for (int n = 0; n < SPACE_SCALE_N_0 - 1; n++) {
                outsin2 << sin_2[m][n] << ", ";
            }
            outsin2 << sin_2[m][SPACE_SCALE_N_0 - 1] << "}, ";
        }
        outsin2 << "{ ";
        for (int n = 0; n < SPACE_SCALE_N_0 - 1; n++) {
            outsin2 << sin_2[TIME_SCALE_M_0 - 1][n] << ", ";
        }
        outsin2 << sin_2[TIME_SCALE_M_0 - 1][SPACE_SCALE_N_0 - 1] << "}}";
    }
    outsin2.close();

    return 0;
}