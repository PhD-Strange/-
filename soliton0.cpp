#include "solution.h"

int main() {
    vector<double> soliton_inits0(SPACE_SCALE_N_0); // solitonic initial conditions
    vector<double> soliton_inits1(SPACE_SCALE_N_0);
    vector<double> soliton_inits2(SPACE_SCALE_N_0);

    double velocity = 0.01; // velocity = c > 0
    double a = 0; // shift
    
    cout << check_stability(velocity/2, TIME_SPLIT_T_0, SPACE_SPLIT_H_0) << endl;
    cout << check_stability(velocity/2, TIME_SPLIT_T_1, SPACE_SPLIT_H_1) << endl;
    cout << check_stability(velocity/2, TIME_SPLIT_T_2, SPACE_SPLIT_H_2) << endl;

    for (int i = 0; i < SPACE_SCALE_N_0; i++) { // init soliton with h = SPACE_SPLIT_0
        soliton_inits0.at(i) = 0.5 * velocity * pow(cosh(sqrt(velocity)*(i*SPACE_SPLIT_H_0 - a)/2), -2);
    }

    for (int i = 0; i < SPACE_SCALE_N_0; i++) { // init soliton with h = SPACE_SPLIT_1
        soliton_inits1.at(i) = 0.5 * velocity * pow(cosh(sqrt(velocity)*(i*SPACE_SPLIT_H_1 - a)/2), -2);
    }


    for (int i = 0; i < SPACE_SCALE_N_0; i++) { // init soliton with h = SPACE_SPLIT_2
        soliton_inits2.at(i) = 0.5 * velocity * pow(cosh(sqrt(velocity)*(i*SPACE_SPLIT_H_2 - a)/2), -2);
    }

    vector<double> rows0(SPACE_SCALE_N_0, 0);

    vector<vector<double>> soliton_0(TIME_SCALE_M_0, rows0);
    vector<vector<double>> soliton_1(TIME_SCALE_M_0, rows0);
    vector<vector<double>> soliton_2(TIME_SCALE_M_0, rows0);

    soliton_0[0] = soliton_inits0; // initial conditions for soliton solution
    soliton_1[0] = soliton_inits1;
    soliton_2[0] = soliton_inits2;

    for (int n = 0; n < SPACE_SCALE_N_0; n++) { //////////// SOLITON 0, AVERAGE ///////////// STEP1
        soliton_0[1][n] = soliton_0[0][n] - TIME_SPLIT_T_0/SPACE_SPLIT_H_0 *(roll(soliton_0[0], n + 1) + soliton_0[0][n] + 
            roll(soliton_0[0], n - 1))*(roll(soliton_0[0], n + 1) - roll(soliton_0[0], n - 1)) - TIME_SPLIT_T_0/(2*pow(SPACE_SPLIT_H_0, 3)) * 
            (roll(soliton_0[0], n + 2) - 2*roll(soliton_0[0], n + 1) + 2*roll(soliton_0[0], n - 1) - roll(soliton_0[0], n - 2));
    }

    for (int m = 1; m < TIME_SCALE_M_0 - 1; m++) { //////////// SOLITON 0, AVERAGE /////////////
        for (int n = 0; n < SPACE_SCALE_N_0; n++) {
            soliton_0[m + 1][n] = soliton_0[m - 1][n] - 2 * TIME_SPLIT_T_0/SPACE_SPLIT_H_0 * (roll(soliton_0[m], n + 1) + soliton_0[m][n] + 
                roll(soliton_0[m], n - 1))*(roll(soliton_0[m], n + 1) - roll(soliton_0[m], n - 1)) - TIME_SPLIT_T_0/pow(SPACE_SPLIT_H_0, 3) *
                (roll(soliton_0[m], n + 2) - 2*roll(soliton_0[m], n + 1) + 2*roll(soliton_0[m], n - 1) - roll(soliton_0[m], n - 2));
        }
    }

    for (int n = 0; n < SPACE_SCALE_N_0; n++) { //////////// SOLITON 1, AVERAGE ///////////// STEP1
        soliton_1[1][n] = soliton_1[0][n] - TIME_SPLIT_T_1/SPACE_SPLIT_H_1 *(roll(soliton_1[0], n + 1) + soliton_1[0][n] + 
            roll(soliton_1[0], n - 1))*(roll(soliton_1[0], n + 1) - roll(soliton_1[0], n - 1)) - TIME_SPLIT_T_1/(2*pow(SPACE_SPLIT_H_1, 3)) * 
            (roll(soliton_1[0], n + 2) - 2*roll(soliton_1[0], n + 1) + 2*roll(soliton_1[0], n - 1) - roll(soliton_1[0], n - 2));
    }

    for (int m = 1; m < TIME_SCALE_M_0 - 1; m++) { //////////// SOLITON 1, AVERAGE /////////////
        for (int n = 0; n < SPACE_SCALE_N_0; n++) {
            soliton_1[m + 1][n] = soliton_1[m - 1][n] - 2 * TIME_SPLIT_T_1/SPACE_SPLIT_H_1 * (roll(soliton_1[m], n + 1) + soliton_1[m][n] + 
                roll(soliton_1[m], n - 1))*(roll(soliton_1[m], n + 1) - roll(soliton_1[m], n - 1)) - TIME_SPLIT_T_1/pow(SPACE_SPLIT_H_1, 3) *
                (roll(soliton_1[m], n + 2) - 2*roll(soliton_1[m], n + 1) + 2*roll(soliton_1[m], n - 1) - roll(soliton_1[m], n - 2));
        }
    }

    for (int n = 0; n < SPACE_SCALE_N_0; n++) { //////////// SOLITON 2, AVERAGE ///////////// STEP1
        soliton_2[1][n] = soliton_2[0][n] - TIME_SPLIT_T_2/SPACE_SPLIT_H_2 *(roll(soliton_2[0], n + 1) + soliton_2[0][n] + 
            roll(soliton_2[0], n - 1))*(roll(soliton_2[0], n + 1) - roll(soliton_2[0], n - 1)) - TIME_SPLIT_T_2/(2*pow(SPACE_SPLIT_H_2, 3)) * 
            (roll(soliton_2[0], n + 2) - 2*roll(soliton_2[0], n + 1) + 2*roll(soliton_2[0], n - 1) - roll(soliton_2[0], n - 2));
    }

    for (int m = 1; m < TIME_SCALE_M_0 - 1; m++) { //////////// SOLITON 2, AVERAGE /////////////
        for (int n = 0; n < SPACE_SCALE_N_0; n++) {
            soliton_2[m + 1][n] = soliton_2[m - 1][n] - 2 * TIME_SPLIT_T_2/SPACE_SPLIT_H_2 * (roll(soliton_2[m], n + 1) + soliton_2[m][n] + 
                roll(soliton_2[m], n - 1))*(roll(soliton_2[m], n + 1) - roll(soliton_2[m], n - 1)) - TIME_SPLIT_T_2/pow(SPACE_SPLIT_H_2, 3) *
                (roll(soliton_2[m], n + 2) - 2*roll(soliton_2[m], n + 1) + 2*roll(soliton_2[m], n - 1) - roll(soliton_2[m], n - 2));
        }
    }

    string line;
    ifstream in("input.txt");
    if (in.is_open()) {
        while (getline(in, line))
            cout << line << endl;
    }
    in.close();

    ofstream out0; // for recording
    out0.open("soliton0 1000*1000 SPLIT 0.00001*0.12.txt"); // record soliton-0 1000*1000 solution
    if (out0.is_open()) {
        out0 << TIME_SPLIT_T_0 << " " << SPACE_SPLIT_H_0 << " " << TIME_SCALE_M_0*TIME_SPLIT_T_0 << " " << SPACE_SCALE_N_0*SPACE_SPLIT_H_0 << endl;
        out0 << "{";
        for (int m = 0; m < TIME_SCALE_M_0 - 1; m++) {
            out0 << "{";
            for (int n = 0; n < SPACE_SCALE_N_0 - 1; n++) {
                out0 << soliton_0[m][n] << ", ";
            }
            out0 << soliton_0[m][SPACE_SCALE_N_0 - 1] << "}, ";
        }
        out0 << "{ ";
        for (int n = 0; n < SPACE_SCALE_N_0 - 1; n++) {
            out0 << soliton_0[TIME_SCALE_M_0 - 1][n] << ", ";
        }
        out0 << soliton_0[TIME_SCALE_M_0 - 1][SPACE_SCALE_N_0 - 1] << "}}";
    }
    out0.close();

    ofstream out1;
    out1.open("soliton1 1000*1000 SPLIT 0.00005*0.15.txt"); // record soliton-1 1000*1000 solution
    if (out1.is_open()) {
        out1 << TIME_SPLIT_T_1 << " " << SPACE_SPLIT_H_1 << " " << TIME_SCALE_M_0*TIME_SPLIT_T_1 << " " << SPACE_SCALE_N_0*SPACE_SPLIT_H_1 << endl;
        out1 << "{";
        for (int m = 0; m < TIME_SCALE_M_0 - 1; m++) {
            out1 << "{";
            for (int n = 0; n < SPACE_SCALE_N_0 - 1; n++) {
                out1 << soliton_1[m][n] << ", ";
            }
            out1 << soliton_1[m][SPACE_SCALE_N_0 - 1] << "}, ";
        }
        out1 << "{ ";
        for (int n = 0; n < SPACE_SCALE_N_0 - 1; n++) {
            out1 << soliton_1[TIME_SCALE_M_0 - 1][n] << ", ";
        }
        out1 << soliton_1[TIME_SCALE_M_0 - 1][SPACE_SCALE_N_0 - 1] << "}}";
    }
    out1.close();

    ofstream out2;
    out2.open("soliton2 1000*1000 SPLIT 0.0001*0.3"); // record soliton-2 1000*1000 soliton
    if (out2.is_open()) {
        out2 << TIME_SPLIT_T_2 << " " << SPACE_SPLIT_H_2 << " " << TIME_SCALE_M_0*TIME_SPLIT_T_2 << " " << SPACE_SCALE_N_0*SPACE_SPLIT_H_2 << endl;
        out2 << "{";
        for (int m = 0; m < TIME_SCALE_M_0 - 1; m++) {
            out2 << "{";
            for (int n = 0; n < SPACE_SCALE_N_0 - 1; n++) {
                out2 << soliton_2[m][n] << ", ";
            }
            out2 << soliton_2[m][SPACE_SCALE_N_0 - 1] << "}, ";
        }
        out2 << "{ ";
        for (int n = 0; n < SPACE_SCALE_N_0 - 1; n++) {
            out2 << soliton_2[TIME_SCALE_M_0 - 1][n] << ", ";
        }
        out2 << soliton_2[TIME_SCALE_M_0 - 1][SPACE_SCALE_N_0 - 1] << "}}";
    }
    out2.close();

    return 0;
}