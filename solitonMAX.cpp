#include "solution.h"

int main() {
    vector<double> soliton_inits0_max(SPACE_SCALE_N_MAX);

    vector<double> sin_inits0_max(SPACE_SCALE_N_MAX);

    double velocity = 0.01; // velocity = c > 0
    double a = 0; // shift

    cout << check_stability(velocity/2, TIME_SPLIT_T_0, SPACE_SPLIT_H_0) << endl;


    for (int i = 0; i < SPACE_SCALE_N_MAX; i++) { // init soliton with h = SPACE_SPLIT_0, MAX
        soliton_inits0_max.at(i) = 0.5 * velocity * pow(cosh(sqrt(velocity)*(i*SPACE_SPLIT_H_0 - a)/2), -2);
    }

    vector<double> rowsMax(SPACE_SCALE_N_MAX, 0);

    vector<vector<double>> soliton_0_max(TIME_SCALE_M_MAX, rowsMax);

    soliton_0_max[0] = soliton_inits0_max;

    for (int n = 0; n < SPACE_SCALE_N_MAX; n++) { //////////// SOLITON 0, MAXIMA ///////////// STEP1
        soliton_0_max[1][n] = soliton_0_max[0][n] - TIME_SPLIT_T_0/SPACE_SPLIT_H_0 *(roll(soliton_0_max[0], n + 1) + soliton_0_max[0][n] + 
            roll(soliton_0_max[0], n - 1))*(roll(soliton_0_max[0], n + 1) - roll(soliton_0_max[0], n - 1)) - TIME_SPLIT_T_0/(2*pow(SPACE_SPLIT_H_0, 3)) * 
            (roll(soliton_0_max[0], n + 2) - 2*roll(soliton_0_max[0], n + 1) + 2*roll(soliton_0_max[0], n - 1) - roll(soliton_0_max[0], n - 2));
    }

    for (int m = 1; m < TIME_SCALE_M_MAX - 1; m++) { //////////// SOLITON 0, MAXIMA /////////////
        for (int n = 0; n < SPACE_SCALE_N_MAX; n++) {
            soliton_0_max[m + 1][n] = soliton_0_max[m - 1][n] - 2 * TIME_SPLIT_T_0/SPACE_SPLIT_H_0 * (roll(soliton_0_max[m], n + 1) + soliton_0_max[m][n] + 
                roll(soliton_0_max[m], n - 1))*(roll(soliton_0_max[m], n + 1) - roll(soliton_0_max[m], n - 1)) - TIME_SPLIT_T_0/pow(SPACE_SPLIT_H_0, 3) *
                (roll(soliton_0_max[m], n + 2) - 2*roll(soliton_0_max[m], n + 1) + 2*roll(soliton_0_max[m], n - 1) - roll(soliton_0_max[m], n - 2));
        }
    }

    string line;
    ifstream in("input.txt");
    if (in.is_open()) {
        while (getline(in, line))
            cout << line << endl;
    }
    in.close();

    ofstream out0_max; // for recording
    out0_max.open("soliton0 10000*10000 SPLIT 0.00001*0.12.txt"); // record soliton-0 10000*10000 solution
    if (out0_max.is_open()) {
        out0_max << TIME_SPLIT_T_0 << " " << SPACE_SPLIT_H_0 << " " << TIME_SCALE_M_MAX*TIME_SPLIT_T_0 << " " << SPACE_SCALE_N_MAX*SPACE_SPLIT_H_0 << endl;
        out0_max << "{";
        for (int m = 0; m < TIME_SCALE_M_MAX - 1; m++) {
            out0_max << "{";
            for (int n = 0; n < SPACE_SCALE_N_MAX - 1; n++) {
                out0_max << soliton_0_max[m][n] << ", ";
            }
            out0_max << soliton_0_max[m][SPACE_SCALE_N_MAX - 1] << "}, ";
        }
        out0_max << "{ ";
        for (int n = 0; n < SPACE_SCALE_N_MAX - 1; n++) {
            out0_max << soliton_0_max[TIME_SCALE_M_MAX - 1][n] << ", ";
        }
        out0_max << soliton_0_max[TIME_SCALE_M_MAX - 1][SPACE_SCALE_N_MAX - 1] << "}}";
    }
    out0_max.close();

    return 0;
}