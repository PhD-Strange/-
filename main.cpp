#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
using namespace std;

const int TIME_SCALE_M_0 = 1000;
const int SPACE_SCALE_N_0 = 1000;

const int TIME_SCALE_M_MAX = 10000;
const int SPACE_SCALE_N_MAX = 10000;

const double TIME_SPLIT_T_0 = 0.00001; // t0
const double SPACE_SPLIT_H_0 = 0.12; // h0

const double TIME_SPLIT_T_1 = 0.00005; // t1
const double SPACE_SPLIT_H_1 = 0.15; // h1

const double TIME_SPLIT_T_2 = 0.0001; // t2
const double SPACE_SPLIT_H_2 = 0.3; // h2

double roll(vector<double> v, unsigned int n) {
    return v[n % v.size()];
}

bool check_stability(double max, const double t, const double h) {
    return t/h * abs(-2 * max + pow(h, -2)) < 2/pow(3, 1.5);
}

int main() {
    vector<double> soliton_inits0(SPACE_SCALE_N_0); // solitonic initial conditions
    vector<double> soliton_inits1(SPACE_SCALE_N_0);
    vector<double> soliton_inits2(SPACE_SCALE_N_0);
    vector<double> soliton_inits0_max(SPACE_SCALE_N_MAX);
    vector<double> soliton_inits1_max(SPACE_SCALE_N_MAX);
    vector<double> soliton_inits2_max(SPACE_SCALE_N_MAX);

    vector<double> sin_inits0(SPACE_SCALE_N_0); // sin initial conditions
    vector<double> sin_inits1(SPACE_SCALE_N_0);
    vector<double> sin_inits2(SPACE_SCALE_N_0);
    vector<double> sin_inits0_max(SPACE_SCALE_N_MAX);
    vector<double> sin_inits1_max(SPACE_SCALE_N_MAX);
    vector<double> sin_inits2_max(SPACE_SCALE_N_MAX);

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

    for (int i = 0; i < SPACE_SCALE_N_MAX; i++) { // init soliton with h = SPACE_SPLIT_0, MAX
        soliton_inits0_max.at(i) = 0.5 * velocity * pow(cosh(sqrt(velocity)*(i*SPACE_SPLIT_H_0 - a)/2), -2);
    }

    for (int i = 0; i < SPACE_SCALE_N_MAX; i++) { // init soliton with h = SPACE_SPLIT_1, MAX
        soliton_inits1_max.at(i) = 0.5 * velocity * pow(cosh(sqrt(velocity)*(i*SPACE_SPLIT_H_1 - a)/2), -2);
    }

    for (int i = 0; i < SPACE_SCALE_N_MAX; i++) { // init soliton with h = SPACE_SPLIT_2, MAX
        soliton_inits2_max.at(i) = 0.5 * velocity * pow(cosh(sqrt(velocity)*(i*SPACE_SPLIT_H_2 - a)/2), -2);
    }

    for (int i = 0; i < SPACE_SCALE_N_0; i++) { // init sin with h = SPACE_SPLIT_0
        sin_inits0.at(i) = sin(M_PI * SPACE_SPLIT_H_0 * i);
    }

    for (int i = 0; i < SPACE_SCALE_N_0; i++) { // init sin with h = SPACE_SPLIT_1
        sin_inits1.at(i) = sin(M_PI * SPACE_SPLIT_H_1 * i);
    }

    for (int i = 0; i < SPACE_SCALE_N_0; i++) { // init sin with h = SPACE_SPLIT_2
        sin_inits2.at(i) = sin(M_PI * SPACE_SPLIT_H_2 * i);
    }

    for (int i = 0; i < SPACE_SCALE_N_MAX; i++) { // init sin with h = SPACE_SPLIT_0, MAX
        sin_inits0_max.at(i) = sin(M_PI * SPACE_SPLIT_H_0 * i);
    }

    for (int i = 0; i < SPACE_SCALE_N_MAX; i++) { // init sin with h = SPACE_SPLIT_1, MAX
        sin_inits1_max.at(i) = sin(M_PI * SPACE_SPLIT_H_1 * i);
    }

    for (int i = 0; i < SPACE_SCALE_N_MAX; i++) { // init sin with h = SPACE_SPLIT_2, MAX
        sin_inits2_max.at(i) = sin(M_PI * SPACE_SPLIT_H_2 * i);
    }

    vector<double> rows0(SPACE_SCALE_N_0, 0);
    vector<double> rowsMax(SPACE_SCALE_N_MAX, 0);

    vector<vector<double>> soliton_0(TIME_SCALE_M_0, rows0);
    vector<vector<double>> soliton_1(TIME_SCALE_M_0, rows0);
    vector<vector<double>> soliton_2(TIME_SCALE_M_0, rows0);
    vector<vector<double>> soliton_0_max(TIME_SCALE_M_MAX, rowsMax);
    vector<vector<double>> soliton_1_max(TIME_SCALE_M_MAX, rowsMax);
    vector<vector<double>> soliton_2_max(TIME_SCALE_M_MAX, rowsMax);
    
    vector<vector<double>> sin_0(TIME_SCALE_M_0, rows0);
    vector<vector<double>> sin_1(TIME_SCALE_M_0, rows0);
    vector<vector<double>> sin_2(TIME_SCALE_M_0, rows0);
    vector<vector<double>> sin_0_max(TIME_SCALE_M_MAX, rowsMax);
    vector<vector<double>> sin_1_max(TIME_SCALE_M_MAX, rowsMax);
    vector<vector<double>> sin_2_max(TIME_SCALE_M_MAX, rowsMax);

    soliton_0[0] = soliton_inits0; // initial conditions for soliton solution
    soliton_1[0] = soliton_inits1;
    soliton_2[0] = soliton_inits2;
    soliton_0_max[0] = soliton_inits0_max;
    soliton_1_max[0] = soliton_inits1_max;
    soliton_2_max[0] = soliton_inits2_max;

    sin_0[0] = sin_inits0; // initial conditions for sin solution
    sin_1[0] = sin_inits1;
    sin_2[0] = sin_inits2;
    sin_0_max[0] = sin_inits0_max;
    sin_1_max[0] = sin_inits1_max;
    sin_2_max[0] = sin_inits2_max;

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

    for (int n = 0; n < SPACE_SCALE_N_MAX; n++) { //////////// SOLITON 1, MAXIMA ///////////// STEP1
        soliton_1_max[1][n] = soliton_1_max[0][n] - TIME_SPLIT_T_1/SPACE_SPLIT_H_1 *(roll(soliton_1_max[0], n + 1) + soliton_1_max[0][n] + 
            roll(soliton_1_max[0], n - 1))*(roll(soliton_1_max[0], n + 1) - roll(soliton_1_max[0], n - 1)) - TIME_SPLIT_T_1/(2*pow(SPACE_SPLIT_H_1, 3)) * 
            (roll(soliton_1_max[0], n + 2) - 2*roll(soliton_1_max[0], n + 1) + 2*roll(soliton_1_max[0], n - 1) - roll(soliton_1_max[0], n - 2));
    }

    for (int m = 1; m < TIME_SCALE_M_MAX - 1; m++) { //////////// SOLITON 1, MAXIMA /////////////
        for (int n = 0; n < SPACE_SCALE_N_MAX; n++) {
            soliton_1_max[m + 1][n] = soliton_1_max[m - 1][n] - 2 * TIME_SPLIT_T_1/SPACE_SPLIT_H_1 * (roll(soliton_1_max[m], n + 1) + soliton_1_max[m][n] + 
                roll(soliton_1_max[m], n - 1))*(roll(soliton_1_max[m], n + 1) - roll(soliton_1_max[m], n - 1)) - TIME_SPLIT_T_1/pow(SPACE_SPLIT_H_1, 3) *
                (roll(soliton_1_max[m], n + 2) - 2*roll(soliton_1_max[m], n + 1) + 2*roll(soliton_1_max[m], n - 1) - roll(soliton_1_max[m], n - 2));
        }
    }

    for (int n = 0; n < SPACE_SCALE_N_MAX; n++) { //////////// SOLITON 2, MAXIMA ///////////// STEP1
        soliton_2_max[1][n] = soliton_2_max[0][n] - TIME_SPLIT_T_2/SPACE_SPLIT_H_2 *(roll(soliton_2_max[0], n + 1) + soliton_2_max[0][n] + 
            roll(soliton_2_max[0], n - 1))*(roll(soliton_2_max[0], n + 1) - roll(soliton_2_max[0], n - 1)) - TIME_SPLIT_T_2/(2*pow(SPACE_SPLIT_H_2, 3)) * 
            (roll(soliton_2_max[0], n + 2) - 2*roll(soliton_2_max[0], n + 1) + 2*roll(soliton_2_max[0], n - 1) - roll(soliton_2_max[0], n - 2));
    }

    for (int m = 1; m < TIME_SCALE_M_MAX - 1; m++) { //////////// SOLITON 2, MAXIMA /////////////
        for (int n = 0; n < SPACE_SCALE_N_MAX; n++) {
            soliton_2_max[m + 1][n] = soliton_2_max[m - 1][n] - 2 * TIME_SPLIT_T_2/SPACE_SPLIT_H_2 * (roll(soliton_2_max[m], n + 1) + soliton_2_max[m][n] + 
                roll(soliton_2_max[m], n - 1))*(roll(soliton_2_max[m], n + 1) - roll(soliton_2_max[m], n - 1)) - TIME_SPLIT_T_2/pow(SPACE_SPLIT_H_2, 3) *
                (roll(soliton_2_max[m], n + 2) - 2*roll(soliton_2_max[m], n + 1) + 2*roll(soliton_2_max[m], n - 1) - roll(soliton_2_max[m], n - 2));
        }
    }

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

    for (int n = 0; n < SPACE_SCALE_N_MAX; n++) { //////////// SIN 0, MAXIMA ///////////// STEP1
        sin_0_max[1][n] = sin_0_max[0][n] - TIME_SPLIT_T_0/SPACE_SPLIT_H_0 *(roll(sin_0_max[0], n + 1) + sin_0_max[0][n] + 
            roll(sin_0_max[0], n - 1))*(roll(sin_0_max[0], n + 1) - roll(sin_0_max[0], n - 1)) - TIME_SPLIT_T_0/(2*pow(SPACE_SPLIT_H_0, 3)) * 
            (roll(sin_0_max[0], n + 2) - 2*roll(sin_0_max[0], n + 1) + 2*roll(sin_0_max[0], n - 1) - roll(sin_0_max[0], n - 2));
    }

    for (int m = 1; m < TIME_SCALE_M_MAX - 1; m++) { //////////// SIN 0, MAXIMA /////////////
        for (int n = 0; n < SPACE_SCALE_N_MAX; n++) {
            sin_0_max[m + 1][n] = sin_0_max[m - 1][n] - 2 * TIME_SPLIT_T_0/SPACE_SPLIT_H_0 * (roll(sin_0_max[m], n + 1) + sin_0_max[m][n] + 
                roll(sin_0_max[m], n - 1))*(roll(sin_0_max[m], n + 1) - roll(sin_0_max[m], n - 1)) - TIME_SPLIT_T_0/pow(SPACE_SPLIT_H_0, 3) *
                (roll(sin_0_max[m], n + 2) - 2*roll(sin_0_max[m], n + 1) + 2*roll(sin_0_max[m], n - 1) - roll(sin_0_max[m], n - 2));
        }
    }

    for (int n = 0; n < SPACE_SCALE_N_MAX; n++) { //////////// SIN 1, MAXIMA ///////////// STEP1
        sin_1_max[1][n] = sin_1_max[0][n] - TIME_SPLIT_T_1/SPACE_SPLIT_H_1 *(roll(sin_1_max[0], n + 1) + sin_1_max[0][n] + 
            roll(sin_1_max[0], n - 1))*(roll(sin_1_max[0], n + 1) - roll(sin_1_max[0], n - 1)) - TIME_SPLIT_T_1/(2*pow(SPACE_SPLIT_H_1, 3)) * 
            (roll(sin_1_max[0], n + 2) - 2*roll(sin_1_max[0], n + 1) + 2*roll(sin_1_max[0], n - 1) - roll(sin_1_max[0], n - 2));
    }

    for (int m = 1; m < TIME_SCALE_M_MAX - 1; m++) { //////////// SIN 1, MAXIMA /////////////
        for (int n = 0; n < SPACE_SCALE_N_MAX; n++) {
            sin_1_max[m + 1][n] = sin_1_max[m - 1][n] - 2 * TIME_SPLIT_T_1/SPACE_SPLIT_H_1 * (roll(sin_1_max[m], n + 1) + sin_1_max[m][n] + 
                roll(sin_1_max[m], n - 1))*(roll(sin_1_max[m], n + 1) - roll(sin_1_max[m], n - 1)) - TIME_SPLIT_T_1/pow(SPACE_SPLIT_H_1, 3) *
                (roll(sin_1_max[m], n + 2) - 2*roll(sin_1_max[m], n + 1) + 2*roll(sin_1_max[m], n - 1) - roll(sin_1_max[m], n - 2));
        }
    }

    for (int n = 0; n < SPACE_SCALE_N_MAX; n++) { //////////// SIN 2, MAXIMA ///////////// STEP1
        sin_2_max[1][n] = sin_2_max[0][n] - TIME_SPLIT_T_2/SPACE_SPLIT_H_2 *(roll(sin_2_max[0], n + 1) + sin_2_max[0][n] + 
            roll(sin_2_max[0], n - 1))*(roll(sin_2_max[0], n + 1) - roll(sin_2_max[0], n - 1)) - TIME_SPLIT_T_2/(2*pow(SPACE_SPLIT_H_2, 3)) * 
            (roll(sin_2_max[0], n + 2) - 2*roll(sin_2_max[0], n + 1) + 2*roll(sin_2_max[0], n - 1) - roll(sin_2_max[0], n - 2));
    }

    for (int m = 1; m < TIME_SCALE_M_MAX - 1; m++) { //////////// SIN 2, MAXIMA /////////////
        for (int n = 0; n < SPACE_SCALE_N_MAX; n++) {
            sin_2_max[m + 1][n] = sin_2_max[m - 1][n] - 2 * TIME_SPLIT_T_2/SPACE_SPLIT_H_2 * (roll(sin_2_max[m], n + 1) + sin_2_max[m][n] + 
                roll(sin_2_max[m], n - 1))*(roll(sin_2_max[m], n + 1) - roll(sin_2_max[m], n - 1)) - TIME_SPLIT_T_2/pow(SPACE_SPLIT_H_2, 3) *
                (roll(sin_2_max[m], n + 2) - 2*roll(sin_2_max[m], n + 1) + 2*roll(sin_2_max[m], n - 1) - roll(sin_2_max[m], n - 2));
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

    ofstream out1_max;
    out1_max.open("soliton1 10000*10000 SPLIT 0.00005*0.15.txt"); // record soliton-1 10000*10000 solution
    if (out1_max.is_open()) {
        out1_max << TIME_SPLIT_T_1 << " " << SPACE_SPLIT_H_1 << " " << TIME_SCALE_M_MAX*TIME_SPLIT_T_1 << " " << SPACE_SCALE_N_MAX*SPACE_SPLIT_H_1 << endl;
        out1_max << "{";
        for (int m = 0; m < TIME_SCALE_M_MAX - 1; m++) {
            out1_max << "{";
            for (int n = 0; n < SPACE_SCALE_N_MAX - 1; n++) {
                out1_max << soliton_1_max[m][n] << ", ";
            }
            out1_max << soliton_1_max[m][SPACE_SCALE_N_MAX - 1] << "}, ";
        }
        out1_max << "{ ";
        for (int n = 0; n < SPACE_SCALE_N_MAX - 1; n++) {
            out1_max << soliton_1_max[TIME_SCALE_M_MAX - 1][n] << ", ";
        }
        out1_max << soliton_1_max[TIME_SCALE_M_MAX - 1][SPACE_SCALE_N_MAX - 1] << "}}";
    }
    out1_max.close();

    ofstream out2_max;
    out2_max.open("soliton2 10000*10000 SPLIT 0.0001*0.3"); // record soliton-2 10000*10000 soliton
    if (out2_max.is_open()) {
        out2_max << TIME_SPLIT_T_2 << " " << SPACE_SPLIT_H_2 << " " << TIME_SCALE_M_MAX*TIME_SPLIT_T_2 << " " << SPACE_SCALE_N_MAX*SPACE_SPLIT_H_2 << endl;
        out2_max << "{";
        for (int m = 0; m < TIME_SCALE_M_MAX - 1; m++) {
            out2_max << "{";
            for (int n = 0; n < SPACE_SCALE_N_MAX - 1; n++) {
                out2_max << soliton_2_max[m][n] << ", ";
            }
            out2_max << soliton_2_max[m][SPACE_SCALE_N_MAX - 1] << "}, ";
        }
        out2_max << "{ ";
        for (int n = 0; n < SPACE_SCALE_N_MAX - 1; n++) {
            out2_max << soliton_2_max[TIME_SCALE_M_MAX - 1][n] << ", ";
        }
        out2_max << soliton_2_max[TIME_SCALE_M_MAX - 1][SPACE_SCALE_N_MAX - 1] << "}}";
    }
    out2_max.close();

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
    outsin2.open("soliton0 1000*1000 SPLIT 0.0001*0.3.txt"); // record sin-2 1000*1000 solution
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

    ofstream outsin0_max; // for recording
    outsin0_max.open("sin0 10000*10000 SPLIT 0.00001*0.12.txt"); // record sin-0 10000*10000 solution
    if (outsin0_max.is_open()) {
        outsin0_max << TIME_SPLIT_T_0 << " " << SPACE_SPLIT_H_0 << " " << TIME_SCALE_M_MAX*TIME_SPLIT_T_0 << " " << SPACE_SCALE_N_MAX*SPACE_SPLIT_H_0 << endl;
        outsin0_max << "{";
        for (int m = 0; m < TIME_SCALE_M_MAX - 1; m++) {
            outsin0_max << "{";
            for (int n = 0; n < SPACE_SCALE_N_MAX - 1; n++) {
                outsin0_max << sin_0_max[m][n] << ", ";
            }
            outsin0_max << sin_0_max[m][SPACE_SCALE_N_MAX - 1] << "}, ";
        }
        outsin0_max << "{ ";
        for (int n = 0; n < SPACE_SCALE_N_MAX - 1; n++) {
            outsin0_max << sin_0_max[TIME_SCALE_M_MAX - 1][n] << ", ";
        }
        outsin0_max << sin_0_max[TIME_SCALE_M_MAX - 1][SPACE_SCALE_N_MAX - 1] << "}}";
    }
    outsin0_max.close();

    ofstream outsin1_max; // for recording
    outsin1_max.open("sin0 10000*10000 SPLIT 0.00005*0.15.txt"); // record sin-1 10000*10000 solution
    if (outsin1_max.is_open()) {
        outsin1_max << TIME_SPLIT_T_1 << " " << SPACE_SPLIT_H_1 << " " << TIME_SCALE_M_MAX*TIME_SPLIT_T_1 << " " << SPACE_SCALE_N_MAX*SPACE_SPLIT_H_1 << endl;
        outsin1_max << "{";
        for (int m = 0; m < TIME_SCALE_M_MAX - 1; m++) {
            outsin1_max << "{";
            for (int n = 0; n < SPACE_SCALE_N_MAX - 1; n++) {
                outsin1_max << sin_1_max[m][n] << ", ";
            }
            outsin1_max << sin_1_max[m][SPACE_SCALE_N_MAX - 1] << "}, ";
        }
        outsin1_max << "{ ";
        for (int n = 0; n < SPACE_SCALE_N_MAX - 1; n++) {
            outsin1_max << sin_1_max[TIME_SCALE_M_MAX - 1][n] << ", ";
        }
        outsin1_max << sin_1_max[TIME_SCALE_M_MAX - 1][SPACE_SCALE_N_MAX - 1] << "}}";
    }
    outsin1_max.close();

    ofstream outsin2_max; // for recording
    outsin2_max.open("soliton0 10000*10000 SPLIT 0.0001*0.3.txt"); // record sin-2 10000*10000 solution
    if (outsin2_max.is_open()) {
        outsin2_max << TIME_SPLIT_T_2 << " " << SPACE_SPLIT_H_2 << " " << TIME_SCALE_M_MAX*TIME_SPLIT_T_2 << " " << SPACE_SCALE_N_MAX*SPACE_SPLIT_H_2 << endl;
        outsin2_max << "{";
        for (int m = 0; m < TIME_SCALE_M_MAX - 1; m++) {
            outsin2_max << "{";
            for (int n = 0; n < SPACE_SCALE_N_MAX - 1; n++) {
                outsin2_max << sin_2_max[m][n] << ", ";
            }
            outsin2_max << sin_2_max[m][SPACE_SCALE_N_MAX - 1] << "}, ";
        }
        outsin2_max << "{ ";
        for (int n = 0; n < SPACE_SCALE_N_MAX - 1; n++) {
            outsin2_max << sin_2_max[TIME_SCALE_M_MAX - 1][n] << ", ";
        }
        outsin2_max << sin_2_max[TIME_SCALE_M_MAX - 1][SPACE_SCALE_N_MAX - 1] << "}}";
    }
    outsin2_max.close();

    return 0;
};