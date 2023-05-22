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