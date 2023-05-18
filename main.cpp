#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
using namespace std;

const int TIME_SCALE_M = 100;
const int SPACE_SCALE_N = 100;
const double TIME_SPLIT_T = 0.00001;
const double SPACE_SPLIT_H = 0.12;

double roll(vector<double>& v, unsigned int n) {
    return v[n % v.size()];
}

bool check_stability(double max) {
    return TIME_SPLIT_T/SPACE_SPLIT_H * abs(-2 * max + pow(SPACE_SPLIT_H, -2)) < 2/pow(3, 1.5);
}

int main() {
    vector<double> inits(SPACE_SCALE_N); // initial conditions
    double velocity = 0.01; // velocity = c > 0
    double a = 0; // shift

    cout << check_stability(velocity/2) << endl;

    for (int i = 0; i < SPACE_SCALE_N; i++) {
        inits.at(i) = 0.5 * velocity * pow(cosh(sqrt(velocity)*(i*SPACE_SPLIT_H - a)/2), -2);
    }

    vector<double> rows(SPACE_SCALE_N, 0);

    vector<vector<double>> u(TIME_SCALE_M, rows);

    u.at(0) = inits;

    for (int n = 0; n < SPACE_SCALE_N; n++) {
        u[1][n] = u[0][n] - TIME_SPLIT_T/SPACE_SPLIT_H *(roll(u[0], n + 1) + u[0][n] + roll(u[0], n - 1))*(roll(u[0], n + 1) - roll(u[0], n - 1)) - TIME_SPLIT_T/(2*pow(SPACE_SPLIT_H, 3)) * (roll(u[0], n + 2) - 2*roll(u[0], n + 1) + 2*roll(u[0], n - 1) - roll(u[0], n - 2));
    }

    for (int m = 1; m < TIME_SCALE_M; m++) {
        for (int n = 0; n < SPACE_SCALE_N; n++) {
            u[m + 1][n] = u[m - 1][n] - 2 * TIME_SPLIT_T/SPACE_SPLIT_H *(roll(u[m], n + 1) + u[m][n] + roll(u[m], n - 1))*(roll(u[m], n + 1) - roll(u[m], n - 1)) - TIME_SPLIT_T/pow(SPACE_SPLIT_H, 3) * (roll(u[m], n + 2) - 2*roll(u[m], n + 1) + 2*roll(u[m], n - 1) - roll(u[m], n - 2));
        }
    }


    string line;
    ifstream in("input.txt");
    if (in.is_open()) {
        while (getline(in, line))
            cout << line << endl;
    }
    in.close();

    ofstream out; // for recording
    out.open("output.txt");
    if (out.is_open())
        out << "Hello World?" << endl;
    out.close();

    return 0;
};