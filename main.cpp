#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
using namespace std;

const int TIME_SCALE_M = 100;
const int SPACE_SCALE_N = 100;
const int TIME_SPLIT_T = 0.01;
const int SPACE_SPLIT_H = 0.03;

double roll(vector<double>& v, unsigned int n) {
    return v[n % v.size()];
}

int main() {
    vector<double> inits(SPACE_SCALE_N); // initial conditions
    double velocity; // velocity = c > 0
    double a; // shift

    for (int i = 0; i < SPACE_SCALE_N; i++) {
        inits.at(i) = 0.5 * velocity * pow(cosh(sqrt(velocity)*(i*SPACE_SPLIT_H - a)/2), -2);
    }

    vector<vector<double>> solution;



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