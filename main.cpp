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

class Node
{
private:
    int m_; // u{m, n} = u(mT, nH)
    int n_;
    double number;
public:
    Node(int m, int n); // constructor
    void setM(int m);
    void setN(int n);
    void setNumber(double z);
    double getNumber();
};

Node::Node(int m, int n) {
    if (n < 0) {
        n_ = n + SPACE_SCALE_N;
    } else if (n > SPACE_SCALE_N) {
        n_ = n - SPACE_SCALE_N;
    } else {
        n_ = n;
    }
    m_ = m;
};

void Node::setM(int m) {
    if (m < 0) {
        m_ = m + TIME_SCALE_M;
    } else if (m > TIME_SCALE_M) {
        m_ = m - TIME_SCALE_M;
    } else {
        m_ = m;
    }
}

void Node::setN(int n) {
    if (n < 0) {
        n_ = n + SPACE_SCALE_N;
    } else if (n > SPACE_SCALE_N) {
        n_ = n - SPACE_SCALE_N;
    } else {
        n_ = n;
    }
}

void Node::setNumber(double z){
    number = z;
}

double Node::getNumber() {
    return number;
}

int main() {
    vector<Node> inits(SPACE_SCALE_N); // initial conditions
    double velocity; // velocity = c > 0
    double a; // shift

    for (int i = 0; i < inits.size(); i++) {
        double z = velocity/(2*pow(cosh(sqrt(velocity)*(i - a)/2), 2)); // solitonic initial conditions
        inits.at(i).setM(0);
        inits.at(i).setN(i);
        inits.at(i).setNumber(z);
    }

    for (const Node& node: inits) {

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