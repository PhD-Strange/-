#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

int TIME_SCALE_M = 100;
int SPACE_SCALE_N = 100;
int TIME_SPLIT_T = 0.01;
int SPACE_SPLIT_H = 0.03;

class Node
{
private:
    int m; // u{m, n} = u(mT, nH)
    int n;
public:
    Node(pair<int, int> u); // constructor
};

Node::Node(pair<int, int> u) {
    if (u.second < 0) {
        u.second = u.second + SPACE_SCALE_N;
    } else if (u.second > SPACE_SCALE_N * SPACE_SPLIT_H) {
        u.second = u.second - SPACE_SCALE_N;
    }
};

int main() {
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