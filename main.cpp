#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

class Node
{
private:
    int m; // time
    int n; // space // u{m, n} = u(mt, nh)
public:
    Node(int m, int n); // constructor
    int GetTime() const; // get time-position
    int GetSpace() const; // get space-position
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