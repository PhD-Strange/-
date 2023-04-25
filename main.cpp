#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main() {
    ofstream out; // for recording
    out.open("hello.txt");
    if (out.is_open())
        out << "Hello World?" << endl;
    out.close();
    cout << "File has been written" << endl;

    return 0;
};