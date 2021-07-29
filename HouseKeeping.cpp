// HouseKeeping.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Parsing.h"

using namespace std;

int main()
{
    cout << "Hello World!\n";
    cout << "write the number of k for kmers: " << endl;
    int in;
    cin >> in;
    Parsing p = Parsing(in);
    cout << "write the filename in FASTA format: " << endl;
    string filename;
    cin >> filename;

    p.Aread(filename);
    cout<< "write the filename in FASTA format : " << endl;
    cin >> filename;
    p.Qread(filename);





}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
