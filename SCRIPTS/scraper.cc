#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <cstring>

using namespace std;

/*
    script for data generating
    extract random fragments from proteins collection
    example line shown below
    1a62A 3 THR 4 MNLTELKNTPV CCHHHHHCCCH 1.284 -5.477 -2.600 99.739 -124.579 91.327 49.439 90.431 53.359 90.948 45.020 87.256 56.852
*/

vector<string> splitBySpace(string line)
{
    vector<string> results;
    istringstream stream(line);
    string newString;
    while (getline(stream, newString, ' '))
    {
        results.push_back(newString);
    }
    return results;
}

vector<char> stringToVector(string str)
{
    vector<char> elements;
    int n = str.length();
    char elementsArray[n + 1];
    strcpy(elementsArray, str.c_str());
    for (int i = 0; i < n; i++)
    {
        elements.push_back(elementsArray[i]);
    }
    return elements;
}

bool isCorrect(string line)
{
    // check if given line contains sequence and secondary structure written in a proper way
    bool correct = true;
    vector<string> elements = splitBySpace(line);
    int AA_ORIDNAL = 4;
    int SS_ORDINAL = 5;
    string aa = elements[AA_ORIDNAL];
    string ss = elements[SS_ORDINAL];
    int aaLength = aa.length();
    int ssLength = ss.length();
    return correct;
}

int main()
{
    return 0;
}
