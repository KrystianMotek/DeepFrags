#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <sstream>
#include <cstring>
#include <fstream>
#include <algorithm>

using namespace std;
namespace fs = filesystem;

/*
    script for data generating
    extract random fragments from proteins collection
    example line shown below
    1a62A 3 THR 4 MNLTELKNTPV CCHHHHHCCCH 1.284 -5.477 -2.600 99.739 -124.579 91.327 49.439 90.431 53.359 90.948 45.020 87.256 56.852
*/

vector<char> aaCodes = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
vector<char> ssCodes = {'H', 'E', 'C'};

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
    vector<char> results;
    int n = str.length();
    char elementsArray[n + 1];
    strcpy(elementsArray, str.c_str());
    for (int i = 0; i < n; i++)
    {
        results.push_back(elementsArray[i]);
    }
    return results;
}

bool isCorrect(string line)
{
    // check if given line contains sequence and secondary structure written in a proper way
    bool correct = true;
    vector<string> items = splitBySpace(line);
    int AA_ORDINAL = 4;
    int SS_ORDINAL = 5;
    vector<char> aa = stringToVector(items[AA_ORDINAL]);
    vector<char> ss = stringToVector(items[SS_ORDINAL]);
    for (char i : aa)
    {
        if (!(count(aaCodes.begin(), aaCodes.end(), i)))
        {
            correct = false;
        }
    }
    for (char j : ss) 
    {
        if (!(count(ssCodes.begin(), ssCodes.end(), j)))
        {
            correct = false;
        }
    }
    return correct;
}

vector<string> readLineByLine(string file) 
{
    vector<string> lines;
    ifstream linesFile;
    string line;
    linesFile.open(file);
    while (getline(linesFile, line))
    {
        lines.push_back(line);
    }
    linesFile.close();
    lines.erase(lines.begin());
    return lines;
}

vector<string> getDatFiles(string directory) 
{
    vector<string> files;
    string const extension = ".dat";
    for (const auto & entry : fs::directory_iterator(directory))
    {
        if (entry.path().extension() == extension)
            files.push_back(entry.path().string());
    }
    return files;
}

void saveToFile(vector<string> lines, string file) 
{   
    ofstream linesFile;
    linesFile.open(file);
    for (string line : lines) 
    {
        cout << line << endl;
    }
}

int main(int argc, char* argv[])
{
    string directory = argv[1];
    vector<string> files = getDatFiles(directory);
    vector<string> allLines;
    for (string file : files)
    {
        vector<string> lines = readLineByLine(file);
        for (string line : lines)
        {
            if (isCorrect(line))
            {
                allLines.push_back(line);
            }
        }
        vector<string>().swap(lines);
    }
    vector<string>().swap(files);
    return 0;
}
