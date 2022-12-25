#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <sstream>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <random>

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

vector<string> randomSample(vector<string> structure, int length)
{
    // select unique elements from given vector
    srand (time(0));
    vector<string> elements;
    int i = 0;
    while (i < length)
    {
        int random = rand() % structure.size(); 
        string newElement = structure[random]; 
        if (!(count(elements.begin(), elements.end(), newElement)))
        {   
            elements.push_back(newElement);
            ++i;
        }
    }
    return elements;
}

vector<string> splitBySpace(string line)
{
    vector<string> results;
    istringstream stream(line);
    string newString;
    while (getline(stream, newString, ' '))
    {
        if (!newString.empty())
        {
            results.push_back(newString);
        }
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
    // amino acids sequence checking
    {
        if (!(count(aaCodes.begin(), aaCodes.end(), i)))
        {
            correct = false;
        }
    }
    for (char j : ss) 
    // secondary structure checking
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
    lines.erase(lines.begin()); // remove headlines 
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
    // save given vector line by line
    ofstream linesFile;
    linesFile.open(file);
    for (string line : lines) 
    {
        linesFile << line << endl;
    }
    linesFile.close();
}

int main(int argc, char* argv[])
{
    // read arguments from console
    string directory = argv[1];
    string output = argv[2];
    int length = stoi(argv[3]);
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
    vector<string> selectedLines = randomSample(allLines, length);
    saveToFile(selectedLines, output);
    return 0;
}
