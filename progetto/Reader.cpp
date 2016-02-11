#include <stdio.h>
#include "reader.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <string>

using namespace std;

ifstream inputFile;

Reader::~Reader()
{
	inputFile.close();
}


Reader::Reader(const char* filename) 	// Funzione per la lettura di un file
{
    inputFile.open(filename);
    if (!inputFile) return;
    string header;
	
	
	// Lettura delle info contenute nel file
    
	getline(inputFile,header,'\n');		
    header = header.substr(1,header.length()-2);
    vector<float> param(2);
    split(header,", ",param);

    Reader::next = 0;
    Reader::problems = param[0];		// Va ad impostare quanti problemi contiene il file
    Reader::nodes = param[1];			// Va ad impostare quanti nodi contiene ogni problema nel file
    
}

void Reader::split(string text, string delimiter, vector<float> & e){
    size_t pos = 0;
    int n = 0;
 
    string token;
    while ((pos = text.find(delimiter)) != std::string::npos) {
        token = text.substr(0, pos);
        float v = atof (token.c_str());
        e[n] = v;
        n++;
        text.erase(0, pos + delimiter.length());
    }
    e[n] = atof(text.c_str());
}

void Reader::getNextProblem(vector<vector<float> > & xy){
    xy.resize(nodes);
    string problem;
    getline(inputFile,problem,'\n');
    problem = problem.substr(1,problem.length()-2);
    //cout << problem << endl << endl;
    for (int i = 0; i < nodes; i++){  //per ogni riga
        xy[i].resize(nodes);
        int pos1 = problem.find("[");
        int pos2 = problem.find("]");
        string row = problem.substr(pos1+1,pos2-pos1-1);
        problem.erase(0,pos2+1);
        //cout << s << endl;
        split(row,", ", xy[i]);
    }
}

