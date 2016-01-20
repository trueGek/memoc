#include <stdio.h>
#include "reader.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <string>

using namespace std;

ifstream XX;

Reader::~Reader() {XX.close();}
Reader::Reader(const char* filename) {
    XX.open(filename);
    if (!XX) return;
    string header;
    //leggo il numero di nodi e problemi nell'intestazione
    getline(XX,header,'\n');
    header = header.substr(1,header.length()-2);
    vector<float> param(2);
    split(header,", ",param);
    //cout << "#nodi:\t" << param[0] << "\t#problemi:\t" << param[1] << endl;
    Reader::next = 0;
    Reader::problems = param[0];
    Reader::nodes = param[1];
    //XX.close();
    
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
    getline(XX,problem,'\n');
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

/*
main(){
    char cc[20] = "data";
    char* c = cc;
    Reader* reader = new Reader(c);
    //delete c;     //  wut??
    vector<vector<float> > xy;
    
    for (int p = 1; p <= reader->problems; p++){
        reader->getNextProblem(xy);
        cout << endl 
             << "  -----------" << endl
             << "  Problem n:" << p << endl
             << "  -----------" << endl;
        for (int row = 0; row < reader->nodes; row++){
            for (int col = 0; col < reader->nodes; col++){
                cout << xy[row][col] << "\t" ;
            }
            cout << endl;
        }
    }
    cout << endl;
    delete reader;
}
*/
