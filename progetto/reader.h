#include <iostream>
#include <vector>
using namespace std;

class Reader{
public:

    int next;       //prossimo problema da restituire
    int problems;   //numero di problemi nel file
    int nodes;        //dimensione matrice del problema

    Reader(const char* filename = new char());
    ~Reader();
    
    void getNextProblem(vector<vector<float> > & xy);

private:
    void getInner(string text, char pi, char po);
    void split(string text, string delimiter, vector<float> & e);
};
