#include <iostream>
#include <vector>

using namespace std;

class S_Annealing{

	public:
		
		S_Annealing(vector<vector<float> > data);
		
		// Restituise la soluzione
		float getSolSA(int multistart = 1, bool random = false );
		float getSolLS();

		// float getLocalSearch(int multistart = 1, bool random = false);
		
		// multistart: numero di riavii relativamente al calcolo dlle soluzioni, restituisce quella migliore
		// random è impostata a false di default, ma la modificherò a true se ho che multistart è > di 1
		
	private:
		
		vector<vector<float> > xy; // Matrice delle distanze
		int n;
		vector<int> sol; // Soluzione attuale
		float f; // Funzione obiettivo
		vector<int> findBestNSA(vector<int>);
		vector<int> getInitialSol(bool random = false); // Soluzione iniziale che calcola
		vector<int> getNeigh(vector<int> node, int k = 2, bool random = false); // 
		float evaluate(vector<int> & node);
		float incEvaluate(vector<int> & node, int lb, int ub);
		float localSearch(vector<int>);
		vector<int> findBestN(vector<int>);

		

};
