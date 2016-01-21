#include <iostream>
#include <vector>

using namespace std;

class project_solver{

	public:
		
		// Metodi pubblici che includono il solver e i metodi per ottenere l'una o l'altra soluzione

		project_solver(vector<vector<float> > data);				// Classe principale
		
		float getSimAnnealing(int multistart = 1, bool random = false );	// Calcolo soluzione Simulated Annealing
		float getLocalSearch();							// Calcolo soluzione Local Search

		
	private:
		
		// Metodi e variabili private utili per il calcolo della soluzione

		vector<vector<float> > xy; 			// Matrice delle distanze
		int n; 						// Numero elementi dataset
		vector<int> sol; 				// Vettore della soluzione attuale
		vector<int> getInitialSol(bool random = false); // Funzione che calcola la soluzione iniziale necessaria per i successivi passi del calcolo
		vector<int> getNeigh(vector<int> node, int k = 2, bool random = false); // Funzione che va ad ottenere il vicinato per un determinato vettore delle soluzioni
		float evaluate(vector<int> & node);		// Funzione per la valutazione del vettore delle soluzioni
		float localSearch(vector<int>);			// Funzione per il calcolo della local search. E' separato rispetto a getLocalSearch perché lo utilizzerò anche in getSimAnnealing
		vector<int> findBestN(vector<int>);		// Trova il vicinato migliore

		// float incEvaluate(vector<int> & node, int lb, int ub); // Funzione scartata 

		

};
