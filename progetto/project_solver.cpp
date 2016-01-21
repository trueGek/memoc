#include <stdio.h>
#include "project_solver.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <math.h>
#include <chrono>
#include "Reader.cpp"
#include <float.h>
#include <time.h>


using namespace std;

// Definizione delle variabili globali utili per il calcolo

int for_random = 0;		// Variabile che mi serve per il calcolo dei valori random nei vari metodi

double SA_sol = 0;		// Soluzione corrente del siumlated annealing
double BSA_sol = DBL_MAX;	// Miglior soluzione trovata durante l'elaborazione per uno specifico dataset
double WSA_sol = DBL_MIN;	// Peggior soluzione trovata durante l'elaborazione per uno specifico dataset
double MSA_sol = 0;		// Valore medio di tutte le soluzioni trovate durante l'elaborazione per uno specifico dataset
vector<double> TSA_sol (10);	// Vettore delle soluzioni trovate per uno specifico dataset, utile per il calcolo della varianza
int SA_count = 1;  		// Mi serve per individuare a quale test del Simulated Annealing per uno specifico dataset sono arrivato


/* ---------------------- Metodi comuni ai solver Simulated Annealing e Local Search ---------------------- */

project_solver::project_solver(vector<vector<float> > data){
	
	// Definizione delle variabili che utilizzerò all'interno del solver

	n = data.size();	// Dimensione del dataset
	xy.resize(n);		// Matrice delle distanze - Definizione dimensioni
	xy = data;		// Matrice delle distanze - Definizione elementi

}


float project_solver::evaluate(vector<int> & node){
	
	float dist = 0;		// Variabile per il calcolo delle distanze nella soluzione trovata
	
	for(int i = 0; i < node.size() -1; i++){
		
		dist += xy[node[i]][node[i+1]];		// La distanza totale è uguale alla distanza fino ad ora calcolata più la distanza tra l'elementi in posizione i e i+1
		
	}
	
	return dist;		// Infine ritorno il valore dist calcolato all'interno del ciclo for
	
}


vector<int> project_solver::getInitialSol(bool random){

	vector<int> sol = vector<int>(n+1); 	// Definisco il vettore delle soluzioni come un vettore di lunghezza n+1
	
	
	// Ciclo for per l'inizializzazione del vettore sol

	for (int i = 0; i < n; i++){
		sol[i] = i;
	}

	sol[n+1]=sol[0];			// Definisco infine la posizione n+1 di sol a 0 poiché il percorso deve ritornare al punto di partenza
	
	if(!random){
		return sol;			
	}
	
	srand(time(NULL)+for_random);		// Definisco il seme per il calcolo dei valori casuali successivi
	for_random = for_random +1;
	

	// Ciclo for che andrà a scambiare n*2 volte delle coppie di elementi casuali all'interno del vettore delle soluzioni

	for (int i = 0; i < 2*n; i++){
		
		int pos1 = int( rand()%(n-1) ) +1;// Calcolo la prima posizione 
		int pos2 = int( rand()%(n-1) ) +1;// Calcolo la seconda posizione
		
		// Operazioni per lo scampio di posto degli elementi situati nelle posizioni pos1 e pos2

		int temp = sol[pos1];
		sol[pos1] = sol[pos2];
		sol[pos2] = temp;
		
	}
}


/* ---------------------- Metodi relativi alla Simulated Annealing ---------------------- */

vector<int> project_solver::getNeigh(vector<int> node, int k, bool random){
	

	// if che mi assicura che non ci siano valori di k < 2, nel qual caso modifica il valore proprio a 2

	if ( k < 2 ){
		k = 2;
	}

	vector<int> newV = node;		// Definisce un nuovo vettore delle soluzioni uguale a quello dato in input

	
	// if che, nel caso in cui rand sia true, va ad inizializzarmi il seme per il calcolo dei valori random

	if (random){
		srand(time(NULL)+for_random);
		for_random = for_random +1;
	}else{
		srand(0);
	}
	

	int pos1, pos2;				// Definisco due nuove variabili che mi serviranno per decidere i punti di scambio


	// Ciclo while che opera finché pos1 e pos2 sono uguali. Quando sono diverse termina

	do{
		pos1 = int( rand()%(n-1) ) +1;	// Primo valore attribuito
		pos2 = int( rand()%(n-1) ) +1;	// Secondo valore attribuito

	}while(pos1 == pos2);
		
	
	// if che mi assicura che pos1 sia minore o uguale di pos2

	if(pos1>pos2){
		int temp = pos2;
		pos2 = pos1;
		pos1 = temp;	
	}
		

	// for che mi va ad effettuare l'inversione di tutti gli elementi contenuti nell'intervallo pos1 - pos2

	for (int i = pos1+1, j = pos2-1; i < j ; i++, j--){
		newV[i] = node[j];
		newV[j] = node[i];			
	}


	return newV;				// Ritorno in output il nuovo vettore calcolato dal nostro metodo

}


float project_solver::getSimAnnealing(int multistart, bool random){
	
	sol = getInitialSol(random);		// Inizializzo la soluzione attuale ottenendo quella iniziale

	vector<int> newSol;			// Vettore che mi servirà per mantenere in memoria il calcolo del vicinato eseguito da getNeigh
	
	double step = 1;			// Indice che mi serve per memorizaare il numero di passi eseguiti fino ad ora
	float temp;				// Valore che utilizzo per memorizzare la temperatura
	double n_passi = 100000.0 * n / 5;	// Valore che indica il numero di passi da eseguire, aumentano all'aumentare del numero di nodi del dataset


	// Ciclo while che opera finché non ha eseguito tutti i passi del caso

	while (step < n_passi){
		
		newSol = getNeigh(sol,2,true);			// Calcolo la nuova soluzione tramite la funzione getNeigh
		float de = evaluate(sol) - evaluate(newSol);	// Vado ad effettuare la valutazione della nuova soluzione


		// if che mi servirà per determinare se la nuova soluzione trovata è migliorativa o peggiorativa

		if (de > 0){
			
			// Entrare in questo ramo significa aver trovato una nuova soluzione migliorativa tramite getNeigh, aggiorno quindi sol a newSol

			sol = newSol;
		}else{ 

			// Se invece sono entrato nel ramo else significa che la soluzione è peggiorativa. A questo punto dovrò capire se accettare o meno tale soluzione

			temp = 1-(step/n_passi);	// Calcolo la temperatura come segue
			double prob = exp((de)/temp);	// Calcolo il valore che utilizzerò per valutare se accetto la soluzione peggiorativa
			srand(time(NULL)+for_random);		// Genero il seme per ottenere il numero casuale
			for_random = for_random + 1;

			// if che mi servirà per determinare se accetto la soluzione peggiorativa o meno

			if (prob*100 > (rand()%100) ){

				// Se la probabilità prob calcolata è maggiore del numero random % 100 allora accetterà la soluzione peggiorativa, altrimenti non eseguo l'azione

				sol = newSol;
			}
			
		}
		
		step ++; // Incremento il valore step per procedere nel ciclo
		
	}

	TSA_sol[SA_count-1] = evaluate(sol);	// Inserisco nel vettore delle soluzioni totali la soluzione corrente
	SA_sol = evaluate(sol);			// Mi calcolo la soluzione corrente
	MSA_sol = MSA_sol + SA_sol;		// Mi calcolo la somma di tutte le soluzioni trovate fino ad ora


	// if che mi serve per capire se la soluzione attualmente calcolata è la migliore fino ad ora o meno

	if(SA_sol < BSA_sol) {
		
		BSA_sol = SA_sol;

	}


	// if che mi serve per capire se la soluzione attualmente calcolata è la peggiore ad ora o meno

	if(SA_sol >= WSA_sol) {
		
		WSA_sol = SA_sol;

	}

	
	// Come ultima operazione vado a restituire il valore restituito dalla funzione di localSearch. Tramite questa operazione vado ad effettuare una fase di intensificazione per 
	// migliorare il risultato ottenuto durante la fase di Simulated Annealing

	return localSearch(sol);		// Alla localSearch do in pasto la soluzione trovata attualmente che verrà eventualmente migliorata 
	
}


/* ---------------------- Metodi relativi alla Local Search ---------------------- */

vector<int> project_solver::findBestN(vector<int> sol){

    vector<int> neigh = sol;

    for (int s1 = 1; s1 < n-1; s1++){  
        for (int s2 = s1+1; s2 < n; s2++){

            neigh = sol;
            for (int i = s1, j = s2; i < j; i++, j--){
                float ppp = neigh[i];
                neigh[i] = sol[j];
                neigh[j] = ppp;
            }

            if (evaluate(neigh) < evaluate(sol)){
                return neigh;
            }
        }
    }
    return neigh;
}


float project_solver::localSearch(vector<int> sol){


	// Ciclo while che continuerà ad operare finché non ritornerò la soluzione

	while (true){

		vector<int> newSol = findBestN(sol);	// Calcolo il vicinato migliore tramite la funzione findBestN
		if (evaluate(newSol) >= evaluate(sol)){	// Se newSol è peggiore o uguale a sol allora ritornerò sol
			return evaluate(sol);
		}else{					// Se newSol è migliore di sol allora aggiornerò sol a newSol
			sol = newSol;
		}
	}
}


float project_solver::getLocalSearch(){

	return localSearch(getInitialSol(true));	// ritorno il risultato calcolato dalla localSearch dandole in pasto la soluzione iniziale

}


/* ---------------------- Metodi scartati ---------------------- */

/*
float project_solver::incEvaluate(vector<int> & node, int ub, int lb){

	float dist = 0;

	// partendo da i = 1 tolgo arco 0__1 e fermandomi a node.size()-2 tolgo arco i-2__i-1	

	dist = evaluate(node);

	int s = node.size();

	dist = dist - xy[node[lb-1]][node[lb]] - xy[node[ub]][node[ub+1]] + xy[node[lb-1]][node[ub]] + xy[node[lb]][node[ub+1]];
	
}
*/


/* ---------------------- MAIN ---------------------- */

main(int argc, char* argv[]){	
        
	if (argc != 2){
		
		cout << "usage: ./s_m filename\n";
		exit(0);

	}

	ofstream SaveFile( "results_new_vers.txt", ios_base::app);
		SaveFile << endl << endl << argv[1] << endl << endl;
		SaveFile.close();
	

	
	char* c = argv[1];
	Reader* initializator = new Reader(c);
	
	int test = 10;
	int x = 0;
	int* p_x = &x;

	std::chrono::system_clock::time_point start, end;
    
    	for (int p = 1; p <= initializator->problems; p++){

		SA_sol = 0;
		BSA_sol = DBL_MAX;
		WSA_sol = DBL_MIN;
		MSA_sol = 0;
		TSA_sol.clear();
		//vector<double> TSA_sol (10);  		
		SA_count = 1;

		int n_nodes = initializator->nodes;
		vector<vector<float> > xy;
     	 	initializator->getNextProblem(xy);
      	  	cout << "problema " << p << endl;
        	float total_time = 0.0f;
		cout << "Numero nodi: " << n_nodes << endl;

		double best_sol = DBL_MAX;
		double worst_sol = DBL_MIN;
		double med_sol = 0;
		vector<double> total_sol (10);  	
        	
		for (int tt = 1; tt <= test; tt++){         		

			project_solver* solver = new project_solver(xy);
		        start = std::chrono::system_clock::now();
			float sol = solver->getSimAnnealing(1,true);
			
			//float sol = solver->getLocalSearch();		        
			med_sol = med_sol + sol;
			end = std::chrono::system_clock::now();  
			float cur_time = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()/1000.0f;
			//float cur_time = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();			
			total_time = total_time + cur_time;
			
			SA_count++;
			
			total_sol[tt-1] = sol;

			if (sol <= best_sol){

				best_sol = sol;

			}
			
			if (sol >= worst_sol){

				worst_sol = sol;
			}

			


		}

		// Ritorno in output i vari valori ottenuti eseguendo l'algoritmo

		double Variance_SA = 0;
		double Variance = 0;

		for (int i = 0; i<test; i++){
				
			Variance_SA += (TSA_sol[i]-(MSA_sol/test))*(TSA_sol[i]-(MSA_sol/test));
			Variance += (total_sol[i]-(med_sol/test))*(total_sol[i]-(med_sol/test));

		}


		MSA_sol /= test;
		Variance_SA /= test;
		med_sol /= test;
		total_time /= test;
		Variance /= test;

		cout << endl;

		cout << "Average_Time " << p << ": " << total_time << endl;
		
		cout << endl;

		cout << "Best_Sol_SA " << p << ": " << BSA_sol << endl;
		cout << "Worst_Sol_SA " << p << ": " << WSA_sol << endl;
		cout << "Average_Sol_SA " << p << ": " << MSA_sol<< endl;
		cout << "Variance_SA " << p << ": " << Variance_SA << endl;
		cout << endl;	

		cout << "Best_Sol " << p << ": " << best_sol << endl;
		cout << "Worst_Sol " << p << ": " << worst_sol << endl;
		cout << "Average_Sol " << p << ": " << med_sol << endl;
		cout << "Variance " << p << ": " << Variance << endl;
		cout << endl;
			
		cout << "*************" << endl;
		cout << endl;

		

		ofstream SaveFile( "results_new_vers.txt", ios_base::app);
		SaveFile << p << "\t"  << total_time << "\t" << BSA_sol << "\t" << WSA_sol << "\t" << MSA_sol<< "\t" << Variance_SA << "\t" <<best_sol << "\t" << worst_sol << "\t" << med_sol << "\t" << Variance << endl;
		SaveFile.close();

	}

}
