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
	
	// Dichiarazione delle variabili che utilizzerò all'interno del solver

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
	

	int pos1, pos2;				// Dichiaro due nuove variabili che mi serviranno per decidere i punti di scambio


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

	vector<int> newSol = sol;			// Definisco newSol impostandola uguale a sol


	// Cicli for che mi serviranno per calcolare il vicinato

	for (int pos1 = 1; pos1 < n-1; pos1++){  
		for (int pos2 = pos1+1; pos2 < n; pos2++){

		        newSol = sol;			// Ad ogni nuovo ciclo for ridefinisco newSol a sol


			// Queso ciclo for serve per invertire gli elementi contenuti all'interno dell'intervallo pos1 - pos2

        		for (int i = pos1, j = pos2; i < j; i++, j--){
        	        	float ppp = newSol[i];
        	        	newSol[i] = sol[j];
        	        	newSol[j] = ppp;
        	    	}


			// if che mi va a valutare se la nuova soluzione newSol è migliore della vecchia soluzione sol

        	   	if (evaluate(newSol) < evaluate(sol)){

				// Se entro nel ciclo significa che la nuova soluzione calcolata è migliore e terminerò così restituendo newSol, altrimenti continuo a ciclare

        	        	return newSol;
        	    	}
        	}
    	}

	return newSol;
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


	// Vado ad inserire nel file di output a che dataset si sta riferendo

	ofstream SaveFile( "results_new_vers.txt", ios_base::app);
	SaveFile << endl << endl << argv[1] << endl << endl;
	SaveFile.close();	

	char* c = argv[1];					// Prendo il nome del file in input
	Reader* dataReader = new Reader(c);			// Do in pasto al reader il nuovo file in input
	
	int test = 10;						// Numero di test da effettuare per ogni dataset

	std::chrono::system_clock::time_point start, end;	// Mi definisco lo strumento per rilevare il tempo
    	

	// Cicli for principali, mi serviranno per andare a calcolare le varie soluzioni che otterrò dai vari dataset

    	for (int curP = 1; curP <= dataReader->problems; curP++){	// Il primo ciclo for opererà per curP volte con curP <= al numero di problemi per ogni dataset

		SA_sol = 0;					// Azzeramento della soluzione corrente del Simulated Annealing
		BSA_sol = DBL_MAX;				// Azzeramento della migliore soluzione del Simulated Annealing
		WSA_sol = DBL_MIN;				// Azzeramento della peggiore soluzione del Simulated Annealing
		MSA_sol = 0;					// Azzeramento della media delle soluzioni del Simulated Annealing
		TSA_sol.clear();				// Azzeramento del vettore delle soluzioni ottenute dal Simulated Annealing
		SA_count = 1;					// Azzeramento della del contatore che utilizzerò nel Simulated Annealing

		int n_nodes = dataReader->nodes;		// Definisco il numero di nodi del problema corrente ottenendolo dal dataReader
		vector<vector<float> > xy;			// Dichiaro il vettore xy
     	 	dataReader->getNextProblem(xy);			// Valorizzo il vettore xy con i valori del nuovo problema ottenuto dall dataReader
		float total_time = 0.0f;      	  		// Dichiarazione e valorizzazione del tempo totale impiegato per la risoluzione di un dataset: somma di tutti i tempi

		double best_sol = DBL_MAX;			// Definizione di best_sol ovvero la migliore soluzione globale trovata fino ad un dato momento 
		double worst_sol = DBL_MIN;			// Definizione di worst_sol ovvero la peggiore soluzione globale trovata fino ad un dato momento
		double med_sol = 0;				// Definizione di med_sol ovvero la media delle soluzioni trovate fino a questo momento
		vector<double> total_sol (10);  		// Dichiarazione del vettore di tutte le soluzioni trovate fino ad ora

		// Stampe per la visualizzazione dei risultati a console

		cout << "problema " << curP << endl;		
		cout << "Numero nodi: " << n_nodes << endl;
        	
		for (int curT = 1; curT <= test; curT++){			// Il secondo ciclo for opererà per curT volte con curT <= al numero di test da effettuare per ogni dataset

			project_solver* solver = new project_solver(xy);	// Definisco un nuovo solver su cui poi opererò per ottenere i risultati dei due metodi implementati	
		        
			start = std::chrono::system_clock::now();		// Inizio a rilevare il tempo

			float sol = solver->getSimAnnealing(1,true);		// Richiamo il metodo per avere la soluzione del Simulated Annealing			
			//float sol = solver->getLocalSearch();			// Richiamo il metodo per avere la soluzione della Local Search
	        
			med_sol = med_sol + sol;				// Aggiorno la somma delle soluzioni

			end = std::chrono::system_clock::now();  		// Termino il rilevamento del tempo

			float cur_time = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()/1000.0f;	// Calcolo del tempo corrente
			total_time = total_time + cur_time;			// Calcolo del tempo totale
			
			SA_count++;						// Incremento il numero delle iterazioni totali da effettuare per un dataset
			total_sol[curT-1] = sol;					// Inserisco nel vettore total_sol la soluzione corrente


			// Vado a calcolare la soluzione migliore trovata fino ad ora

			if (sol <= best_sol){	

				best_sol = sol;

			}
			

			// Vado a calcolare la soluzione peggiore trovata fino ad ora

			if (sol >= worst_sol){

				worst_sol = sol;
			}

			


		}


		double Variance_SA = 0;						// Definizione della varianza del Simulated Annealing
		double Variance = 0;						// Definizione della varianza generale


		// Ciclo for per il calcolo della varianza

		for (int i = 0; i<test; i++){
				
			Variance_SA += (TSA_sol[i]-(MSA_sol/test))*(TSA_sol[i]-(MSA_sol/test));
			Variance += (total_sol[i]-(med_sol/test))*(total_sol[i]-(med_sol/test));

		}


		// Aggiornamento dei vari valori necessari per la raccolta dati

		MSA_sol /= test;
		Variance_SA /= test;
		med_sol /= test;
		total_time /= test;
		Variance /= test;


		// Operazioni di output a console per visualizzare i dati

		cout << endl;
		cout << "Average_Time " << curP << ": " << total_time << endl;
		
		cout << endl;

		cout << "Best_Sol_SA " << curP << ": " << BSA_sol << endl;
		cout << "Worst_Sol_SA " << curP << ": " << WSA_sol << endl;
		cout << "Average_Sol_SA " << curP << ": " << MSA_sol<< endl;
		cout << "Variance_SA " << curP << ": " << Variance_SA << endl;
		cout << endl;	

		cout << "Best_Sol " << curP << ": " << best_sol << endl;
		cout << "Worst_Sol " << curP << ": " << worst_sol << endl;
		cout << "Average_Sol " << curP << ": " << med_sol << endl;
		cout << "Variance " << curP << ": " << Variance << endl;
		cout << endl;
			
		cout << "*************" << endl;
		cout << endl;

		
		// Operazioni per la scrittura su file dei dati calcolati

		ofstream SaveFile( "results_new_vers.txt", ios_base::app);
		SaveFile << curP << "\t"  << total_time << "\t" << BSA_sol << "\t" << WSA_sol << "\t" << MSA_sol<< "\t" << Variance_SA << "\t" <<best_sol << "\t" << worst_sol << "\t" << med_sol << "\t" << Variance << endl;
		SaveFile.close();

	}

}
