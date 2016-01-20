#include <stdio.h>
#include "S_Annealing.h"
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

int for_random = 0;

double SA_sol = 0;
double BSA_sol = DBL_MAX;
double WSA_sol = DBL_MIN;
double MSA_sol = 0;
vector<double> TSA_sol (10);
int SA_count = 1;  

S_Annealing::S_Annealing(vector<vector<float> > data){
	
	n = data.size();
	xy.resize(n);
	xy = data;

}

float S_Annealing::getSolSA(int multistart, bool random){
	
	sol = getInitialSol(random);

	vector<int> next;
	

	double step = 1;
	float T;
	double n_passi = 1000000.0;
	
	//cout << "Start_Sol: " << evaluate(sol) << endl;	

	while (step < n_passi){
		

		//next = findBestNSA(sol);
		next = getNeigh(sol,2,true);
		float de = evaluate(sol) - evaluate(next);
		
		// delta è maggiore di 0 significa che è migliore next e quindi aggiorno

		if (de > 0){
			sol = next;
		}else{ // entro in questo ramo solo se de è minore di 0

			//T = 1-(step 		

			T = 1-(step/n_passi);
			double p = exp((de)/T);
			srand(time(NULL));
			if (p*100 > (rand()%100) ){ // mi serve per dire se accetto la mossa peggiorativa secondo la probabilità p casuale
				sol = next;
			}
			
		}
		
		step ++;
		
	}

	

	TSA_sol[SA_count-1] = evaluate(sol);
	SA_sol = evaluate(sol);
	MSA_sol = MSA_sol + SA_sol;

	if(SA_sol < BSA_sol) {
		
		BSA_sol = SA_sol;

	}

	if(SA_sol >= WSA_sol) {
		
		WSA_sol = SA_sol;

	}

	
	return localSearch(sol);
	
}


float S_Annealing::getSolLS(){

	return localSearch(getInitialSol(random));

}




float S_Annealing::localSearch(vector<int> sol){


	while (true){

		vector<int> b_n = findBestN(sol);
		if (evaluate(b_n) >= evaluate(sol)){	
			return evaluate(sol);
		}else{
			sol = b_n;
		}
	}

}


vector<int> S_Annealing::findBestN(vector<int> sol){

    //cout << "Allerta DUE!!!" << endl; 

    int e = n+1;
    vector<int> neigh = sol;

    for (int s1 = 1; s1 < e-2; s1++){  //non serve che scambio il primo e l'ultimo elemento
        for (int s2 = s1+1; s2 < e-1; s2++){
            //valuto il vicino scambiando s1s2
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


vector<int> S_Annealing::findBestNSA(vector<int> sol){

    //cout << "Allerta!!!" << endl; 

    int e = n+1;
    vector<int> neigh = sol;
    //cout << "prima della modifica\t";
    //printSol(neighbor);
    int ntest = 0;
    for (int s1 = 1; s1 < e-2; s1++){  //non serve che scambio il primo e l'ultimo elemento
        for (int s2 = s1+1; s2 < e-1; s2++){
            //valuto il vicino scambiando s1s2
            neigh = sol;
            for (int i = s1, j = s2; i < j; i++, j--){
                float ppp = neigh[i];
                neigh[i] = sol[j];
                neigh[j] = ppp;
            }
            ntest ++;
            //cout << s1 << ' ' << s2 << ' ' << evaluate(neighbor)<< endl;
            //if (s1 ==2) return neighbor;
            //printSol(neighbor);
            if (evaluate(neigh) < evaluate(sol)){
                //bestNeighbor = neighbor;
                //cout << "vicini analizzati: " << ntest << endl;
                return neigh;
            }
        }
    }
    return getNeigh(neigh, 2, true);
}

vector<int> S_Annealing::getInitialSol(bool random){

	int n2 = n;
	vector<int> sol = vector<int>(n2+1);
	
	for (int i = 0; i < n2; i++){
		sol[i] = i;
	}
	sol[n2+1]=sol[0];
	
	if(!random){
		return sol;
	}
	
	srand(time(NULL)+for_random);

	for_random = for_random +1;
	
	for (int i = 0; i < 2*n2; i++){
		
		int s1 = int( rand()%(n2-1) ) +1;
		int s2 = int( rand()%(n2-1) ) +1;
	
//		cout << s1 << " " << s2 << endl;
		
		int temp = sol[s1];
		
		sol[s1] = sol[s2];
		sol[s2] = temp;
		
	}

	
}

vector<int> S_Annealing::getNeigh(vector<int> node, int k, bool random){
	

	if ( k < 2 ){
		k = 2;
	}
	
	vector<int> newV = node;

	if (random){
		srand(time(NULL)+for_random);
		for_random = for_random +1;
		//srand(time(NULL));
	}else{
		srand(0);
	}
	

		int s1, s2;
		do{
			s1 = int( rand()%(n-1) ) +1;
			s2 = int( rand()%(n-1) ) +1;

		}while(s1 == s2);
		
		// Se l'indice s1 è maggiore dell'indice s2 vado a scambiarli in modo che siano in ordine
	
		if(s1>s2){
			int temp = s2;
			s2 = s1;
			s1 = temp;	
		}
		
		for (int i = s1+1, j = s2-1; i < j ; i++, j--){
			newV[i] = node[j];
			newV[j] = node[i];			
		}


	return newV;
	
// 2-opt
	
	
}


float S_Annealing::evaluate(vector<int> & node){
	
	float dist = 0;
	
	for(int i = 0; i < node.size() -1; i++){
		
		dist += xy[node[i]][node[i+1]];
		
	}
	
	return dist;
	
}


/*
float S_Annealing::incEvaluate(vector<int> & node, int ub, int lb){

	float dist = 0;

	// partendo da i = 1 tolgo arco 0__1 e fermandomi a node.size()-2 tolgo arco i-2__i-1	

	dist = evaluate(node);

	int s = node.size();

	dist = dist - xy[node[lb-1]][node[lb]] - xy[node[ub]][node[ub+1]] + xy[node[lb-1]][node[ub]] + xy[node[lb]][node[ub+1]];
	
}
*/


main(int argc, char* argv[]){	
        
	if (argc != 2){
		
		cout << "usage: ./s_m filename\n";
		exit(0);

	}

	ofstream SaveFile( "results.txt", ios_base::app);
		SaveFile << endl << endl << argv[1] << endl << endl;
		SaveFile.close();
	

	
	char* c = argv[1];
	Reader* initializator = new Reader(c);
	
	int test = 10;
	
	int x = 0;

	int* p_x = &x;

	// float f = solver->getSolution(1,true);
	
	// cout << "Risultato ==> " << f << endl;

	std::chrono::system_clock::time_point start, end;
    
	// cout << "Risultati parziali: " << endl;

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

			S_Annealing* solver = new S_Annealing(xy);
		        start = std::chrono::system_clock::now();
			float sol = solver->getSolSA(1,true);
			
			//float sol = solver->getSolLS();		        
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

		

		ofstream SaveFile( "results.txt", ios_base::app);
		SaveFile << p << "\t"  << total_time << "\t" << BSA_sol << "\t" << WSA_sol << "\t" << MSA_sol<< "\t" << Variance_SA << "\t" <<best_sol << "\t" << worst_sol << "\t" << med_sol << "\t" << Variance << endl;
		SaveFile.close();

	}

}
