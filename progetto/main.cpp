#include <cstdio>
#include <iostream>
#include <vector>
#include "cpxmacro.h"
#include "Reader.cpp"
#include <chrono>

using namespace std;

int status;
char errmsg[BUF_SIZE];
	
const int NAME_SIZE = 512;
char name[NAME_SIZE];
	
int n;						// Dichiarazione del valore n che indica il numero totale di elementi nel mio vettore xy
vector<vector<float> > xy;			// Dichiarazione del vettore xy
	
void setupLP(CEnv env, Prob lp)
{
		
	int cVar = 0;				// Dichiarazione dell'indice che mi servirà per individuare in che posizione si trovano le variabili nel vettore di output
	vector<vector<int> > xM;		// Matrice che memorizzerà le posizioni delle variabili x
	vector<vector<int> > yM;		// Matrice che memorizzerà le posizioni delle variabili y

	
	// Aggiunta delle variabili x

	{
		xM.resize(n);			

		double obj = 0;			
		double lowerB = 0.0;
		double upperB = CPX_INFBOUND;

		for (int i = 0; i < n; i++)
		{
			xM[i].resize(n);

			for (int j = 0; j < n; j++)
			{
				char xtype = 'I';	
				snprintf(name, NAME_SIZE, "x_%c_%c", i,j);
				char* xname = (char*)(&name[0]);
				CHECKED_CPX_CALL( CPXnewcols, env, lp, 1, &obj, &lowerB, &upperB, &xtype, &xname );
				xM[i][j]=cVar;
				cVar++;
			}
		}
	}
	

	// Aggiunta delle variabili y

	{
		yM.resize(n);

		double lowerB = 0;
		double upperB = 1;

		for (int i = 0; i < n; i++)
		{
			yM[i].resize(n);

			for (int j = 0; j< n; j++)
			{
				char ytype = 'I';	
				double obj = xy[i][j];
				snprintf(name, NAME_SIZE, "x_%c_%c", i,j);	
				char* yname = (char*)(&name[0]);
				CHECKED_CPX_CALL( CPXnewcols, env, lp, 1, &obj, &lowerB, &upperB, &ytype, &yname );
				yM[i][j]=cVar;
				cVar++;
			}
		} 
	}   


	// 1) Definizione del primo vincolo -- Flusso uscente da x0j deve essere massimo, cioè |N|

	{
		char sense = 'E';
		int matbeg = 0;
		double fluxVal = n;
		vector<int> x0j(n);
		vector<double> c0j(n,1.0);
		

		// Tramite queste operazioni faccio in modo che il flusso uscente da x0 diretto verso il nodo successivo del percorso xj sia pari ad n
		// Deve esserci chiaramente un solo nodo xj verso cui il flusso viene spostato

		for (int j = 0; j < n; j++)	
		{
			x0j[j] = xM[0][j];	// In x0j mi troverò verso quale nodo j vi è un arco che parte da x0
		}

		CHECKED_CPX_CALL( CPXaddrows, env, lp, 0, 1, x0j.size(), &fluxVal, &sense, &matbeg, &x0j[0], &c0j[0], NULL, NULL );
	}

	
	// 2) Definizione del secondo vincolo -- Ogni nodo spreca al massimo una unità di flusso, tranne il nodo di partenza

	{
		char sense = 'E';
		int matbeg = 0;
		double fluxUse = 1;


		// Tramite queste operazioni faccio in modo che ogni nodo utilizzi al massimo una unità di flusso, non può utilizzarne di più
		// Percorrere un arco da k a i ha costo -1, il percorso contrario cioè da i a k ha chiaramente costo 1 

		for (int k=1; k<n; k++) 
		{
			vector<int> id(2*n);		// id ha dimensione 2*n perché conterrà al suo interno sia gli archi k--j che gli archi i--k
			vector<double> coeff(2*n);	// coeff ha dimensione 2*n perché anch'esso conterrò i coefficienti di "uso" del flusso

			for (int i=0; i<n; i++)
			{
				id[i] = xM[i][k];	// In id[0,..,n] troverò al suo interno gli archi che vi sono tra i e k
				coeff[i] = 1;		// Il costo di un arco da i a k, qual'ora vi sia, è pari a 1 (starei tornando indietro)
			}

			for (int j=0; j<n; j++)
			{
				id[n+j] = xM[k][j];	// In id[n+1,..,2*n] troverò al suo interno gli archi che vi sono tra k e j
				coeff[n+j] = -1;	// Il costo di un arco da k a j, qual'ora vi sia, è pari a -1 (starei proseguendo in avanti)
			}

			CHECKED_CPX_CALL( CPXaddrows, env, lp, 0, 1, id.size(), &fluxUse, &sense, &matbeg, &id[0], &coeff[0], NULL, NULL );
		}
	}
	

	// 3) Definizione del terzo vincolo -- Ogni nodo ha un solo arco in entrata

	{
		char sense = 'E';
		int matbeg = 0;
		double exitArc = 1;


		// Tramite queste operazioni faccio in modo che ogni ogni nodo abbia in entrata verso di sé uno e un solo arco proveniente da un altro nodo

		for (int i=0; i<n; i++)
		{
			vector<int> idy(n);
			vector<double> coeff(n,1);

			for (int j=0; j<n; j++)
			{
				idy[j] = yM[i][j];	// In idy[0,..,n] troverò al suo interno gli archi entranti che vi sono tra gli n nodi
			}

			CHECKED_CPX_CALL( CPXaddrows, env, lp, 0, 1, idy.size(), &exitArc, &sense, &matbeg, &idy[0], &coeff[0], NULL, NULL );
		}
	}
	
	
	// 4) Definizione quarto vincolo -- Ogni nodo ha un solo arco in uscita

	{
		char sense = 'E';
		int matbeg = 0;
		double enterArc = 1;


		// Tramite queste operazioni faccio in modo che ogni ogni nodo abbia in uscita verso un altro nodo uno e un solo arco proveniente

		for (int j=0; j<n; j++)
		{
			vector<int> idy(n);
			vector<double> coeff(n,1);

			for (int i=0; i<n; i++)
			{
				idy[i] = yM[i][j];	// In idy[0,..,n] troverò al suo interno gli archi uscenti che vi sono tra gli n nodi
			}

			CHECKED_CPX_CALL( CPXaddrows, env, lp, 0, 1, idy.size(), &enterArc, &sense, &matbeg, &idy[0], &coeff[0], NULL, NULL );
		}
	}

		
	// 5) Definizione quinto vincolo -- Se vi è un'unità di flusso trasportata da i a j deve di conseguenza esserci un arco che va da i a j
	{
		char sense = 'L';
		int matbeg = 0;
		double rhs = 0;

		for (int i=0; i<n; i++)
		{
			for (int j=0; j<n; j++)
			{
				vector<int> id(2);
				vector<double> coeff(2);

				id[0] = xM[i][j];	// in id[0] troverò se è stata utilizzata un'unità di flusso dal nodo i al nodo j
				id[1] = yM[i][j];	// in id[1] troverò se vi è un arco dal nodo i al nodo j
				coeff[0] = 1;
				coeff[1] = -n;

				CHECKED_CPX_CALL( CPXaddrows, env, lp, 0, 1, 2, &rhs, &sense, &matbeg, &id[0], &coeff[0], NULL, NULL );
			}
		}
	}

	// print (debug)
	CHECKED_CPX_CALL( CPXwriteprob, env, lp, "pr.lp", NULL );
	
}

int main (int argc, char *argv[])
{
	
    if (argc != 2)
    {
	cout << "usage: ./main filename\n";
	exit(0);
    }


    // Vado ad inserire nel file di output le informazioni relative a quale dataset sto analizzando

    ofstream SaveFile("results_TEST.txt", ios_base::app);
    SaveFile << endl << endl << argv[1] << endl << endl;
    SaveFile.close();


    // char cc[40] = "istanze/100_rand_dataset";	// Comando per dare in input uno specifico dataset

    char* c = argv[1];					// Prendo il nome del file in input
    Reader* dataReader = new Reader(c);			// Do in pasto al reader il nuovo file in input
    
    int test = 10;					// Numero di test da effettuare per ogni dataset

    std::chrono::system_clock::time_point start, end;	// Mi definisco lo strumento per rilevare il tempo
    

    // Cicli for principali, mi servono per andare a calcolare le varie soluzioni che otterrò dai vari dataset utilizzando le funzioni di CPLEX

    for (int curP = 1; curP <= dataReader->problems; curP++)	// Questo ciclo for serve per scorrere tutti i problemi per ogni tipologia di dataset
    {
        dataReader->getNextProblem(xy);			// Valorizzo il vettore xy con i valori del nuovo problema ottenuto da dataReader
        float total_time = 0.0f;			// Definizione del tempo totale impiegato per la risoluzione di un dataset: somma tutti i tempi di esecuzione. All'inizio è 0
	double timelimit = 1200;			// Definizione del tempo massimo che la funzione CPLEX può ciclare per trovare una soluzione
	double cur_sol;					// Dichiarazione di una variabile che mi servirà per tenere traccia di qual'è la soluzione del dataset corrente
	n = dataReader->nodes;				// Definisco il numero di nodi del problema corrente ottenendolo dal dataReader
	

	// Stampe per la visualizzazione dei risultati a console

	cout << "problema " << curP << endl;	
	cout << "Numero nodi: " << n << endl;
		
        for (int curT = 1; curT <= test; curT++)	// Questo secondo ciclo for serve per eseguire tutti i test previsti per ogni problema di ogni dataset
	{
        
           
	        try
	        {
		        DECL_ENV( env );
		        DECL_PROB( env, lp );		
		        setupLP(env, lp);

		        CPXsetdblparam(env, CPX_PARAM_TILIM, timelimit);	// Sto dando al CPLEX un limite di tempo massimo per trovare la soluzione del problema corrente. Superato tale limite restituirà quel che ha trovato
			
		        start = std::chrono::system_clock::now();		// Inizio a rilevare il tempo

		        CHECKED_CPX_CALL( CPXmipopt, env, lp );			// Richiamo la funzione di CPLEX per cui trovare l'ottimo del problema corrente

		        end = std::chrono::system_clock::now();			// Interrompo la rilevazione del tempo

		        double sol;						// Dichiaro variabile che definirò poi con l'ottimo trovato da CPLEX
		        CHECKED_CPX_CALL( CPXgetobjval, env, lp, &sol );	// Estraggo l'ottimo trovato da CPLEX

		        float cur_time = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()/1000.0f;	// Calcolo del tempo corrente
		        total_time += cur_time;					// Calcolo del tempo totale
				
			cur_sol = sol;						// Imposto la soluzione corrente a sol
				
		        CHECKED_CPX_CALL( CPXwriteprob, env, lp, "log_file.lp", NULL);	// Mi salvo il file di log
		
		        CPXfreeprob(env, &lp);					
		        CPXcloseCPLEX(&env);					 
		
	        }
	        catch(std::exception& e)
	        {
		        // cout << ">>>EXCEPTION: " << e.what() << endl;
	        }
	

	  }


	  // Stampa a console dei valori trovati

	  cout << "Av_Time " << curP << ": " << total_time / test << endl;
	  cout << "Best_Val " << curP << ": " << cur_sol << endl;
          cout << endl;


	  // Salvataggio su results_cplex.txt dei risultati del solver

	  ofstream SaveFile("results_TEST.txt", ios_base::app);
          SaveFile << curP << "\t" << total_time/test << "\t" << cur_sol << endl;
          SaveFile.close();

	}
	return 0;
}
