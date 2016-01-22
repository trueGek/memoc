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
	
    if (argc != 2){
	
	cout << "usage: ./main filename\n";
	exit(0);

    }

    ofstream SaveFile("results_TEST.txt", ios_base::app);
    SaveFile << endl << endl << argv[1] << endl << endl;
    SaveFile.close();

    // char cc[40] = "istanze/100_rand_dataset";
    char* c = argv[1];
    Reader* initializator = new Reader(c);
    
    int test = 10;

    std::chrono::system_clock::time_point start, end;
    
    for (int p = 1; p <= initializator->problems; p++){
        initializator->getNextProblem(xy);
        cout << "problema " << p << endl;
        
        float tc = 0.0f;

	double timelimit = 1200;
    
	double cur_obj;
		
	n = initializator->nodes;
		
	cout << "Numero nodi: " << n << endl;
		
        for (int tt = 1; tt <= test; tt++){
        
           
	        try
	        {
		        DECL_ENV( env );
		        DECL_PROB( env, lp );
		
		        //n=2;
		        //xy[][] = {{1.0,1.0},{1.0,1.0}};
		
		        setupLP(env, lp);

			// Serve per impostare un limite di tempo al solver così che non cicli all'infinito per trovare una soluzione
		        CPXsetdblparam(env, CPX_PARAM_TILIM, timelimit);
			
		        start = std::chrono::system_clock::now();
		        // optimize
		        CHECKED_CPX_CALL( CPXmipopt, env, lp );
		        end = std::chrono::system_clock::now();
		
		        // print objval
		        double objval;
		        CHECKED_CPX_CALL( CPXgetobjval, env, lp, &objval );
		        float elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()/1000.0f;
		        tc += elapsed_time;
		        // std::cout << "problem: " << p << "   test: " << tt << "   Objval: " << objval << "   time: " << elapsed_time << "s" << std::endl;
				
			cur_obj = objval;
				
		        CHECKED_CPX_CALL( CPXwriteprob, env, lp, "lel.lp", NULL);
		
		        //print solution (var values)
		        int n1 = CPXgetnumcols(env, lp);
		        if (n1 != n) { throw std::runtime_error(std::string(__FILE__) + ":" + STRINGIZE(__LINE__) + ": " + "different number of variables"); }
		        std::vector<double> varVals;
		        varVals.resize(n);
		        CHECKED_CPX_CALL( CPXgetx, env, lp, &varVals[0], 0, n - 1 );
		        /// status = CPXgetx (env, lp, x, 0, CPXgetnumcols(env, lp)-1);
		        for ( int i = 0 ; i < n ; ++i ) {
			        std::cout << "var in position " << i << " : " << varVals[i] << std::endl;
		        /// to get variable name, use the RATHER TRICKY "CPXgetcolname"
		        /// status = CPXgetcolname (env, lp, cur_colname, cur_colnamestore, cur_storespace, &surplus, 0, cur_numcols-1);
		        }
		        CHECKED_CPX_CALL( CPXsolwrite, env, lp, "pr.sol" );
		        // free
		
		        CPXfreeprob(env, &lp);
		        CPXcloseCPLEX(&env);
		
	        }
	        catch(std::exception& e)
	        {
		        // cout << ">>>EXCEPTION: " << e.what() << endl;
	        }
	
	    //cout << endl;
	    }//chiudo il for dei test
	  cout << "Av_Time " << p << ": " << tc / test << endl;
	  cout << "Best_Val " << p << ": " << cur_obj << endl;
          cout << endl;

	  ofstream SaveFile("results_TEST.txt", ios_base::app);
          SaveFile << p << "\t" << tc/test << "\t" << cur_obj << endl;
          SaveFile.close();

	}//chiudo il for dei problemi
	return 0;
}
