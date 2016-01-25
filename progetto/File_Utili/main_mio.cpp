/**
 * @file antenne.cpp
 * @brief 
 */

#include <cstdio>
#include <iostream>
#include <vector>
#include "cpxmacro.h"

using namespace std;

// error status and messagge buffer
int status;
char errmsg[BUF_SIZE];

// data
const int N = 5;    
const int A = 6;
const char nameN[N] = { 'A', 'B', 'C', 'D' , 'E' }; // Insieme dei nodi, rappresentano le posizioni dei fori
const char nameA[A] = { 'AB', 'BC', 'CD', 'DE', 'AC', 'CE' }; // Insieme degli archi che collegano i vari nodi, rappresentano il tragitto percorso dalla punta

const double C[N*N] = {	 // Tempo impiegato per lo spostamento della punta da un nodo ad un altro per ogni coppia di nodi i e j appartenenti ad N
  1.0, 1.0, 1.0, 1.0, 1.0,
  1.0, 1.0, 1.0, 1.0, 1.0,
  1.0, 1.0, 1.0, 1.0, 1.0,
  1.0, 1.0, 1.0, 1.0, 1.0,
  1.0, 1.0, 1.0, 1.0, 1.0,
};

const char root = 'A'; // Nodo di partenza 

const int flux = N; // Flusso iniziale uscente dalla radice

const int NAME_SIZE = 512;
char name[NAME_SIZE];
	
void setupLP(CEnv env, Prob lp, int & numVars )
{
  //intialize parameters
  
  std::cout << endl;
	// add x vars
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			char xtype = 'I';
			double obj = 0;
			double lb = 0;
			double ub = N;
			snprintf(name, NAME_SIZE, "x_%c_%c", nameN[i], nameN[j]);
			char* xname = (char*)(&name[0]);
			CHECKED_CPX_CALL( CPXnewcols, env, lp, 1, &obj, &lb, &ub, &xtype, &xname );
		}
	}
	// add y vars
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0, j < N, ++j)
		{
			char ytype = 'B';
			double obj = 1.0;
			double lb = 0.0;
			double ub = 1.0;
			snprintf(name, NAME_SIZE, "y_%c_%c", nameN[i], nameN[j]);
			char* yname = (char*)(&name[0]);
			CHECKED_CPX_CALL( CPXnewcols, env, lp, 1, &obj, &lb, &ub, &ytype, &yname );
		}
	}
	
	numVars = CPXgetnumcols(env, lp);
	
	// Aggiungo vincoli variabili
	
	// Vincoli variabili X
	
	
	
	
}

int main (int argc, char const *argv[])
{
	try
	{
		// init
		DECL_ENV( env );
		DECL_PROB( env, lp );
		// setup LP
		int numVars;
		setupLP(env, lp, numVars);
		// optimize
		CHECKED_CPX_CALL( CPXmipopt, env, lp );
		// print
		double objval;
		CHECKED_CPX_CALL( CPXgetobjval, env, lp, &objval );
		std::cout << "Objval: " << objval << std::endl;
		int n = CPXgetnumcols(env, lp);
		cout << n << " " << numVars << endl;
		if (n != numVars) { throw std::runtime_error(std::string(__FILE__) + ":" + STRINGIZE(__LINE__) + ": " + "different number of variables"); }
	  std::vector<double> varVals;
	  varVals.resize(n);
  	CHECKED_CPX_CALL( CPXgetx, env, lp, &varVals[0], 0, n - 1 );
  	for ( int i = 0 ; i < n ; ++i ) {
  	  std::cout << "var in position " << i << " : " << varVals[i] << std::endl;
  	}
		CHECKED_CPX_CALL( CPXsolwrite, env, lp, "antenne.sol" );
		// free
		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(&env);
	}
	catch(std::exception& e)
	{
		std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
	}
	return 0;
}
