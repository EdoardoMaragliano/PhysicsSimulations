/* AGGIORNATA AL 9 APRILE 2021
USO UN SOLO GENERATORE CON SEME INIZIALIZZATO A TIME(0)*/


#ifndef _Reticolo
#define _Reticolo

#include <iostream>
#include <ctime>
#include <iomanip>
#include <vector>
#include <fstream>
#include <map>

#include "TRandom3.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include <functional>
#include <unistd.h>
#include "TH2I.h"

using namespace std;


class Reticolo{
	
	public:
		Reticolo(int rows=10, int cols=10, double filled=0):m_rows(rows),m_cols(cols),m_Npart(filled){
			m_Npart=0;
			m_rand.SetSeed(time(0));
			
			cout << endl<< "created matrix "<< m_rows << " x " << m_rows << endl;
			cout << "random variables set"	<<endl;
			this->Init();
			this->RandomFill(filled);
			cout << "matrix filled with " << m_Npart << "particles" << endl;
		};

		vector<vector<bool>> GetMatrix();

		void Print();								//stampa in terminale
		void Print(string);							//stampa su file
		void PrintGr(TH2I &);

		void Init();				//inizializes the lattice to 0
		void RandomFill(int);		//given a % of fillment, randomly fills the lattice
		void Diffusione(int, int);	//makes a diffusion move from (x,y)	
		void Deposizione();			//fills an empty place of the lattice
		void ComputeParOrdine();
		void CreaClassi();			//divides points into classes depending on first neighbors
		void Crescita(double,double, double, double, double); //random growth MC algorithm

		int NBounds(int, int );		//given a position, counts the bonds
		int ContaParticelle();		//counts the filled places of the lattice
		int GetNPart();				//gives the value of the prv member m_nPart
		
		double GetParOrdine();
		vector<vector<int>>* GetClassi();	//gets m_classi

	private:
		int m_rows, m_cols;
		double m_Npart;
		vector<vector<bool>> m_matr;		//matrix
		vector<vector<int>> m_classi[5];	//m_classi[i] contains a vector of coordinates couples (x,y)
		TRandom m_rand;
		double m_p; 					//order parameter
	
};


#endif