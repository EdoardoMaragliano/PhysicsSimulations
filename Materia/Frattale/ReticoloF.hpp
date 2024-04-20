/* AGGIORNATA AL 26 MAGGIO 2021
USO UN SOLO GENERATORE CON SEME INIZIALIZZATO A TIME(0)*/


#ifndef _Reticolo
#define _Reticolo

#include <iostream>
#include <ctime>
#include <iomanip>
#include <vector>
#include <fstream>

#include "TRandom3.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TMatrix.h"
#include <functional>
#include <unistd.h>
#include "TH2I.h"

using namespace std;


class Reticolo{
	
	public:
		Reticolo(int rows=10, int cols=10, double piene=0):m_rows(rows),m_cols(cols),m_Npart(piene){
			m_Npart=0;
			m_rand.SetSeed(time(0));
			
			cout << endl<< "creata matrice "<< m_rows << " x " << m_rows << endl;
			cout << "variabili casuali settate"	<<endl;
			this->Init();
			this->RandomFill(piene);
			cout << "matrice riempita con " << m_Npart << " particelle" << endl;
		};

		void Print();								//stampa in terminale
		void Print(string);							//stampa su file
		void PrintGr(TH2I &);

		void Init();								//inizializza a 0
		void RandomFill(int);						//riempie a caso dato un riempimento
		void Diffusione(int, int);					//fa una mossa di diffusione	
		void Deposizione();							//riempie un sito vuoto
		void ComputeParOrdine();
		void CreaClassi();							//divide i siti in classi
		void Crescita(double, double, double, bool flag);

		int NBounds(int, int );			//conta i legami, dato un sito
		int ContaParticelle();			//fa un ciclo sulla matrice e conta i siti pieni
		int GetNPart();					//restituisce il valore di m_Npart
		int PC(int num);

		double GetParOrdine();
		vector<vector<int>>* GetClassi();
		vector<vector<int>> GetMatrix();
		double GetRandom();

	private:
		int m_rows, m_cols;
		double m_Npart;
		vector<vector<int>> m_matr;			//matrice
		vector<vector<int>> m_classi[5];	//m_classi[i] contiene le coord dei siti della classe
		TRandom m_rand;
		double m_p; 					//parametro d'ordine
	
};


#endif