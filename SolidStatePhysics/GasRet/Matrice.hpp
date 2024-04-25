#ifndef _MATRICE
#define _MATRICE

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
#include <functional>
#include <unistd.h>
#include "TH2I.h"

using namespace std;


class Matrice{
	
	public:
		Matrice(int rows=10, int cols=10, double BetaJ=2.8, double piene=50):m_rows(rows),m_cols(cols),m_BetaJ(BetaJ){
			m_Npart=piene;
			m_lx.SetSeed(65789036);
			m_ly.SetSeed(987654321);
			m_dir.SetSeed(89786);
			cout << endl<< "creata matrice "<< m_rows << " x " << m_rows << endl;
			cout << "variabili casuali settate"	<<endl;
			this->Init();
			this->RandomFill(m_Npart);
		};

		vector<vector<bool>> GetMatrix();

		void Print();								//stampa in terminale
		void Print(string);							//stampa su file
		void PrintGr(TH2I &);

		void Init();				//inizializza a 0
		void RandomFill(int);		//riempie a caso dato un riempimento
		void RandomMove(int &);		//fa una mossa	

		int NBounds(int, int );		//conta i legami, dato un sito
		int ContaParticelle();
		
		void ComputeParOrdine();
		double GetParOrdine();

	private:
		int m_rows, m_cols;
		double m_Npart;
		vector<vector<bool>> m_matr;	//matrice
		TRandom m_lx,m_ly,m_dir;
		double m_p; 					//parametro d'ordine
		double m_BetaJ;					//agitazione termica
};


#endif