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
#include <unistd.h>

#include "Reticolo.hpp"

using namespace std;

	int L=64;							//side of the square lattice
	int M=L*L;							//total number of sites in the lattice
	//double kb=1.380649E-23;

	double theta=0.1; 					//max fillment %
		
	double T=200;						//interazione attrattiva
	double nu=1E12;
	double E0=0.4*11604; //eV->K 		//E0 barrier for an atom without First Neighbours
	double Eb=0.2*11604; //eV->K 		//binding energy
	double F=1./60;						//deposition flux (particles per sec per site)

	double pDep=double(M)*F;	//peso della deposizione costante
	TRandom3 rnd;

//ALCUNE FUNZIONI
int approx(float a){				//non trovavo una funzione che approssimasse alla parte
	float appo= a-int(a);			//intera piu' vicina cosi' l'ho scritta
	if(appo<0.5)
		return floor(a);
	else
		return ceil(a);
}


//MAIN
int main(int argc, char **argv){

	int N_max = approx(M*theta); 				//max number of filled sites
	cout << "n max= " << N_max << endl;

	TH2I grEvol("Evolution: step 0","",L,0,L,L,0,L);
	TApplication app("app",NULL, NULL);
	
    gStyle->SetOptStat(0);

	Reticolo matrice(L,L,2);						//creates a lattice LxL with n particles
	
	//Graphics
	TCanvas can_dyn;
	can_dyn.SetTitle("Evolution");
	cout << "loading dynamical evolution of the system" << endl;
	matrice.PrintGr(grEvol);
	grEvol.SetTitle("Temporal evolution:step 0");
	grEvol.Draw("COL");
	matrice.PrintGr(grEvol);


	int conta=0;									//counter for the graph printing
	while(matrice.GetNPart()<N_max){				//while number of particles does not exceed the max fillment

		//deposizione oppure diffusione -> estraggo una classe 0 1 2 3 4
		//se diffusione estraggo unif in 1-n_mosse(C) e scelgo la mossa corrispondente
		//se deposizione matrice.Deposition()
		matrice.CreaClassi();						//divides the sites in classes according to #FN
		matrice.Crescita(nu,T, E0, Eb, pDep);		//growth algorithm: diffusion or deposition
		
		if(conta%100==0){							//updates the graph every 100 steps
			matrice.PrintGr(grEvol);
			string title="Evoluzione Temporale: step "+to_string(conta);
			grEvol.SetTitle(title.c_str());
		}

		//cout << matrice.GetNPart()<< "\t" << matrice.ContaParticelle()<<endl;
		++conta;
	}
	cout << "n iterations="<<conta<<endl;
	matrice.PrintGr(grEvol);
	cout << matrice.GetNPart()<< "\t" << matrice.ContaParticelle()<<endl;
	cout <<"app.Run(true)"<< endl;
	app.Run(true);


	return 0;
}