// AGGIORNATA AL 26 MAGGIO 2021

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

#include "ReticoloF.hpp"


using namespace std;

	int L=64;
	int M=L*L;							//numero siti totali
	double theta=0.1; 					//riempimento percentuale massimo
		
	//double T=159;						
	double nu=1E12;						//interazione attrattiva
	double E0=0.4*11604; //eV->K 		//E0 barriera per un atomo senza primi vicini	
	double pDep=double(M);	//peso della deposizione costante
	TRandom3 rnd;

//ALCUNE FUNZIONI
int approx(float a){				//non trovavo una funzione che approssimasse alla parte
	float appo= a-int(a);			//intera piu' vicina cosi' l'ho scritta
	if(appo<0.5)
		return floor(a);
	else
		return ceil(a);
}

/*DovF=	1000 tante isole
		1E09 una sola
*/


//MAIN
int main(int argc, char **argv){
	if(argc==1) {
		cout <<"inserire D/F"<<endl;
		return -1;
	}
	double DovF=atof(argv[1]);
	

	cout << "D/F=" << DovF << endl;


	int N_max = approx(M*theta); 						//numero massimo siti pieni
	cout << "n max= " << N_max << endl;

	TH2I grEvol("Evoluzione: step 0","",L,0,L,L,0,L);	//grafico
	TApplication app("app",NULL, NULL);					
    gStyle->SetOptStat(0);
    TCanvas can_dyn;
	can_dyn.SetTitle("Evoluzione");

	Reticolo matrice(L,L,2);							//creo reticolo con 2 particelle
	
	cout << "carico evoluzione dinamica del sistema" << endl;
	matrice.PrintGr(grEvol);
	grEvol.SetTitle("Evoluzione temporale:step 0");
	grEvol.Draw("COL");
	matrice.PrintGr(grEvol);


	int conta=0;
	while(matrice.GetNPart()<N_max+1){				//finchè non raggiungo il riempimento massimo

		matrice.CreaClassi();
		matrice.Crescita(nu,DovF, pDep, true);
		
		if(conta%1000==0){							//stampa ogni tot cicli
			//usleep(100000);
			matrice.PrintGr(grEvol);
			string title="Evoluzione Temporale: step "+to_string(conta);
			grEvol.SetTitle(title.c_str());
		}

		//cout << matrice.GetNPart()<< "\t" << matrice.ContaParticelle()<<endl;
		++conta;
	}
	
	int extra=10000;
	cout << "continuo a diffondere ancora un poco"<< endl;
	for(int i=0;i<extra;++i){
		matrice.Crescita(nu,DovF,pDep,false);
	}
	cout << "n iterazioni="<<conta<<endl;
	matrice.PrintGr(grEvol);
	string title="D/F="+to_string(int(DovF));
	grEvol.SetTitle(title.c_str());
	can_dyn.Modified();
	can_dyn.Update();

	matrice.Print("matriceFin.dat");				//esporto matrice in file di testo
	cout <<"app.Run(true)"<< endl;
	
	app.Run(true);


	return 0;
}