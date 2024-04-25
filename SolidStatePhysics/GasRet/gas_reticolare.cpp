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

#include "Matrice.hpp"

using namespace std;

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

	int L=64;
	double theta=0.5; //riempimento percentuale
	int N_mcs=150000;
	double BetaJ=5.0;
	int dyn=1;
	cout << "argc="<<argc<<endl;
/*	if(argc!=1){
		int L=atoi(argv[1]);
		double theta=atoi(argv[2]);
		int N_mcs = atoi(argv[3]);
		double BetaJ = atoi(argv[4]);
		int dyn=atoi(argv[5]);
	}
	*/

/* con L=128 ci mette tantissimo, ridurre frequenza di stampa
	con L=64 e seed(time(0)) servono un po' più di 100k MCS per avere convergenza
	*/

	int spostamenti=0;	//contatore
	TGraph gr_p;
	TH2I gr0("Reticolo iniziale","",L,0,L,L,0,L);
	TH2I gr1("Reticolo intermedio","",L,0,L,L,0,L);
	TH2I gr2("Reticolo finale","",L,0,L,L,0,L);
	TH2I grEvol("Evoluzione: step 0","",L,0,L,L,0,L);
	TApplication app("app",NULL, NULL);
	
	gr0.SetTitle("Reticolo iniziale");
	gr1.SetTitle("Reticolo intermedia");
	gr2.SetTitle("Reticolo finale");

	gr_p.SetMarkerStyle(7);
	gr_p.SetTitle("Parametro d'ordine");
	gr_p.GetXaxis()->SetTitle("n");
	gr_p.GetYaxis()->SetTitle("p");

    gStyle->SetOptStat(0);
	TCanvas can;
	can.SetTitle("Riassunto");
	can.Divide(2,2);
	can.cd(1);
	gr0.Draw("COL");
	can.cd(2);
	gr1.Draw("COL");
	can.cd(3);
	gr2.Draw("COL");
	can.cd(4);
	gr_p.Draw("AL");

	while(L%2!=0){
		cout << "L deve essere pari, inserire L" << endl;
		cin >> L;
	}
	while(dyn!=1 && dyn!=0){
		cout << "vuoi l'evoluzione dinamica? 1/0"<<endl;
		cin >> dyn;
	}		

	int M=L*L;						//numero siti totali
	int N_part = approx(M*theta); 	//numero siti pieni
	cout << "n particelle= " << N_part << endl;

	Matrice matrice(L,L, BetaJ, N_part);

	matrice.Print("output0.dat");
	matrice.PrintGr(gr0);

	TCanvas can_dyn;
	can_dyn.SetTitle("Evoluzione");

	if(dyn==1){		//set del grafico dell'evoluzione
		cout << "carico evoluzione dinamica del sistema" << endl;
		matrice.PrintGr(grEvol);
		grEvol.SetTitle("Evoluzione temporale:step 0");
		grEvol.Draw("COL");
	}	

	for(int i=0;i<N_mcs;i++){
		for(int k=0;k<N_part;k++)
			matrice.RandomMove(spostamenti);					//sposta una particella
		if(i%250==0){		//ogni 500 mcs
			matrice.ComputeParOrdine();
			gr_p.SetPoint(gr_p.GetN(),i,matrice.GetParOrdine());
			gr_p.Draw("COL"); 		
		}
		if (dyn==1 && i%250==0){
			matrice.PrintGr(grEvol);
			string title="Evoluzione Temporale: step "+to_string(i);
			grEvol.SetTitle(title.c_str());
			grEvol.Draw("COL");
		}	
		if(i==13000){
			matrice.Print("output1.dat");
			matrice.PrintGr(gr1);
			gr1.Draw("COL");
		}		
	}
		
	//cout << "spostamenti= " << spostamenti << endl << endl;

	matrice.Print("output2.dat");
	matrice.PrintGr(gr2);
	gr2.Draw("COL");
	
	cout << "p medio= " << gr_p.GetMean(2) << endl;
	cout << "RMS p medio= " << gr_p.GetRMS(2)/sqrt(gr_p.GetN()) << endl;

	app.Run(true);


	return 0;
}