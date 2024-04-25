//g++ -g -o GasVideo GasVideo.cpp `root-config --cflags --glibs`

#include "TRandom3.h"
#include "TMatrix.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TSystem.h"                                  
#include "TStyle.h"
#include <cmath>
#include <iostream>
#include "unistd.h"  
#include <fstream>

using namespace std;
int MS = 100000;
double BetaJ = 3;
int L = 64;

int Np = 64*64*0.5;
double Theta = 0.5;


int pc(int num){                 //periodic conditions           
  if (num == L){num = 0;}
  else if (num == -1){num = L-1;}
  return num;
}

double Prob(int vicini){           
  return exp(-vicini*BetaJ);
}

int main(){
  TMatrix C(2,Np);								//cosa ci salva?
  TMatrix A(L,L);								//matrice LxL
  ofstream ofile1; 
  ofile1.open("Ordine.dat");
  TH2I Fine ("RETICOLO","",L,0,L,L,0,L);
  TApplication app("app",NULL,NULL);
  TCanvas c;

  Fine.Draw("COLZ");
  TRandom3 rnd;
  //rnd.SetSeed(1585222552);
  rnd.SetSeed(time(0));
  cout << "Seed: " << time(0) << endl;
  
  
  int n = 0; 
  int b,i,j,iold,jold,d; 
  double B,I,J,D, DN, val, NA,NB;

  //RIEMPIE LA MATRICE A CASO
  while(n < Np){
    I = rnd.Rndm()*L; i = (int) I; 
    J = rnd.Rndm()*L; j = (int) J;
    if(A(i,j)==0){
      A(i,j)=1; 
      C(0,n)= i; 
      C(1,n)= j;
       n++;}
  }

  cout << "Sto calcolando, mettiti comodo..." << endl;
  
  for(int c=0; c<MS; c++){		//MCS
		
	    for(int r=0; r<Np; r++){		//ciclo sulle particelle
			B = rnd.Rndm()*Np; b = (int) B;
			i = C(0,b); j = C(1,b);
			D = rnd.Rndm()*4; d = (int) D;	//estrae la direzione
			iold = i; jold = j; A(iold,jold)=0;
			
			if(d==0){i = pc(i+1);}			//scelta la direzione determina la casella successiva
			if(d==1){j = pc(j+1);}
			if(d==2){i = pc(i-1);}
			if(d==3){j = pc(j-1);}
			
			DN = 0;							//nfin-nin
			if (A(i,j)==1){
				A(iold,jold)=1; //se la casella successiva é piena lascia la precedente piena
				continue;
			}
			else {				//se è vuota conta i vicini
			  DN = (A(pc(i-1),j) + A(pc(i+1),j) + A(i,pc(j-1)) + A(i,pc(j+1)));		//basta sommare i valori di A(i,j)
			  DN = DN - (A(pc(iold-1),jold) + A(pc(iold+1),jold) + A(iold,pc(jold-1)) + A(iold,pc(jold+1)));
			  if (DN <= 0){	//se DeltaN è negativo exp(-DN)>1 quindi accetto sempre
			  	A(i,j)=1; 	//riempio la nuova casella
			  	C(0,b)=i; 	
			  	C(1,b)=j;
			  }
			  else {		//se exp(-DN)<0 devo estrarre valore casuale in (0,1)
			    val = rnd.Rndm();
			    if(val < Prob(DN)){	
			    	A(i,j)=1;
			    	C(0,b)=i;
			    	C(1,b)=j;
			    }
			    else {A(iold,jold)=1;}
			  }
			}
	    }	//FINE CICLO INTERNO (PARTICELLE)
	    
	    
	    NA=0;NB=0;							//contatori per parametro d'ordine
	    for(int k=0; k<L/2; k++){			//genialata: fa un ciclo su k/2, quindi fa k**2/4 operazioni 
	      for(int q=0; q<L/2; q++){			//invece di k**2 operazioni (direi)	
		NA += A(2*k,2*q) + A(2*k+1,2*q+1);	//conta solo colonna e riga entrambe pari o dispari	
		NB += A(2*k,2*q+1) + A(2*k+1,2*q);	//conta solo colonna pari e riga dispari e viceversa
	      }
	    }
	    ofile1 << c  << " " << (NA-NB)/Np << endl; //esporta su file il par d'ordine
	    
	    Fine.Reset();						//fa il reset del TH2I (comodo, per TGraph sta cosa non c'e')
	    for(int k=0; k<L; k++){
	      for(int q=0; q<L; q++){					//ciclo su tutta la matrice
	      	if(A(k,q)==1){							//se casella piena
	      		Fine.SetBinContent(k+1,q+1,2);		//la riempie(il +1 presumo dipenda da come TH2I indicizza)
	      	}										//
	      }
	    }
	    
	    for(int k=0; k<L/2; k++){			//colora diversamente i due reticoli (per la prima volta)
	      for(int q=0; q<L/2; q++){
			if(A(2*k,2*q)==0){Fine.SetBinContent(2*k+1,2*q+1,8);}
			if(A(2*k+1,2*q+1)==0){Fine.SetBinContent(2*k+1+1,2*q+1+1,8);}
			if(A(2*k,2*q+1)==0){Fine.SetBinContent(2*k+1,2*q+1+1,6);}
			if(A(2*k+1,2*q)==0){Fine.SetBinContent(2*k+1+1,2*q+1,6);}   
	      }
	    }
	    
	    Fine.GetXaxis()->SetLabelOffset(999); Fine.GetYaxis()->SetLabelOffset(999);	//non so esattamente, make-up
	    Fine.GetXaxis()->SetLabelSize(0);Fine.GetYaxis()->SetLabelSize(0);			
	    gStyle->SetOptStat("n");                                                                   
	    Fine.Draw("COLZ");												//stampa il grafico
	    /*A box is drawn for each cell with a color scale varying with contents. 
  		All the none empty bins are painted. Empty bins are not painted */
	    gPad->Modified();	gPad->Update();  gSystem->ProcessEvents();                                                     
			       
	    /*
	    BJ=1.7
	    usleep(100000);
	    */

	}	//FINE CICLO ESTERNO (MCS)
  
  TCanvas c2;
  c2.Divide(2,1);

  c2.cd(1);

  Fine.Reset();
  	for(int k=0; k<L; k++){
	    for(int q=0; q<L; q++){
	      	if(A(k,q)==1){
	      		Fine.SetBinContent(k+1,q+1,2);
	      	}
	    }
  	}
  
  	for(int k=0; k<L/2; k++){			//stessa cosa di prima, colora i due reticoli
    	for(int q=0; q<L/2; q++){
	      if(A(2*k,2*q)==0){Fine.SetBinContent(2*k+1,2*q+1,8);}
	      if(A(2*k+1,2*q+1)==0){Fine.SetBinContent(2*k+1+1,2*q+1+1,8);}
	      if(A(2*k,2*q+1)==0){Fine.SetBinContent(2*k+1,2*q+1+1,6);}
	      if(A(2*k+1,2*q)==0){Fine.SetBinContent(2*k+1+1,2*q+1,6);}   
    	}
  	}

  Fine.GetXaxis()->SetLabelOffset(999); Fine.GetYaxis()->SetLabelOffset(999);
  Fine.GetXaxis()->SetLabelSize(0);Fine.GetYaxis()->SetLabelSize(0);
  gStyle->SetOptStat("n");
  Fine.Draw("COL");
  /*A box is drawn for each cell with a color scale varying with contents. 
  All the none empty bins are painted. Empty bins are not painted */
  gPad->Modified();gPad->Update(); gSystem->ProcessEvents();
  
  
  c2.cd (2);
  TGraph gr("Ordine.dat");		//ha esportato su file il p ordine e costruisce da file il gr
  gr.Draw();
  gr.SetTitle("Parametro d'ordine");gr.GetXaxis()->SetTitle("MCS");
  
  ofile1.close();
  gPad->Modified();gPad->Update(); gSystem->ProcessEvents();
  app.Run(false);
  return 0;
}
