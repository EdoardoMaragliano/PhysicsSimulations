/* AGGIORNATA AL 26 MAGGIO 2021
USO UN SOLO GENERATORE CON SEME INIZIALIZZATO A TIME(0)*/
/*MODELLO DDA DIFFONDO SOLO SE HO 0 VICINI*/

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

#include "ReticoloF.hpp"


using namespace std;

//PRINTERS

void Reticolo::Print(){
	for(int i=0;i<m_rows;i++){
		for(int j=0;j<m_cols;j++){
			cout << m_matr[i][j] << " ";
		}
		cout << endl;
	}
}

void Reticolo::Print(string nomefile){
	ofstream fout(nomefile);
	for(int i=0;i<m_rows;i++){
		for(int j=0;j<m_cols;j++){ 
			//fout << i << "\t" << j << "\t" << m_matr[i][j] << endl;
			fout << m_matr[i][j] << " ";
		}
		fout <<endl;
	}
	fout.close();
}


void Reticolo::PrintGr(TH2I &grP){

	grP.Reset();
			
	for(int i=0;i<m_rows;++i){
		for(int j=0;j<m_cols;++j){
			if(m_matr[i][j]==1){
				grP.SetBinContent(i+1,j+1,2);		//TH1I numera i bin da 1 invece che da 0
			}
		}
	}
	gPad->Modified(); gPad->Update(); gSystem->ProcessEvents(); 
}

void Reticolo::Init(){
  
    // Inserting elements into vector 
    for (int i = 0; i < m_rows; i++) { 
        vector<int> v_appo;
  
        for (int j = 0; j < m_cols; j++) { 
            v_appo.push_back(0); 
        } 
        m_matr.push_back(v_appo); 
    } 

	cout << "Reticolo inizializzato a 0" << endl;
}

int Reticolo::ContaParticelle(){
	int conta=0;
	for(int i=0;i<m_rows;i++){
		for(int j=0;j<m_cols;j++){
			if(m_matr[i][j]==1)
				conta++;
		}
	}
	return conta;
}


int Reticolo::NBounds(int x, int y){
	int N=0; //do per scontato che x,y sia un sito pieno

	/*int nx=x+1,ny=y+1,px=x-1,py=y-1;
	//implemento condiz al contorno cicliche(si puo' fare meglio, ma funziona)
	if(x==m_rows-1)
		nx=0;
	if(x==0)
		px=m_rows-1;
	if(y==m_cols-1)
		ny=0;
	if(y==0)
		py=m_cols-1;*/

	//uso fne PC per le cond periodiche (26.05.21)
	//conto i legami della particella nello stato  xy
	if(m_matr[PC(x-1)][y]==1)
		N++;
	if(m_matr[x][PC(y-1)]==1)
		N++;
	if(m_matr[PC(x+1)][y]==1)
		N++;
	if(m_matr[x][PC(y+1)]==1)
		N++;
	return N;
}

int Reticolo::PC(int num){           
  if (num == m_rows){num = 0;}
  else if (num == -1){num = m_rows-1;}
  return num;
}


void Reticolo::CreaClassi(){
	for(int cl=0;cl<4;++cl)								//ogni volta che la chiamo 
		m_classi[cl].clear();							//ripulisco le classi

	for(int i=0;i<m_rows;i++){							//ciclo sul reticolo
		for(int j=0;j<m_cols;j++){
			if(m_matr[i][j]==1){						//se il sito è pieno
				vector<int> appo={i,j};					//salvo le coordinate del sito	
				if(NBounds(i,j)==0){					//se ho 0 legami		
					m_classi[0].push_back(appo);		//metto il punto nella classe 0	
				}			
			}
		}
	}
}

void Reticolo::Crescita(double nu, double DovF, double pDep, bool flag){

	double p0=(4*m_classi[0].size()*DovF);																
	
	double pTot=0;				//peso totale P(C)
	pTot=p0+pDep;
	
	//scelgo la classe
	double r=m_rand.Rndm()*pTot;	//0-p0-p0+pdep


	if(r<p0){
		int k=m_rand.Integer(m_classi[0].size());	//numero di processi nella classe c scelta

		int x=m_classi[0][k][0];					//scelgo un elemento della classe
		int y=m_classi[0][k][1];					
		this->Diffusione(x,y);
	} else if (flag==true){
		//depongo
		if(m_Npart<0.1*m_rows*m_cols+1)
			this->Deposizione();
		else{cout << "chiudo tutto" << endl;
			return;}
	}
}



void Reticolo::Diffusione(int x, int y){	

	//ho rimosso i vincoli sull'energia interni alle classi
	//uso la funzione PC invece delle variabili ausiliarie(26.05.21)
	
	/*int nextx=x+1,nexty=y+1,prevx=x-1,prevy=y-1;
	//implemento condiz al contorno cicliche(si puo' fare meglio, ma funziona)
	if(x==m_rows-1)
		nextx=0;
	if(x==0)
		prevx=m_rows-1;
	if(y==m_cols-1)
		nexty=0;
	if(y==0)
		prevy=m_cols-1;*/
	
	if(m_matr[x][y]==1){						//se il sito e' pieno

		int dir=m_rand.Integer(4); 				//estraggo una direzione 
		if(dir==0){								//mi sposto in su
			if (m_matr[PC(x-1)][y]==0){ 			//se è vuoto
				m_matr[x][y]=0;
				m_matr[PC(x-1)][y]=1;
			}		
		}
		if(dir==1){								//mi sposto a dx
			if (m_matr[x][PC(y+1)]==0){			
				m_matr[x][y]=0;
				m_matr[x][PC(y+1)]=1;
			}			
		}
		if(dir==2){								//mi sposto in giu
			if (m_matr[PC(x+1)][y]==0){ 
				m_matr[x][y]=0;
				m_matr[PC(x+1)][y]=1;
			}		
		}
		if(dir==3){								//mi sposto a sx
			if (m_matr[x][PC(y-1)]==0){	
				m_matr[x][y]=0;
				m_matr[x][PC(y-1)]=1;	
			}		
		}
	}
	//cout << "diffusione effettuata" << endl;
}

void Reticolo::Deposizione(){
	int x=m_rand.Integer(m_rows); 							//estraggo un sito a caso	
	int y=m_rand.Integer(m_cols);
	while(m_matr[x][y]==1){								//finchè il sito estratto è pieno
		x=m_rand.Integer(m_rows); 							//continuo ad estrarre		
		y=m_rand.Integer(m_cols);							
	}
	m_matr[x][y]=1;					//quando esce dal ciclo(cioe' il sito è vuoto), lo riempio
	m_Npart++;						//se deposito, aggiorno il numero di particelle
}

void Reticolo::RandomFill(int Nfilled){
	int conta=0;
	int x=0, y=0;
	
	while(conta<Nfilled){
		x=m_rand.Integer(m_rows);		//estaggo indice tra 0 e L-1
		y=m_rand.Integer(m_cols);	

		if(m_matr[x][y]==0){  			//se la coordinata è vuota
			m_matr[x][y]=1;	   			//riempila	
			conta++;				  	//e aggiorna il conto
			m_Npart++;					//aggiorna il contatore interno alla classe
		}	
	}
	cout << "Reticolo riempito correttamente, nPart="<< m_Npart << endl;
}

void Reticolo::ComputeParOrdine(){
	float N_A=0, N_B=0;
	//fisso la riga
	for(int i=0;i<m_rows;i++){
		for(int j=0;j<m_cols;j++){
			if(m_matr[i][j]==1){
				if((i+j)%2==0)		//entrambi pari o entrambi dispari
					N_A++;
				if((i+j)%2==1)		//un pari e un dispari
					N_B++;
			}
		}
		m_p=(N_A-N_B)/(N_A+N_B);
	}	
}



//GETTERS

vector<vector<int>> Reticolo::GetMatrix(){
	return m_matr;
}

int Reticolo::GetNPart(){
	return m_Npart;
}

double Reticolo::GetParOrdine(){
	return m_p;
}

vector<vector<int>>* Reticolo::GetClassi(){
	return m_classi;
}

double Reticolo::GetRandom(){
	return m_rand.Integer(m_cols);
}