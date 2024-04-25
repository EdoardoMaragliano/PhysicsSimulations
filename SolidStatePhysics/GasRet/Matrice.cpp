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

vector<vector<bool>> Matrice::GetMatrix(){
	return m_matr;
}

int Matrice::ContaParticelle(){
	int conta=0;
	for(int i=0;i<m_rows;i++){
		for(int j=0;j<m_cols;j++){
			if(m_matr[i][j]==true)
				conta++;
		}
	}
	return conta;
}

void Matrice::Print(){
	for(int i=0;i<m_rows;i++){
		for(int j=0;j<m_cols;j++){
			cout << m_matr[i][j] << " ";
		}
		cout << endl;
	}
}

void Matrice::Print(string nomefile){
	ofstream fout(nomefile);
	for(int i=0;i<m_rows;i++){
		for(int j=0;j<m_cols;j++){ 
			fout << i << "\t" << j << "\t" << m_matr[i][j] << endl;
		}
		fout << endl;
	}
	fout.close();
}


void Matrice::PrintGr(TH2I &grP){

	grP.Reset();
			
	for(int i=0;i<m_rows;i++){
		for(int j=0;j<m_cols;j++){
			
			if(m_matr[i][j]==true && (i+j)%2==0){
				grP.SetBinContent(i,j,8);
			}
			if(m_matr[i][j]==true && (i+j)%2==1){
				grP.SetBinContent(i,j,2);
			}
		}
	}
	gPad->Modified();gPad->Update(); gSystem->ProcessEvents(); 
}

void Matrice::Init(){
  
    // Inserting elements into vector 
    for (int i = 0; i < m_rows; i++) { 
        vector<bool> v_appo;
  
        for (int j = 0; j < m_cols; j++) { 
            v_appo.push_back(false); 
        } 
        m_matr.push_back(v_appo); 
    } 

	cout << "matrice inizializzata a 0" << endl;
}

void Matrice::RandomFill(int Nfilled){
	int conta=0;
	int x=0, y=0;
	m_lx.SetSeed(time(0));
	//m_ly.SetSeed(time(0));
	
	while(conta<Nfilled){
		x=m_lx.Integer(m_rows);		//estaggo indice tra 0 e L-1
		y=m_ly.Integer(m_cols);	

		if(m_matr[x][y]==false){  		//se la coordinata è vuota
			m_matr[x][y]=true;	   		//riempila	
			conta++;				  	//e aggiorna il conto	
		}	
	}
	cout << "matrice riempita correttamente" << endl;
}


int Matrice::NBounds(int x, int y){
	int N=0; //do per scontato che x,y sia un sito pieno

	int nx=x+1,ny=y+1,px=x-1,py=y-1;
	//implemento condiz al contorno cicliche(si puo' fare meglio, ma funziona)
	if(x==m_rows-1)
		nx=0;
	if(x==0)
		px=m_rows-1;
	if(y==m_cols-1)
		ny=0;
	if(y==0)
		py=m_cols-1;
	//conto i legami della particella nello stato  xy
	if(m_matr[px][y]==true)
		N++;
	if(m_matr[x][py]==true)
		N++;
	if(m_matr[nx][y]==true)
		N++;
	if(m_matr[x][ny]==true)
		N++;
	return N;
}

void Matrice::RandomMove(int &successi){	
	//andrà iterata N_mosse volte
	//la variabile successi è aggiornata per referenza
	m_dir.SetSeed(time(0)); //0 su 1 dx 2 giu 3 sx
	
	
	int x=m_lx.Integer(m_rows); 							//estraggo un sito a caso	
	int y=m_ly.Integer(m_cols);
	int nextx=x+1,nexty=y+1,prevx=x-1,prevy=y-1;
	//implemento condiz al contorno cicliche(si puo' fare meglio, ma funziona)
	if(x==m_rows-1)
		nextx=0;
	if(x==0)
		prevx=m_rows-1;
	if(y==m_cols-1)
		nexty=0;
	if(y==0)
		prevy=m_cols-1;
	

	if(m_matr[x][y]==true){						//se il sito e' pieno

		int dir=m_dir.Integer(4); 				//estraggo una direzione 
		if(dir==0){								//mi sposto in su
			if (m_matr[prevx][y]==false){ 		//se è vuoto
				if(m_lx.Rndm()<exp(-m_BetaJ*double(this->NBounds(prevx,y)-this->NBounds(x,y)))){
					//cout << "up" << endl;
					//cout << "exp="<<exp(-m_BetaJ*double(this->NBounds(prevx,y)-this->NBounds(x,y)))<<endl;
					m_matr[x][y]=false;
					m_matr[prevx][y]=true;
					successi++;
				}
			}		
		}
		if(dir==1){								//mi sposto a dx
			if (m_matr[x][nexty]==false){
				if(m_lx.Rndm()<exp(-m_BetaJ*double(this->NBounds(x,nexty)-this->NBounds(x,y)))){
				//cout << "dx" << endl; 
				//cout << "exp="<<exp(-m_BetaJ*double(this->NBounds(x,nexty)-this->NBounds(x,y)))<<endl;			
					m_matr[x][y]=false;
					m_matr[x][nexty]=true;
					successi++;
				}	
			}			
		}
		if(dir==2){								//mi sposto in giu
			if (m_matr[nextx][y]==false){ 
				if(m_lx.Rndm()<exp(-m_BetaJ*double(this->NBounds(nextx,y)-this->NBounds(x,y)))){
					//cout << "down" << endl;
					//cout << "exp="<<exp(-m_BetaJ*double(this->NBounds(nextx,y)-this->NBounds(x,y)))<<endl;	
					m_matr[x][y]=false;
					m_matr[nextx][y]=true;
					successi++;
				}	
			}		
		}
		if(dir==3){								//mi sposto a sx
			if (m_matr[x][prevy]==false){
				if(m_lx.Rndm()<exp(-m_BetaJ*double(this->NBounds(x,prevy)-this->NBounds(x,y)))){ 
				//cout << "exp="<<exp(-m_BetaJ*double(this->NBounds(x,prevy)-this->NBounds(x,y)))<<endl;	
				//cout << "sx" << endl;		
					m_matr[x][y]=false;
					m_matr[x][prevy]=true;
					successi++;
				}	
			}		
		}
	}
}

void Matrice::ComputeParOrdine(){
	float N_A=0, N_B=0;
	//fisso la riga
	for(int i=0;i<m_rows;i++){
		for(int j=0;j<m_cols;j++){
			if(m_matr[i][j]==true){
				if((i+j)%2==0)		//entrambi pari o entrambi dispari
					N_A++;
				if((i+j)%2==1)		//un pari e un dispari
					N_B++;
			}
		}
		m_p=(N_A-N_B)/m_Npart;
	}	
}

double Matrice::GetParOrdine(){
	return m_p;
}