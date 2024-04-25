/* AGGIORNATA AL 9 APRILE 2021
USO UN SOLO GENERATORE CON SEME INIZIALIZZATO A TIME(0)*/

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

#include "Reticolo.hpp"


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
			fout << i << "\t" << j << "\t" << m_matr[i][j] << endl;
		}
		fout << endl;
	}
	fout.close();
}


void Reticolo::PrintGr(TH2I &grP){			//fills and prints the graph

	grP.Reset();
			
	for(int i=0;i<m_rows;++i){
		for(int j=0;j<m_cols;++j){
			if(m_matr[i][j]==true){
				grP.SetBinContent(i+1,j+1,50);		//TH1I starts labeling bin from 1, not from 0
			}
		}
	}
	gPad->Modified(); gPad->Update(); gSystem->ProcessEvents(); 
}

void Reticolo::Init(){
  
    // Inserting elements into vector 
    for (int i = 0; i < m_rows; i++) { 
        vector<bool> v_appo;
  
        for (int j = 0; j < m_cols; j++) { 
            v_appo.push_back(false); 
        } 
        m_matr.push_back(v_appo); 
    } 

	cout << "Lattice initialized to 0" << endl;
}

int Reticolo::ContaParticelle(){
	int conta=0;
	for(int i=0;i<m_rows;i++){
		for(int j=0;j<m_cols;j++){
			if(m_matr[i][j]==true)
				conta++;
		}
	}
	return conta;
}



int Reticolo::NBounds(int x, int y){
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

void Reticolo::CreaClassi(){
	for(int cl=0;cl<4;++cl)								//each time I call the method 
		m_classi[cl].clear();							//I clean the classes

	for(int i=0;i<m_rows;i++){							//loop on the lattice
		for(int j=0;j<m_cols;j++){
			if(m_matr[i][j]==true){						//if the site is filled
				vector<int> appo={i,j};					//save the coords of the site
				for(int cl=0;cl<4;++cl){	
					if(NBounds(i,j)==cl){				//if it has cl bonds		
						m_classi[cl].push_back(appo);	//save the site in the cl-bonds class	
					}	
				}	
			}
			/*else{
				vector<int> appo2={i,j};					//salvo le coordinate del sito
				m_classi[4].push_back(appo2);
			}*/
		}
	}
	/*cout << "processi divisi in classi" << endl;
	for(int i=0;i<5;++i)
		cout << "dim classe " << i << "=" << m_classi[i].size() << endl;	*/
}

void Reticolo::Crescita(double nu, double T, double E0, double Eb, double pDep){

	vector<double> p;													//each diffusion class has a weight p[i]
	for(int i=0;i<4;++i){												//loop on the 4 diffusion classes
		p.push_back(m_classi[i].size()*(4-i)*nu*exp(-(E0+i*Eb)/T));		
		//cout <<	"classe"<<i<<"=" <<m_classi[i].size() << endl;
	}			
	p.push_back(pDep);											
	//cout << "p size=" << p.size() << endl;						//deposition weight is contant
	
	double pTot=0;													// total weight of the configuration P(C)
	for(int j=0; j<p.size();++j){
		//cout << "p[" <<j <<"]=" << p[j] << endl;
		pTot+=p[j];	
	}
	//cout << "pTot=" << pTot << endl;

		double p0,p1,p2,p3;
		p0=(p[0]); p1=(p0+p[1]);p2=(p1+p[2]);p3=(p2+p[3]);
		
		
	//random choice of the class
	int c=0;
	double r=m_rand.Rndm()*pTot;
	//cout << "estratto="<<r<<endl;

	if(0<r<p0){c=0;}					//0 First neighbors
	else if(p0<r<p1){c=1;}				//1 First neighbors
	else if(p1<r<p2){c=2;}				//2 First neighbors	
	else if(p2<r<p3){c=3;}				//3 First neighbors	
	else if(r>p3){c=4;}					//deposition

	//random choice of the move among possible moves of the chosen class
	if(c<4){ 										//if the process is diffusive
		int k=m_rand.Integer(m_classi[c].size());	//picks a rnd integer between 0 and the size of the chosen

		//cout << "c="<<c << "classi[c].size()="<<m_classi[c].size() << endl;
		//cout << "processo k=" << k << endl;

		int x=m_classi[c][k][0];					//the number labels a particular element of the class(thus a lattice site)
		int y=m_classi[c][k][1];					//saves the x and y coord of the site
		this->Diffusione(x,y);						//makes a diffusion move from the site
	}
	else{
		this->Deposizione();						//deposition process
	}
}

void Reticolo::Diffusione(int x, int y){			//algorithm for the diffusion process

	//0 su 1 dx 2 giu 3 sx
	
	int nextx=x+1,nexty=y+1,prevx=x-1,prevy=y-1;
	//implements cyclical boundary conditions (may be done better, still it works)
	if(x==m_rows-1)
		nextx=0;
	if(x==0)
		prevx=m_rows-1;
	if(y==m_cols-1)
		nexty=0;
	if(y==0)
		prevy=m_cols-1;
	
	if(m_matr[x][y]==true){						//if the site is filled

		int dir=m_rand.Integer(4); 				//randomly picks a direction 
		if(dir==0){								//mi sposto in su
			if (m_matr[prevx][y]==false){ 		//se è vuoto
				m_matr[x][y]=false;
				m_matr[prevx][y]=true;
			}		
		}
		if(dir==1){								//mi sposto a dx
			if (m_matr[x][nexty]==false){			
				m_matr[x][y]=false;
				m_matr[x][nexty]=true;
			}			
		}
		if(dir==2){								//mi sposto in giu
			if (m_matr[nextx][y]==false){ 
				m_matr[x][y]=false;
				m_matr[nextx][y]=true;
			}		
		}
		if(dir==3){								//mi sposto a sx
			if (m_matr[x][prevy]==false){	
				m_matr[x][y]=false;
				m_matr[x][prevy]=true;	
			}		
		}
	}
	//cout << "diffusione effettuata" << endl;
}

void Reticolo::Deposizione(){
	int x=m_rand.Integer(m_rows); 							//randomly picks a site	
	int y=m_rand.Integer(m_cols);
	while(m_matr[x][y]==true){								//while the picked site is filled
		x=m_rand.Integer(m_rows); 							//keep picking		
		y=m_rand.Integer(m_cols);							
	}
	m_matr[x][y]=true;					//when out of the while loop(thus the site is empty), fill it
	m_Npart++;							//update the counter of m_Npart
}

void Reticolo::RandomFill(int Nfilled){		//randomly fills the lattice until Nfilled sites are filled
	int conta=0;
	int x=0, y=0;
	
	while(conta<Nfilled){
		x=m_rand.Integer(m_rows);		//estaggo indice tra 0 e L-1
		y=m_rand.Integer(m_cols);	

		if(m_matr[x][y]==false){  		//se la coordinata è vuota
			m_matr[x][y]=true;	   		//riempila	
			conta++;				  	//e aggiorna il conto
			m_Npart++;					//aggiorna il contatore interno alla classe
		}	
	}
	cout << "lattice correctly filled, nPart="<< m_Npart << endl;
}

void Reticolo::ComputeParOrdine(){
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
		m_p=(N_A-N_B)/(N_A+N_B);
	}	
}



//GETTERS

vector<vector<bool>> Reticolo::GetMatrix(){
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