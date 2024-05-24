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

#include "Lattice.hpp"

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

	double pDep=double(M)*F;	//deposition probability
	TRandom3 rnd;

//approximation function
int approx(float a){				
	float appo= a-int(a);			
	if(appo<0.5)
		return floor(a);
	else
		return ceil(a);
}


//main
int main(int argc, char **argv){

	int N_max = approx(M*theta); 				//max number of filled sites
	cout << "n max= " << N_max << endl;

	TH2I grEvol("Evolution: step 0","",L,0,L,L,0,L);
	TApplication app("app",NULL, NULL);
	
    gStyle->SetOptStat(0);

	Lattice matrix(L,L,2);						//creates a lattice LxL with n particles
	
	//Graphics
	TCanvas can_dyn;
	can_dyn.SetTitle("Evolution");
	cout << "loading dynamical evolution of the system" << endl;
	matrix.PrintGr(grEvol);
	grEvol.SetTitle("Temporal evolution:step 0");
	grEvol.Draw("COL");
	matrix.PrintGr(grEvol);


	int counter=0;									//counter for the graph printing
	while(matrix.GetNPart()<N_max){				//while number of particles does not exceed the max fillment

		//deposition or diffusion -> draw a class 0 1 2 3 4
		//if diffusion draw unif in 1-n_mosse(C) and pick the corresponding move
		//if deposition matrix.Deposition()
		matrix.CreateClasses();						//divides the sites in classes according to #FN
		matrix.Growth(nu,T, E0, Eb, pDep);		//growth algorithm: diffusion or deposition
		
		if(counter%100==0){							//updates the graph every 100 steps
			matrix.PrintGr(grEvol);
			string title="Evoluzione Temporale: step "+to_string(counter);
			grEvol.SetTitle(title.c_str());
		}

		//cout << matrix.GetNPart()<< "\t" << matrix.CountNParticles()<<endl;
		++counter;
	}
	cout << "n iterations="<<counter<<endl;
	matrix.PrintGr(grEvol);
	cout << matrix.GetNPart()<< "\t" << matrix.CountNparticles()<<endl;
	cout <<"app.Run(true)"<< endl;
	app.Run(true);


	return 0;
}