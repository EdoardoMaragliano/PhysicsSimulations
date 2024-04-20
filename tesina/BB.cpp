#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TApplication.h>
#include <iostream>
#include <iomanip>
#include <THStack.h>


using namespace std;

double Z=29;//82;                                      //numero atomico materiale, ora è Cu //Pb
  double A=63.546;//208;      
  double L=10;                                         //in cm
  double rho=8.94;//11.29;                             //density g/cm3
  double c=TMath::Ccgs();
  int nEv=1000;                                   //ciascuna misura
  int nSamp=20;                                   //viene da 20 campionamenti (e.g. 20 strati)
  double pbeam = 4;//rnd.Uniform(3,5);           //impulso iniziale del fascio in GeV

TDatabasePDG pdg;
TParticlePDG *pip = pdg.GetParticle("pi+");   //code 211

class Particle: public TLorentzVector{
public:
  using TLorentzVector::TLorentzVector;
  Particle(TParticlePDG p, TLorentzVector t): TLorentzVector(t), m_code(p.PdgCode()),m_charge(p.Charge()/3){
    SetVectM(t.Vect(),p.Mass());};
  Particle(TParticlePDG *p, TLorentzVector *t): TLorentzVector(*t),m_code(p->PdgCode()),m_charge(p->Charge()/3){
    SetVectM(t->Vect(),p->Mass());};
  int Charge(){return m_charge;}
  int PdgCode(){return m_code;}
  double Mass(){return this->M();}
private:
  int m_code, m_charge;
};


double BetheBloch(Particle p, double Z, double A ){
  //NB masse DPG in GeV
  double K=0.307;                       //Mev g^-1 cm^2    //da Taiuti
  double I=10;
  if(Z>=13){
   I=9.76*Z+58.8*pow(Z,-0.19);    //eV per Z>13 da Taiuti
  } else{
   I=12*Z+7; //eV per Z<13
  }
  //K=K*1e-3; //in GeV g^-1 cm^2
  I=I*1e-6;   //in MeV
  double argLog=2*0.511*pow(TMath::Ccgs()*p.Beta()*p.Gamma(),2)/I;
  //cout << p.Beta() << setw(10) <<  p.Gamma() << endl;
  double delta=0; //non trovo espressione per delta ma se uso basse energie no prob I guess
  return rho*K*Z/A*pow(p.Charge()/p.Beta(),2)*(log(argLog)-pow(p.Beta(),2)-delta/2);
}
//In MeV/cm

int main(){

	TApplication app("app",0,NULL);
	TCanvas can;
	can.SetTitle("BetheBloch");
	can.SetLogx();
	TGraph gr;
	gr.SetTitle("Bethe Bloch");
	gr.GetXaxis()->SetTitle("p");
	gr.GetYaxis()->SetTitle("dEdx");
	for(int i=1;i<1000;++i){
		double pbeam=0.1*i;
		TLorentzVector *in_p=new TLorentzVector(0,0,pbeam,sqrt(pbeam*pbeam+pip->Mass()*pip->Mass()));
		Particle pion(pip,in_p);
		double BB=BetheBloch(pion,Z,A);
		gr.SetPoint(gr.GetN(),i,BB);
	}

	gr.Draw("ACP");

	app.Run(true);
	return 0;
}