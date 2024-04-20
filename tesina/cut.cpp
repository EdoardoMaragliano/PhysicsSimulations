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

/*dE/dx + TOF
      * simulazione dEdx e TOF ad impulso fissato della particella
		studio le medie troncate
		*/

TRandom3 rnd;

//////////////////////////////////////////////////////////////////////
TDatabasePDG pdg;
TParticlePDG *kp  = pdg.GetParticle(321);      //k+ 
TParticlePDG *k0  = pdg.GetParticle(311);      //k0
TParticlePDG *km  = pdg.GetParticle(-321);     //k-
TParticlePDG *pip = pdg.GetParticle("pi+");   //code 211
TParticlePDG *pim = pdg.GetParticle("pi-");
TParticlePDG *pi0 = pdg.GetParticle("pi");
TParticlePDG *Dp  = pdg.GetParticle(413);     //D*+
TParticlePDG *D0  = pdg.GetParticle(421);     //D0
TParticlePDG *Dm  = pdg.GetParticle(-413);    //D*-

//////////////////////////////////////////////////////////////////////

using namespace std;

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



  double Z=29;//82;                                      //numero atomico materiale, ora è Cu //Pb
  double A=63.546;//208;      
  double L=10;                                         //in cm
  double rho=8.94;//11.29;                             //density g/cm3
  double c=TMath::Ccgs();
  int nEv=1000;                                   //ciascuna misura
  int nSamp=20;                                   //viene da 20 campionamenti (e.g. 20 strati)
  double pbeam = 4;//rnd.Uniform(3,5);           //impulso iniziale del fascio in GeV



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

  rnd.SetSeed(1627479);
	TLorentzVector *in_p=new TLorentzVector(0,0,pbeam,sqrt(pbeam*pbeam+pip->Mass()*pip->Mass()));
	Particle part(pip,in_p);

	 TApplication aptest("aptest",0,NULL);
  TH1D *Emeas = new TH1D("E(last measurment)","E(last measurment)",100,50,100);
  Emeas     ->  GetXaxis()->SetTitle("dEdx [MeV/cm]");
  //GRAFICI MEDIE TRONCATE
  TH1D *h1=new TH1D("Ecut 90%","Ecut 90%",100,0,0);
  TH1D *h2=new TH1D("Ecut 80%","Ecut 80%",100,0,0);
  TH1D *h3=new TH1D("Ecut 70%","Ecut 70%",100,0,0);
  TH1D *h4=new TH1D("Ecut 60%","Ecut 60%",100,0,0);
  TH1D *h5=new TH1D("Ecut 55%","Ecut 55%",100,0,0);
  TH1D *h6=new TH1D("Ecut 50%","Ecut 40%",100,0,0);
  TH1D *h7=new TH1D("Ecut 47.5%","Ecut 447%",100,0,0);
  TH1D *h8=new TH1D("Ecut 45%","Ecut 45%",100,0,0);
  TH1D *h9=new TH1D("Ecut 42.5%","Ecut 42.5%",100,0,0);
  TH1D *h10=new TH1D("Ecut 40%","Ecut 40%",100,0,0);


  TH1D hCut[10]={*h1,*h2,*h3,*h4,*h5,*h6,*h7,*h8,*h9,*h10};
  for(auto h:hCut){
  	h.GetXaxis()->SetTitle("dEdx [MeV/cm]");
  }


  TGraph *grRMS=new TGraph();
  grRMS->SetTitle("RMS");
  TGraph *grDiff=new TGraph();
  grDiff->SetTitle("|E_{reco}-E_{gen}|");
  grDiff->GetXaxis()->SetTitle("troncamento");
  grRMS->GetXaxis()->SetTitle("troncamento");
  grDiff->SetMarkerStyle(7);
  grRMS->SetMarkerStyle(7);
 // grDiff->SetMarkerSize(5);
 // grRMS->SetMarkerSize(5);

  //NB in PDG le masse sono in GeV
  //NB in PDG le cariche sono in abs(e)/3
 
  double dEdx;
  double cut[10]={0.9,0.80,0.70,0.6,0.55,0.50,0.475,0.45,0.425,0.40};            //percentuali troncamento
  vector <double> data;
  double BB=BetheBloch(part,Z,A);
  cout <<endl<<endl<<"BetheBloch="<<BB<<endl<<endl;

	for(int n=0; n<nEv;++n){
	    data.clear();
	    
	    for(int i=0;i<nSamp;++i){                     //genero 20 campionamenti e ottengo la Landau
	      dEdx=rnd.Landau(BB,1);
	      data.push_back(dEdx); 
	    }
	    sort(data.begin(),data.end());               //ordino i dati

	    for(int j=0;j<10;j++){                       //per ciascun troncamento
	      for(int i=0;i<cut[j]*data.size();++i){    
	        Emeas->Fill(data[i]);                       //faccio il grafico troncato
	      }
	      double mean=Emeas->GetMean();          //estraggo la media dal gr troncato
	      Emeas->Reset();                        //resetto il grafico di Emeas                   
	      hCut[j].Fill(mean);                    //riempio hCut con la media
	    }
  	}
    for(int i=0;i<10;++i)
      hCut[i].Fit("gaus","L");
  for(auto mis:data)
    Emeas->Fill(mis);     //riempio con l'ultima misura
  for(int j=0;j<10;++j){
    grRMS->SetPoint(grRMS->GetN(),cut[j],hCut[j].GetStdDev());    //grafico rms vs troncamento
    grDiff->SetPoint(grDiff->GetN(),cut[j],abs(BB-hCut[j].GetMean()));
   // grDiff->SetPointError(grDiff->GetN(),)
  }

  TCanvas *can1=new TCanvas("Emeas","Emeas",500,500);
  Emeas->Draw();
  TCanvas *can2=new TCanvas("RMS","RMS",800,600);
  can2->Divide(2);
  can2->cd(1);
  grDiff->Draw("AP");
  can2->cd(2);
  grRMS->Draw("AP");
  TCanvas *can3=new TCanvas("Medie troncate","Medie Troncate",800,600);
  can3->Divide(4);
  for(int i=0;i<4;i++){
    can3->cd(i+1);
    hCut[2*i].Draw();
  }
  TCanvas *can4=new TCanvas("55cut","55cut",500,500);
  hCut[4].Draw();

  cout << "uso cut off al 55%" << endl;

  double TOF=L/TMath::Ccgs()*sqrt(1+pow(part.Mass()*TMath::Ccgs()/part.P(),2));
  cout << "TOF=" <<TOF<<endl;
  double tres=100e-12; //sec
  
  TOF=rnd.Gaus(TOF,tres);
  cout << "TOF measurment="<<TOF<<endl;
  aptest.Run(true);


  return 0;
}
