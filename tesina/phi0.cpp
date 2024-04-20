#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include <TRandom3.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TApplication.h>
#include <iostream>
#include <cmath>
#include <TGraph.h>



using namespace std;


TRandom3 rnd;

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
double c =3e10;
double Z = 14; //Cu
double A = 28.0855; //Cu massa atomica
double rho = 2.33; //densità
TH1D hE("E","dE/dx Landau->Gaussian",100,0,0);

class Particle: public TLorentzVector{
public:
  using TLorentzVector::TLorentzVector;
  Particle(TParticlePDG p, TLorentzVector t): TLorentzVector(t), m_code(p.PdgCode()),m_charge(p.Charge()){
    SetVectM(t.Vect(),p.Mass());};
  Particle(TParticlePDG *p, TLorentzVector *t): TLorentzVector(*t),m_code(p->PdgCode()),m_charge(p->Charge()){
    SetVectM(t->Vect(),p->Mass());};
  int Charge(){return m_charge;}
  int PdgCode(){return m_code;}
private:
  int m_code, m_charge;
};

double BetheBloch(Particle p, double Z, double A, double rho){
  double K = 0.307;
  double I = 173e-6;
  double arg = 2*0.5*pow(c*p.Beta()*p.Gamma(),2)/I;
  double delta = 0;
  return rho*K*Z/A*pow(p.Charge()/p.Beta(),2)*(log(arg)-pow(p.Beta(),2)-delta/2)*300e-4;
}

bool PId(Particle trk, TLorentzVector *t){
  
  
  bool id = false;
  
  double energy = BetheBloch(trk,Z,A,rho);
  energy = rnd.Gaus(energy,0.5);
  Particle pi(pip, t);
  Particle kaon(kp,t);
  
  

  double energyp = BetheBloch(pi,Z,A,rho);
  double energyk = BetheBloch(kaon,Z,A,rho);
  //  cout << "energy " << energyp << " " << energyk << endl;
  //  cout << "mass " << pi.M() << " " << kaon.M() << endl;
  double ka = TMath::Gaus(energy, energyk,0.5);
  double pion = TMath::Gaus(energy, energyp,0.5);
  //  cout <<"ka/pi" << ka/pion << endl;
  
  
  if(ka/pion>0.5)
    id = true;
  return id;
}

    

int main(){
  rnd.SetSeed(123456);
  TApplication app("app",0,NULL);
 

  TH1D hM("hM","",100,0,0);
  
  TGenPhaseSpace tgps;

  double pphi = 1e-2;
  int nEv = 1e5;
  double w;
  
  for(int i = 0; i<nEv;i++){
    
  
   
    vector<Particle> decayParticles;
    vector<TLorentzVector> decayP;

        //decadimento Jpsi
    TLorentzVector in(0.,0.,0.,4.040);       //ferma nel CM
    //cout <<"Jpsi -> D*p + D*m"<<endl;
    double masses[2]={Dp->Mass(),Dm->Mass()};
    tgps.SetDecay(in,2,masses);                       
    w=tgps.Generate();                              
    TLorentzVector *v_Dp=tgps.GetDecay(0);          //cinematica di D*p
    TLorentzVector *v_Dm=tgps.GetDecay(1);          //cinematica di D*m
    double prob = rnd.Rndm();
    if(prob<pphi){
      double masses[2]={D0->Mass(),pip->Mass()};
      tgps.SetDecay(*v_Dp,2,masses);
      w = w*tgps.Generate();

      Particle pionp(pip,tgps.GetDecay(1));
      decayParticles.push_back(pionp);
      TLorentzVector* a = tgps.GetDecay(1);
      decayP.push_back(*a);

      TLorentzVector* vphi = tgps.GetDecay(0);
      double mk[2] = {km->Mass(), pip->Mass()};
      tgps.SetDecay(*vphi,2,mk);
      w = w*tgps.Generate();

      Particle kaonp(kp,tgps.GetDecay(0));
      a = tgps.GetDecay(1);
      decayP.push_back(*a);
      Particle kaonm(km,tgps.GetDecay(1));
      a = tgps.GetDecay(1);
      decayP.push_back(*a);
      decayParticles.push_back(kaonp);
      decayParticles.push_back(kaonm);
    } else {
      double masses[3] = {pip->Mass(),pip->Mass(),pim->Mass()};
      tgps.SetDecay(in,3,masses);
      w = tgps.Generate();
      Particle pionp1(pip,tgps.GetDecay(0));
      Particle pionp2(pip,tgps.GetDecay(1));
      Particle pionm(pim, tgps.GetDecay(2));
      decayParticles.push_back(pionp1);
      TLorentzVector* a = tgps.GetDecay(0);
      decayP.push_back(*a);
      decayParticles.push_back(pionp2);
      a = tgps.GetDecay(1);
      decayP.push_back(*a);
      decayParticles.push_back(pionm);
      a = tgps.GetDecay(2);
      decayP.push_back(*a);
    }


    vector<vector<TVector3>> reco;
     
    for(int i = 0;i<decayParticles.size();i++){
      for(int j = i+1;j<decayParticles.size();j++){
	if(decayParticles[i].Charge()*decayParticles[j].Charge()<0){
	  bool pidi = PId(decayParticles[i],&decayP[i]);
	  bool pidj = PId(decayParticles[j],&decayP[j]);
	  
	   if ((pidi | pidj)){
	    vector<TVector3> def;
	    def.push_back(decayParticles[i].Vect());
	    def.push_back(decayParticles[j].Vect());
	    reco.push_back(def);
	  }
	}
      }	
    }
    for(int i = 0; i<reco.size();i++){
      Particle finale[2];
      finale[0].SetVectM(reco[i][0],kp->Mass());
      finale[1].SetVectM(reco[i][1],kp->Mass());
      
  auto vect = finale[0]+finale[1];
  hM.Fill(vect.M(),w);
    }
  }
  
   hM.Draw();
  app.Run(true);
  return 0;
}






   
 
