#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TApplication.h>
#include <iostream>

//g++ --std=c++11 -lEG JpsiExp_sol.cpp -o JpsiExp_sol.x `root-config --cflags --glibs`

TRandom3 rnd;
using namespace std;

TDatabasePDG pdg;
TParticlePDG *pr  = pdg.GetParticle("proton");
TParticlePDG *ep  = pdg.GetParticle("e+");
TParticlePDG *em  = pdg.GetParticle("e-");    //code 11
TParticlePDG *pip = pdg.GetParticle("pi+");   //code 211
TParticlePDG *pim = pdg.GetParticle("pi-");
TParticlePDG *jp  = pdg.GetParticle("J/psi");

//creo classe particella costruita da TLorentzVector e da TParticlePDG
class Particle: public TLorentzVector{
public:
  using TLorentzVector::TLorentzVector;
  Particle(TParticlePDG p, TLorentzVector t):TLorentzVector(t),m_code(p.PdgCode()),m_charge(p.Charge()){};
  Particle(TParticlePDG *p, TLorentzVector *t):TLorentzVector(*t),m_code(p->PdgCode()),m_charge(p->Charge()){};
  int Charge() {return m_charge;}
  int PdgCode(){return m_code;}
private:
  int m_code;
  int m_charge;
};


//particle identification trk=traccia. Identifico gli elettroni perche J/psi->ep+em
bool PId(Particle trk){
  double pmis =  0.1;             //prob di scambiare un pi per un e
  double peff =  1.0;             //efficienza nel rivelare un e
  bool id = false;
  if (abs(trk.PdgCode())==211)    //se la traccia è un pi ma sbaglio
    if (rnd.Rndm()<pmis)
      id = true;
  if (abs(trk.PdgCode())==11)     //se la traccia è un e
    if (rnd.Rndm()<peff)
      id = true;
  return id;
}

//risoluzione del rivelatore (prende una particella e aggiunge gli errori di misura)
void applyReso(Particle& trk){
  
  double etheta = 0.015;
  double ephi   = 0.015;
  double ep     = 0.0013;
  
  trk.SetTheta(rnd.Gaus(trk.Theta(),etheta));       //gaussiana attorno al valore dato 
  trk.SetPhi(rnd.Gaus(trk.Phi(),ephi));        

  TVector3    p = trk.Vect();                       //estraggo triimpulso
  double   erel = rnd.Gaus(1,ep*p.Mag());           //ep è propto p    
  p = erel*p;

  //  trk.SetVectM(p,em->Mass());
  trk.SetRho(p.Mag());
  return;
}

vector<Particle>  recoTrks(vector<Particle> vtrk){
  double dang=20/180.*M_PI, theta=40/180.*M_PI;           //accettanza angolare 20 deg
  double pmin=0.1;                                        //accettanza impulso
  int nup=0,ndown=0;
  vector <Particle> sel;

  for (int i=0;i<vtrk.size();i++){                        //ciclo sul vettore di particelle
    if (vtrk[i].P()<pmin) continue;
    if (abs(vtrk[i].Theta()-theta)<dang){                 // accettanza theta
      if (abs(vtrk[i].Phi()-M_PI/2)<dang){                // accettanza phi up
      	if (!PId(vtrk[i])) continue;                      //se non e' ne' pi ne' e continua
      	sel.push_back(vtrk[i]);
      	nup++;                                            //ramo up
      } else if (abs(vtrk[i].Phi()+M_PI/2)<dang){         // accettanza phi down
      	if (!PId(vtrk[i])) continue;
      	sel.push_back(vtrk[i]);
      	ndown++;                                          //ramo down
      }
    }
  }
  //richiedo che ci sia solo una traccia carica per ogni braccio (le due traccie hanno q opposta)
  if (ndown==1 && nup==1){                      
    if (sel[0].Charge()*sel[1].Charge()<0){
      for (int i=0;i<sel.size();i++){
	     applyReso(sel[i]);                  //applico la risoluzione del rivelatore
      }
    } else {
      sel.clear();                         //svuoto il vettore se non rispetto la condizione
    }
  }
  return sel;
}

////////////////////////////////////////////////////////////////////////////////////

int main(){

  double pbeam = 30;              //impulso iniziale del fascio
  double mBe   =  9*0.9315;       //massa del Berillio9 in MeV
  
  TApplication app("app",0,NULL);

  TGenPhaseSpace tgps;                        //genero spazio delle fasi
  rnd.SetSeed(1234578);

  double pjpsi = 1e-4;                        //prob Jpsi
  int nEv      = 1e7;                         //numero eventi

  double w;
  TH1D htheta("htheta","",100,0,M_PI);
  TH1D hphi("hphi","",100,-M_PI,M_PI);
  TH1D hp("hp","",60,0,30);
  TH1D hM("hm","",30,2.5,4.5);
  for (int i=0;i<nEv;i++){

    TLorentzVector in(0,0,pbeam,mBe+sqrt(pbeam*pbeam+pr->Mass()*pr->Mass())); //quadrimpulso tot in
    vector<Particle> decayParticles;

    double prob = rnd.Rndm();

    if (prob<pjpsi){                                  //se ho Jpsi

      double masses[3]={mBe,jp->Mass(),pip->Mass()};
      tgps.SetDecay(in,3,masses);                     //Be9-> Jpsi+pip
      w = tgps.Generate();
      Particle pion(pip,tgps.GetDecay(2));            //salvo il pione nel vettore di prodotti
      decayParticles.push_back(pion);
      
      TLorentzVector* vjp  = tgps.GetDecay(1);        //salvo la cinematica di Jpsi
      double me[2]={ep->Mass(),em->Mass()};
      tgps.SetDecay(*vjp,2,me);                       //Jpsi -> ep+em
      w = w*tgps.Generate();                          //peso finale
      
      Particle ele1(ep,tgps.GetDecay(0));             //salvo il positrone
      Particle ele2(em,tgps.GetDecay(1));             //salvo l'elettrone
      decayParticles.push_back(ele1); 
      decayParticles.push_back(ele2);
      
    } else {                                          //se non ho Jpsi
      if (rnd.Rndm()<(1-pjpsi)){                      //3 pioni

      	double masses[4]={mBe,pip->Mass(),pip->Mass(),pim->Mass()};  //Be9->pip+pip+pim
      	tgps.SetDecay(in,4,masses);
      	w = tgps.Generate();
      	Particle pion1(pip,tgps.GetDecay(1));            //salvo i tre pioni
      	Particle pion2(pip,tgps.GetDecay(2));
      	Particle pion3(pim,tgps.GetDecay(3));
      	decayParticles.push_back(pion1);
      	decayParticles.push_back(pion2);
      	decayParticles.push_back(pion3);
      } else {                                        //5 pioni
      	double masses[6]={mBe,pip->Mass(),pip->Mass(),pip->Mass(),pip->Mass(),pip->Mass()};
      	tgps.SetDecay(in,6,masses);
      	w = tgps.Generate();
      	Particle pion1(pip,tgps.GetDecay(1));
      	Particle pion2(pip,tgps.GetDecay(2));
      	Particle pion3(pim,tgps.GetDecay(3));
      	Particle pion4(pip,tgps.GetDecay(4));
      	Particle pion5(pim,tgps.GetDecay(5));
      	decayParticles.push_back(pion1);
      	decayParticles.push_back(pion2);
      	decayParticles.push_back(pion3);
      	decayParticles.push_back(pion4);
      	decayParticles.push_back(pion5);
      }
    }

    // vector<Particle> reco = rec...
    // for (int i=0;i<reco.size();i++){
    //   htheta.Fill(reco[i].Theta(),w);
    // }
    auto reco = recoTrks(decayParticles);   //identificazione e accettanza
    for (auto part:reco){
      htheta.Fill(part.Theta(),w);          //grafici 
      hphi.Fill(part.Phi(),w);
      hp.Fill(part.P(),w);
    }
    if (reco.size()==2){                              //massa invariante
      reco[0].SetVectM(reco[0].Vect(),em->Mass());
      reco[1].SetVectM(reco[1].Vect(),em->Mass());
      auto cand = reco[0]+reco[1];                    //candidato Jpsi ricostruito
      hM.Fill(cand.M(),w);
    }
  }
  TCanvas c("c","",5,5,700,300);
  c.Divide(4,1);
  c.cd(1);
  htheta.Draw();
  c.cd(2);
  hphi.Draw();
  c.cd(3);
  hp.Draw();
  c.cd(4);
  hM.Draw();
  app.Run(true);
  return 0;
}
