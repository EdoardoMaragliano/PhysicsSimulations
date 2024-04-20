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
#include <THStack.h>
#include <iomanip>
#include <TLatex.h>
#include <TH1.h>
#include <TStyle.h>
#include <functional>
#include <string>
#include <time.h>


//g++ --std=c++11 -lEG tesina.cpp -o tesina.x `root-config --cflags --glibs`


using namespace std;

////GLOBAL VARIABLES

TRandom3 rnd;

  double Z=14;                                  //numero atomico materiale, ora è Cu//Pb
  double A=28.08;      
  double L=50;                                      //distanza in out in cm
  double l=300e-04;                                 //spessore dello strato in cm
  double rho=2.33;                                  //density g/cm3
  double c=TMath::Ccgs();
  int nEv=1e6; 
  double pD0=0.1; 
  double tres=100e-12;                                //100ps                              

  double nK_gen=0;
  double nK_reco=0;
  double nPI_reco=0;
  double nMis=0;
  double nLost=0;
  double nTest=0;

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

//creo classe particella costruita da TLorentzVector e da TParticlePDG
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

double TimeOfFlight(Particle part){
  return L/TMath::Ccgs()*sqrt(1+pow(part.Mass()/part.P(),2)); 
}

double BetheBloch(Particle p, double Z, double A ){
  //NB masse DPG in GeV
  double K=0.307;                       //Mev g^-1 cm^2    //da Taiuti
  double I=10;
  if(Z>=13){
   I=9.76*Z+58.8*pow(Z,-0.19);          //eV per Z>13 da Taiuti
  } else{
   I=12*Z+7; //eV per Z<13
  }
  //K=K*1e-3; //in GeV g^-1 cm^2
  I=I*1e-6; //in MeV
  double argLog=2*0.511*pow(TMath::Ccgs()*p.Beta()*p.Gamma(),2)/I;
  double delta=0; //non trovo espressione per delta 
  return rho*K*Z/A*pow(p.Charge()/p.Beta(),2)*(log(argLog)-pow(p.Beta(),2)-delta/2)*l;
}
//In MeV

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

  trk.SetRho(p.Mag());
  return;
}

bool PId(Particle trk,bool flagTOF,bool flagEn){

  bool id = false;
  TLorentzVector *v=new TLorentzVector(trk.X(),trk.Y(),trk.Z(),trk.T());
  Particle pi(pip, v);
  Particle kaon(kp,v);

  double Etrk=BetheBloch(trk,Z,A);   //valore nominale dell'energia persa

  //poiché il mio esperimento gaussianizza i dati troncando la Landau al 55%, assumo a questo punto
  //che l'energia persa sia distribuita secondo una Gaussiana centrata in BB e con StDev 
  //ottenuta dalla distribuzione delle medie troncate studiata in precedenza

  Etrk = rnd.Gaus(Etrk,0.5);                //questo é il valore ricostruito
  double TOF_trk=TimeOfFlight(trk);

  TOF_trk=rnd.Gaus(TOF_trk,tres);

  //calcolo i valori attesi per pione e kaone
  double E_p = BetheBloch(pi,Z,A);
  double E_k = BetheBloch(kaon,Z,A);
  double TOF_p  = TimeOfFlight(pi);
  double TOF_k  = TimeOfFlight(kaon);

  //calcolo la prob che sia pione e quella che sia kaone

 double pK_E=TMath::Gaus(Etrk,E_k,0.5);
 double pPI_E=TMath::Gaus(Etrk,E_p,0.5);
 double pK_TOF=TMath::Gaus(TOF_trk,TOF_k,tres);
 double pPI_TOF=TMath::Gaus(TOF_trk,TOF_p,tres);

  double pK_tot=(pK_TOF*pK_E)/(pK_TOF*pK_E+pPI_TOF*pPI_E);
  double pPI_tot=(pPI_TOF*pPI_E)/(pK_TOF*pK_E+pPI_TOF*pPI_E);

  if(flagTOF==true && flagEn==false){   //solo tof
    pK_tot=pK_TOF/(pK_TOF+pPI_TOF);
    pPI_tot=pPI_TOF/(pK_TOF+pPI_TOF);
    //cout << setw(10) << pK_tot << setw(10) << pPI_tot << setw(10) << pK_tot+pPI_tot << endl;
  }
  if(flagEn==true && flagTOF==false){   //solo energia
    //cout << "only energy" << endl;
    pK_tot=pK_E/(pK_E+pPI_E);
    pPI_tot=pPI_E/(pK_E+pPI_E);
  }

  /*
  cout << endl <<endl <<endl;
  if(abs(trk.PdgCode())==321){
  cout << "it is a kaon" <<endl;}
  
  else {cout << "it is a pion" << endl;
  }
  cout << "TOF PI" << setw(30) << "TOF K" << setw(30) << "TOF reco" << endl;
  cout << TOF_p << setw(30) << TOF_k << setw(30) << TOF_trk << endl;

  cout << "pK|E" << setw(20) << "pK|TOF" << setw(20) << "pPI|E" << setw(20) << "pPI|TOF" << endl;
  cout << pK_E << setw(20) << pK_TOF << setw(20) << pPI_E << setw(20) << pPI_TOF << endl << endl;
  cout << "pK_tot= " << pK_tot << endl;
  cout << "pPI_tot=" << pPI_tot<<endl;*/
  
  if(pK_tot/pPI_tot>2) {
    id=true;
    //cout << "identified as a kaon" << endl;
    if(abs(trk.PdgCode())==abs(kp->PdgCode())) ++nK_reco;      //é k e scelgo k
    if(abs(trk.PdgCode())==abs(pip->PdgCode())) ++nMis;         //é pi ma scelgo k
  }
  else {
    //cout << "identified as a pion" << endl;
    if(abs(trk.PdgCode())==abs(kp->PdgCode())) {++nLost;}    //ho perso un k scambiandolo per pi
    if(abs(trk.PdgCode())==abs(pip->PdgCode())) {++nPI_reco;}
  }
  return id;
}



///////////////////////////////////   APPLICAZIONE    /////////////////////////////////////////


int main(int argc, char **argv){
  cout << endl << endl <<endl << setw(20)<< "DEFAULT SETTINGS " << endl;
  cout << "nEv=" << nEv << endl;
  cout << "pD0=" << pD0 << endl;
  cout << "particle id by dEdx and TOF" << endl << endl;
  cout << "to change settings execute as:" << endl << endl;
  cout << "./tesina.x nEv pD0 TOF(on/off) dEdx(on/off)" << endl;
  cout << "WARNING: either TOF or dEdx MUST be on"  << endl;
  cout << "a graph without PID is automatically generated" << endl;
  cout << "******************************************************" << endl << endl;

  bool flagTOF=true;
  bool flagEn=true;

  if(argc!=1 & argc!=5){cout << "wrong number of main parameters, abort" << endl;
    return -1;
  }
  
  if(argc==5) {
    cout << endl << endl <<endl << setw(20)<< "CUSTOM SETTINGS " << endl;
    nEv=atoi(argv[1]);
    pD0=atof(argv[2]);
    cout << "nEv=" << nEv << endl;
    cout << "pD0=" << pD0 << endl;
    if(string(argv[3])=="off"){flagTOF=false;
     cout << "Pid by TOF is off" << endl;}
    if(string(argv[3])=="on"){
     cout << "Pid by TOF is on" << endl;}

    if(string(argv[4])=="off"){flagEn=false;
      cout << "Pid by dEdx is off" << endl;}
    if(string(argv[4])=="on"){
     cout << "Pid by dEdx is on" << endl;}

    if (string(argv[3])=="off" & string(argv[4])=="off"){
      cout << "WARNING: either TOF or dEdx MUST be on"  << endl;
      cout << "a graph without PID is automatically generated" << endl << endl;
      cout << "proceeding with both TOF and dEdx on " << endl;

    }
  } 

  cout << "starting simulation" << endl;
  rnd.SetSeed(time(0));
  TApplication app("app",0,NULL);
  TCanvas *can=new TCanvas("Massa Invariante","Massa Invariante",800,1000);
  TH1D *histo= new TH1D("PID ON","PID ON",100,0,2);
  histo->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
  TH1D *histoOFF=new TH1D("PID OFF","PID OFF",100,0,2);
  histoOFF->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
  gPad->Modified(); gPad->Update();


  TGenPhaseSpace tgps;                        //genero spazio delle fasi
                             
  double w;
  cout <<endl << endl<< setw(10)<< "masse" << endl;
  cout <<endl << setw(10) << "D*p" << setw(10) << "D0" <<setw(10) << "D*m"<< endl;
  cout <<setw(10) << Dp->Mass() << setw(10) << D0->Mass() <<setw(10) << Dm->Mass()<< endl<<endl;



  for(int i=0;i<nEv;++i){
    //URTO e- e+ nel CM (fasci opposti=>P=0) Edisponibile=5Gev
   
    vector <Particle> products;

    //decadimento Jpsi
    TLorentzVector in(0.,0.,0.,4.040);       //ferma nel CM
    //cout <<"Jpsi -> D*p + D*m"<<endl;
    double masses[2]={Dp->Mass(),Dm->Mass()};
    tgps.SetDecay(in,2,masses);                       
    w=tgps.Generate();                              
    TLorentzVector *v_Dp=tgps.GetDecay(0);          //cinematica di D*p
    TLorentzVector *v_Dm=tgps.GetDecay(1);          //cinematica di D*m
    
    //decadimento D*p
    if(rnd.Rndm()<pD0){
      //cout << "D*p->D0+pip"<<endl;
      double masses[2]={D0->Mass(),pip->Mass()};  
      tgps.SetDecay(*v_Dp,2,masses);              
      w=w*tgps.Generate();                        
      Particle pion1(pip,tgps.GetDecay(1));       
      products.push_back(pion1);                   //aggiungo pione 1 ai prodotti 
      TLorentzVector* vD0=tgps.GetDecay(0);
      //cout << "getdecay" << endl;
      masses[0]=km->Mass();                        //cambio le masse dei prodotti
      masses[1]=pip->Mass();       
      //decadimento D0 in k e pi
      //cout << "d0->k pi"<<endl;
      tgps.SetDecay(*vD0,2,masses);
      w=w*tgps.Generate();            
    
      Particle pion2(pip,tgps.GetDecay(1));
      Particle kaon(km,tgps.GetDecay(0));
      products.push_back(pion2);                  //aggiungo pione 2 ai prodotti
      products.push_back(kaon);                   //aggiungo kaone ai prodotti
      ++nK_gen;

    } else {  

      //decadimento Dstar in 3 pi
      //cout << "D*p->pip+pip+pim" << endl;
      double masses[3]={pip->Mass(),pip->Mass(),pim->Mass()};
      tgps.SetDecay(*v_Dp,3,masses);
      w=w*tgps.Generate();
      Particle pion1(pip,tgps.GetDecay(0));
      Particle pion2(pip,tgps.GetDecay(1));
      Particle pion3(pim,tgps.GetDecay(2));
      products.push_back(pion1);
      products.push_back(pion2);
      products.push_back(pion3);
    }
    //applico risoluzione in impulso ai prodotti
    for(auto trk:products)
      applyReso(trk);

    //ho tutti i prodotti salvati nel vector<Particle> products
    //devo identificarli attraverso BB e TOF 
    //cerco tutte le coppie +- in cui uno é k e plotto la massa invariante

    for(int i = 0;i<products.size();++i){
      for(int j = i+1;j<products.size();++j){
        if(products[i].Charge()*products[j].Charge()<0){    //coppie (+,-)

          bool pId_i = PId(products[i], flagTOF, flagEn);                 
          bool pId_j = PId(products[j], flagTOF, flagEn);   
          if ((pId_i==true && pId_j==false)|(pId_i==false && pId_j==true)){ //se solo uno é un k
            TVector3 v1=products[i].Vect();
            TVector3 v2=products[j].Vect();
            Particle p1; p1.SetVectM(v1,km->Mass());
            Particle p2; p2.SetVectM(v2,pip->Mass());

            auto Ptot=p1+p2;
            double s=Ptot.M();
            histo->Fill(s,w);
          }
        }
      } 
    } 

    for(int i = 0;i<products.size();++i){
      for(int j = i+1;j<products.size();++j){
        //if(products[i].Charge()*products[j].Charge()<0){      //coppie (+,-)
          TVector3 v1=products[i].Vect();
          TVector3 v2=products[j].Vect();
          Particle p1; p1.SetVectM(v1,km->Mass());
          Particle p2; p2.SetVectM(v2,pip->Mass());
          
          auto Ptot=p1+p2;
          double s=Ptot.M();
          histoOFF->Fill(s,w);
        //}
      } 
    }
  }

  gStyle->SetOptStat("ne");
  TF1 f("f","[5]*TMath::Gaus(x,[0],[1])+[2]*x**2+[3]*x+[4]",1.7,2);
  TF1 g("g","[5]*TMath::Gaus(x,[0],[1])+[2]*x**2+[3]*x+[4]",1.7,2);
  TF1 psi("psi","[0]*x**2+[1]*x+[2]",1.2,1.8);
  f.SetParameter(0,D0->Mass());
  f.SetParameter(1,0.01);
  f.SetParameter(2,-14450);
  f.SetParameter(3,50000);
  f.SetParameter(4,-500);
  f.SetParameter(5,1000);
  g.SetParameter(0,D0->Mass());
  g.SetParameter(1,0.01);

  psi.SetParameter(0,-1.595e04);
  psi.SetParameter(1,3.973e04);
  psi.SetParameter(2,-1.355e04);

  if(pD0<0.03){

  f.SetParameter(0,D0->Mass());
  f.SetParameter(1,0.01);
  f.SetParameter(2,-14450);
  f.SetParameter(3,50000);
  f.SetParameter(4,-500);
  f.SetParameter(5,1000);
  g.SetParameter(0,D0->Mass());
  g.SetParameter(1,0.01);

  psi.SetParameter(0,-1.595e04);
  psi.SetParameter(1,3.973e04);
  psi.SetParameter(2,-1.355e04);
  }

  //cout << "arrivo qui" << endl;
  can->Divide(2);
  can->cd(1);
  //histo->Fit("f","LR");
  histo->Draw();
  
  can->cd(2);
  histoOFF->Fit("psi","LR");
  g.SetParameter(2,psi.GetParameter(0));
  g.SetParameter(3,psi.GetParameter(1));
  g.SetParameter(4,psi.GetParameter(2));
  g.SetParameter(5,3e04);
  //histoOFF->Fit("g","LR");
  histoOFF->Draw();

  cout << endl << endl << endl;
  cout << "nK_gen=" << nK_gen << endl;
 // cout << "nEv=" << nEv<<endl;
  cout << "nK_reco=" << nK_reco << endl;
  cout << "nMis="   << nMis << endl;
  cout << "nLost=" << nLost << endl;
  cout << "nPI_reco="<<nPI_reco<<endl;
  //cout << "nTest=" << nTest << endl

  cout << endl << setw(30)<< "signal/noise ratio" << endl <<endl;
  cout << setw(20) << "PID ON "<< setw(20) << "PID OFF" << endl ;
  /*cout << setw(20) << f.Eval(f.GetParameter(0))/f.Eval(f.GetParameter(0)-3*f.GetParameter(1)) 
  <<setw(20)<< g.Eval(g.GetParameter(0))/g.Eval(g.GetParameter(0)-3*g.GetParameter(1)) << endl;*/
  cout << setw(20) << f.Eval(f.GetParameter(0))/f.Eval(1.8) 
  <<setw(20)<< g.Eval(g.GetParameter(0))/g.Eval(1.8) << endl;
  app.Run(true);

  return 0;
}
