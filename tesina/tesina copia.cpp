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
#include <functional>
#include <string>


//g++ --std=c++11 -lEG tesina.cpp -o tesina.x `root-config --cflags --glibs`


using namespace std;

////GLOBAL VARIABLES

TRandom3 rnd;

  double Z=29;//82;                                  //numero atomico materiale, ora è Cu//Pb
  double A=63.546;//208;      
  double L=10;                                              //distanza in out in cm
  double l=300e-04;                                         //spessore dello strato in cm
  double rho=8.96;//11.29;                             //density g/cm3
  double c=TMath::Ccgs();
  int nEv=1e5; 
  double pD0=0.30; 
  double tres=200e-12;                                  

  double nK_gen=0;
  double nK_reco=0;
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

bool PId(Particle trk){

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

  /*cout << "TOF PI" << setw(30) << "TOF K" << setw(30) << "TOF reco" << endl;
  cout << TOF_p << setw(30) << TOF_k << setw(30) << TOF_trk << endl;*/

  //calcolo la prob che sia pione e quella che sia kaone

 double pK_E=TMath::Gaus(Etrk,E_k,0.5);
 double pK_TOF=TMath::Gaus(TOF_trk,TOF_k,tres);
 double pPI_E=TMath::Gaus(Etrk,E_p,0.5);
 double pPI_TOF=TMath::Gaus(TOF_trk,TOF_p,tres);

  double pK_tot=(pK_TOF*pK_E)/(pK_TOF*pK_E+pPI_TOF*pPI_E);
  double pPI_tot=(pPI_TOF*pPI_E)/(pK_TOF*pK_E+pPI_TOF*pPI_E);
  //pK_tot=pK_E/(pPI_E+pPI_E);
  
  if(abs(trk.PdgCode())==321){//cout << "it is a kaon"<<endl;
  }
  else {//cout << "it is a pion" << endl;
  }

  //cout << "pK|E" << setw(20) << "pK|TOF" << setw(20) << "pPI|E" << setw(20) << "pPI|TOF" << endl;
  //cout << pK_E << setw(20) << pK_TOF << setw(20) << pPI_E << setw(20) << pPI_TOF << endl << endl;
  //cout << "pK_tot= " << pK_tot << endl;
  //cout << "nK=" << nK << endl << endl;
  
  if(pK_tot/pPI_tot>0.9) {
    id=true;
    //cout << "identified as a kaon" << endl;
  }
  else {//cout << "identified as a pion" << endl;
  }


  return id;
}



///////////////////////////////////   APPLICAZIONE    /////////////////////////////////////////
/*
Un collisionatore e+e- a energie tra 3 e 5 GeV produce per la maggior parte coppie di adroni charm
supponiamo per semplicità che produca solo mesoni D*+ (in realtà produce D*+,D+,D0,Ds e barioni con charm)
e supponiamo che questi possano decadere solo in 3 pioni o D0pi con seguente D0 -> k pi
In tutto ci saranno 6 particelle cariche nello stato finale. Cercando le combinazioni +-+ (pi/k/pi)
si individuano i decadimenti del D*+ in D0pi.
Per fare questo
  - si applica l'identificazione per il kaone
  - si plotta la massa inv. D*+
  - si plotta la diff tra la mass inv. D*+ e D0 (questa seconda è più discriminante perché seleziona un piccolo spazio
     delle fasi per il fondo combinatorio
*/

int main(int argc, char **argv){
  bool ID=true;
  if(argc!=2){
    cout << "inserire pId ON/OFF" << endl << endl;
    return 0;
  }
  if(string(argv[1])=="OFF") ID=false;

  rnd.SetSeed(12);
  TApplication app("app",0,NULL);
  TCanvas *can=new TCanvas("Massa Invariante","Massa Invariante",800,800);
  TH1D *histo= new TH1D("PID ON","PID ON",100,0,0);
  histo->GetXaxis()->SetTitle("sqrt(s)");
  histo->Draw();
  TH1D *histoOFF=new TH1D("PID OFF","PID OFF",100,0,0);
  histoOFF->GetXaxis()->SetTitle("sqrt(s)");
  histoOFF->Draw();
  gPad->Modified(); gPad->Update();
  //app.Run(true);

  TGenPhaseSpace tgps;                        //genero spazio delle fasi
                             
  double w;
  /*cout << "masse" << endl;
  cout <<endl << setw(10) << "D*p" << setw(10) << "D0" <<setw(10) << "D*m"<< endl;
  cout <<setw(10) << Dp->Mass() << setw(10) << D0->Mass() <<setw(10) << Dm->Mass()<< endl<<endl;
*/
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

   // cout << v_Dp->X() << setw(10) << v_Dp->Y() << setw(10) << v_Dp->Z() << setw(10) << v_Dp->T() << endl;
    
    //decadimento D*p
    if(rnd.Rndm()<pD0){
      //cout << "D*p->D0+pip"<<endl;
      double masses[2]={D0->Mass(),pip->Mass()};  //D*p->D0+pip
      tgps.SetDecay(*v_Dp,2,masses);              
      w=w*tgps.Generate();                        //cout << "genero" << endl;
      Particle pion1(pip,tgps.GetDecay(1));       
      products.push_back(pion1);                   //aggiungo pione 1 ai prodotti 
      //cout << "pushback" << endl;
      TLorentzVector* vD0=tgps.GetDecay(0);
      //cout << "getdecay" << endl;
      masses[0]=km->Mass();                        //cambio le masse dei prodotti
      masses[1]=pip->Mass();       
      //decadimento D0 in k e pi
      //cout << "d0->k pi"<<endl;
      tgps.SetDecay(*vD0,2,masses);
      w=w*tgps.Generate();            //cout << "genero" << endl;
    
      Particle pion2(pip,tgps.GetDecay(1));
      Particle kaon(km,tgps.GetDecay(0));
      products.push_back(pion2);                  //aggiungo pione 2 ai prodotti
      products.push_back(kaon);                   //aggiungo kaone ai prodotti
      nTest++;
      //cout << pion1.P() << setw(10) << pion2.P() << setw(10) << kaon.P() << endl;

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

     // cout << pion1.P() << setw(10) << pion2.P() << setw(10) << pion3.P() << endl;

    }

    //ho tutti i prodotti salvati nel vector<Particle> products
    //devo identificarli attraverso BB e TOF 
    //cout << "nProd="<<products.size() << endl;      //dovrebbe essere 3 
    //mi interessano i kaoni per ricostruire D0
    //cerco tutte le coppie +- in cui uno é k e plotto la massa invariante
    if(ID){
      for(int i = 0;i<products.size();++i){
        //conteggio per verifica
        for(int j = i+1;j<products.size();++j){
          if(products[i].Charge()*products[j].Charge()<0){      //coppie (+,-)
            bool pId_i = PId(products[i]);                      //true == km, false == pip
            bool pId_j = PId(products[j]);   
            if ((pId_i==true && pId_j==false)|(pId_i==false && pId_j==true)){ //se solo uno é un k
              TVector3 v1=products[i].Vect();
              TVector3 v2=products[j].Vect();
              Particle p1; p1.SetVectM(v1,km->Mass());
              Particle p2; p2.SetVectM(v2,pip->Mass());
              ++nK_reco;
              if(abs(products[i].PdgCode())==321|abs(products[j].PdgCode())==321){
              ++nK_gen;
              }
              //se é pi ma l'ho identificato come k
              if((abs(products[i].PdgCode())==221 && pId_i==true) | (abs(products[j].PdgCode())==221 && pId_j==true)) {
                ++nMis;
                }
              //se é un kaone ma lo identifico pi
              if((abs(products[i].PdgCode())==321 && pId_i==false)|(abs(products[j].PdgCode())==321 && pId_j==false)){
                ++nLost;
              }

              auto Ptot=p1+p2;
              double s=Ptot.M();
            //  cout << "s=" << s << endl;
            //  cout << "w=" << w << endl << endl;
              histo->Fill(s,w);
              /*histo->Draw();
              gPad->Modified();
              gPad->Update();*/
            }
          }
        } 
      }
    } 
    else {
      for(int i = 0;i<products.size();++i){
        for(int j = i+1;j<products.size();++j){
          if(products[i].Charge()*products[j].Charge()<0){      //coppie (+,-)
            bool pId_i = PId(products[i]);                      //true == km, false == pip
            bool pId_j = PId(products[j]); 
            TVector3 v1=products[i].Vect();
            TVector3 v2=products[j].Vect();
            Particle p1; p1.SetVectM(v1,km->Mass());
            Particle p2; p2.SetVectM(v2,pip->Mass());
            ++nK_reco;
            if(abs(products[i].PdgCode())==321|abs(products[j].PdgCode())==321){
            ++nK_gen;
            }
              //se é pi ma l'ho identificato come k
              if((abs(products[i].PdgCode())==221 && pId_i==true) | (abs(products[j].PdgCode())==221 && pId_j==true)) {
                ++nMis;
                }
              //se é un kaone ma lo identifico pi
              if((abs(products[i].PdgCode())==321 && pId_i==false)|(abs(products[j].PdgCode())==321 && pId_j==false)){
                ++nLost;
              }

            auto Ptot=p1+p2;
            double s=Ptot.M();
          //  cout << "s=" << s << endl;
          //  cout << "w=" << w << endl << endl;
            histo->Fill(s,w);
            /*histo->Draw();
            gPad->Modified();
            gPad->Update();*/
          }
        } 
      }
    }
  }

  //cout << "arrivo qui" << endl;
  histo->Draw();
  cout << "nK_gen=" << nK_gen << endl << endl;
  cout << "nK_reco=" << nK_reco << endl;
  cout << "nMis="   << nMis << endl<<endl;
  cout << "nLost=" << nLost << endl;
  cout << "nTest=" << nTest << endl;
  app.Run(true);

  return 0;
}
