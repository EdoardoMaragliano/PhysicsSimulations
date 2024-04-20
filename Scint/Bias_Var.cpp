#include <TH1D.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TArc.h>
#include <TApplication.h>
#include <iostream>
#include <TGraph.h>
#include <TF1.h>
#include <TStopwatch.h>
#include <TFile.h>

// g++ Bias_var.cpp -o Bias_var.x `root-config --cflags --glibs`

using namespace std;



double me   = 0.511;   
double n    = 3.67*1e6/(149.89)*6.022e23;

double r    = 0.02;
double dz   = 0.025;
double dist = 0.10;
double E    = 1.27;
int    nev  = 1e6;
int    ngr  = 100;

TRandom3 rnd;
TCanvas   *c;

//Sezione d'urto effetto fotoelettrico
double fPE(double *x, double *p){     
  double ub=1e-28;
  double sig;
  double e = *x;
  if (e<0.033){
    sig = 0.06 * pow(e,-2.89);
  } else {
    sig = 0.92 * pow(e,-2.65);
  }
  return sig*ub;
}

//Sezione d'urto effetto Compton
double fComp(double *x, double *p){
  double ub=1e-28;
  double sig;
  double e = *x;
  sig    = 18*TMath::Exp(-0.3*e)+24*TMath::Exp(-3*e);
  return   sig*ub;
}


// Disegno della geometria
void DrawGeom(double r, double dz){
  double framew = r*1.5;
  // Disegno il volume dello scintillatore
  c = new TCanvas("c","",10,10,900,320);
  c->Divide(3,1);
  //Vista X-Y
  c->cd(1);
  gPad->Range(-framew,-framew,framew,framew);
  TArc *arc = new TArc;
  arc->SetFillColor(11);
  arc->DrawArc(0.,0.,r);
  // Vista Z-Y
  c->cd(2);
  gPad->Range(-framew/2,-framew,framew+framew/2,framew);
  TGraph *grp = new TGraph;
  grp->SetPoint(0,0 ,-r);
  grp->SetPoint(1,dz,-r);
  grp->SetPoint(2,dz,r);
  grp->SetPoint(3,0 ,r);
  grp->SetFillColor(11);
  grp->Draw("F");
}

double Edep(TVector3 x0, TVector3 d0, double E, TGraph*& grxy, TGraph*& grzy, bool fillgraph){

  double Ed=0;
  
  do {
    TVector3 x,d;
    
    // calcolo sezioni d'urto
    double sigma_PE   = fPE(&E,0);      
    double sigma_Comp = fComp(&E,0);
    double sigma_tot  = sigma_PE+sigma_Comp;
    double lambda     = 1/(n*sigma_tot);         

    // calcolo libero cammino
    double s = -lambda*log(rnd.Rndm());
    x  = x0 + s*d0;
    x0 = x;

    // verifico che non sia uscito dal rivelatore
    double Dxy = sqrt(pow(x.X(),2)+pow(x.Y(),2));
    if ( Dxy > r || x.Z() < 0 || x.Z() > dz)
      break;

    //Grafica
    if (fillgraph){
      grxy->SetPoint(grxy->GetN(),x.X(),x.Y());
      grzy->SetPoint(grzy->GetN(),x.Z(),x.Y());
    }
    
    // decido se la prossima interazione e' Compton o Fotoelettrico
    double PE_ov_tot = sigma_PE/(sigma_PE+sigma_Comp);
    
    if(rnd.Rndm()<PE_ov_tot){   //FOTOELETTRICO
      Ed+=E;
      E=0;
    } else{                     //COMPTON
      double b = E/me;
      double ntheta = acos(1 - (1/b)*( pow(1+2*b,rnd.Rndm()) - 1 )   );
      double nphi   = 2*M_PI*rnd.Rndm();
      double Eprim  = E/(1+b*(1-cos(ntheta)));    //energia fotone scatterato compton  
      Ed += (E-Eprim);                            //aggiorno E depositata
      E   = Eprim;                                //aggiorno E fotone
      if(d0.Z()==1){
      	d(0) = sin(ntheta)*cos(nphi);
      	d(1) = sin(ntheta)*sin(nphi);
      	d(2) = cos(ntheta);
      } else{
      	double sq = sqrt(1 - pow(d0.Z(),2));
      	d(0) = d0.X()*cos(ntheta)+sin(ntheta)/sq*(d0.X()*d0.Z()*cos(nphi)-d0.Y()*sin(nphi));
      	d(1) = d0.Y()*cos(ntheta)+sin(ntheta)/sq*(d0.Y()*d0.Z()*cos(nphi)+d0.X()*sin(nphi));
      	d(2) = d0.Z()*cos(ntheta)-sin(ntheta)*sq*cos(nphi);
      }
      d0=d;
    }
  } while (E!=0);

  return Ed;
}


void Bias_Var(string modus){

  rnd.SetSeed(time(NULL));

  // Disegno geometria e istogramma
  DrawGeom(r,dz);
  TStopwatch tstop;
  int naccept=0;

  TH1D *hE = new TH1D("hEdep","",100,0,2);
  c->cd(3);
  hE->Draw();
  double coslim = dist/sqrt(dist*dist+r*r);
  
  // Itero sui fotoni
  for (int iev=0;iev<nev;iev++){
    
    double   th[3];   //theta 
    double   ph[3];   //phi
    TVector3 u[3];    //direzioni
    TVector3 x0,d;
    double   e[3]={1.275,0.511,0.511};    //energie
    int biasmask[3]={0,0,0};              //per riduzione varianza
    double w=1;                           //peso iniziale

    if(modus=="detsimple"){
      x0=TVector3(0,0,1);
      d=TVector3(0,0,1);
    } else if(modus=="det"){
      th[0] = acos(2*rnd.Rndm()-1);
      th[1] = acos(2*rnd.Rndm()-1);
      ph[0] = 2*M_PI*rnd.Rndm();
      ph[1] = 2*M_PI*rnd.Rndm();
      for (int i=0;i<2;i++){
        u[i].SetX(sin(th[i])*cos(ph[i]));
        u[i].SetY(sin(th[i])*sin(ph[i]));
        u[i].SetZ(cos(th[i]));
      }
      u[2] = -u[1];                     //fotoni back to back
    } else if(modus=="bias"){
      double pgeo=(1-coslim)/2;         //probabilità geometrica (sono nel giusto angolo solido?)
      double p1=pgeo*(1-pgeo);          //fotone 1 nel rilevatore
      double p2=2*pgeo*(1-pgeo);        //fotone 2 nel rilevatore
      double p12=2*pgeo*pgeo;           //entrambi nel rilevatore
      w=p1+p2+p12;
      double test=rnd.Rndm()*w;
      if (test<p1){
        biasmask[0]=1;                  //solo il primo
      } else if(test<p1+p2){
        biasmask[1]=1;                  //solo il secondo
      } else {
        biasmask[0]=1;
        biasmask[1]=1;                  //entrambi
      }
    }
    
    double Ed=0;
    
    for (int k=0;k<3;k++){
      if(modus=="detsimple"){                         //sorgente non isotropa
        double theta=acos(2*rnd.Rndm()-1);
        if(cos(theta)<coslim) continue;
        naccept++;
      } else if(modus=="det"){                        //sorgente isotropa
      if (u[k].Z()<coslim) continue;
      //th < th_lim => cos(th)> cos(th_lim)
      //perciò se cos(th)<cos(th_lim) non vado avanti
      //e passo alle iterazioni successive (continue)
      
      double   R = dist/u[k].Z();
      x0=TVector3(R*u[k].X(),R*u[k].Y(),0);
      d=u[k];
      } else if(modus=="bias"){
        if(biasmask[k]==0) continue;  //se biasmask è 0 salto all'iter succ
        naccept++;
        double th=acos((1-coslim)*rnd.Rndm()+coslim);
        double ph=2*M_PI*rnd.Rndm();
        d=TVector3(sin(th)*cos(ph),sin(th)*sin(ph),cos(th));
        double R=dist/d.Z();
        x0=TVector3(R*d.X(),R*d.Y(),0);
      }

      TGraph *grxy,*grzy;
      bool fillgraph=false;
      if (iev<=ngr){
      	grxy = new TGraph();
      	grxy->SetMarkerColor(iev+1);
      	grxy->SetMarkerStyle(20);
      	grxy->SetMarkerSize(0.40);
      	grzy = new TGraph(*grxy);
      	grxy->SetPoint(grzy->GetN(),x0.X(),x0.Y());
      	grzy->SetPoint(grzy->GetN(),x0.Z(),x0.Y());
      	fillgraph = true;
      }

      Ed += Edep(x0,d,e[k],grxy,grzy,fillgraph);
    
      // Partendo dall'energia depositata estraggo l'energia misurata.
      // assumo un errore relativo sull'energia di
      // (sigma(E)/E)^2 = somma in quad(0.02 + 0.04/sqrt(E))
    
      if (fillgraph){
      	c->cd(1);
      	grxy->Draw("LP");
      	c->cd(2);
      	grzy->Draw("LP");
      }
    }
    if (Ed!=0){
      double sigE = Ed*sqrt(pow(0.02,2)+pow(0.04,2)/Ed);
      Ed = rnd.Gaus(Ed,sigE);
      hE->Fill(Ed);
    }
  }
  double sig=hE->GetBinContent(50);   //prendo un bin stabile (una spalla del grafico)
  double err=hE->GetBinError(50);     //prendo il suo errore
  double FOM=1/(tstop.CpuTime()*pow(err/sig,2));
  cout <<"FOM= "<<FOM<<endl;

}

#ifndef __CINT__
int main(int argc, char **argv){
  if(argc!=2){
    cout << "devi selezionare una delle seguenti" << endl;
    cout << "detsimple"<< endl << "det" << endl <<"bias" << endl;
    return 1;
  }
  TApplication app("app",0,NULL);
  Bias_Var(string(argv[1]));
  gPad->Update();
  app.Run(true);
}
#endif
