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


// g++ SimNeuVar.cpp -o SimNeuVar.x `root-config --cflags --glibs`
using namespace std;


double A     = 10;
double n     = 2460/A*6.022e23;

double  dz    = 1;           //spessore rivelatore
double  detr  = 0.05;         //raggio del rivelatore
TVector3 det(0.0,0.0,0.5);    //posizione del rivelatore
  
double E     = 0.1;
int    nneu  = 1e5;
int    ngr   = 100;

TRandom3 rnd;
TCanvas   *c;

double sige = 18e-28;  //sezione d'urto elastica


// Disegno della geometria
void DrawGeom(){
  double framew = dz;
  // Disegno il volume dello scintillatore
  c = new TCanvas("c","",10,10,900,320);
  c->Divide(3,1);
  // Vista X-Y
  c->cd(1);
  gPad->Range(-framew,-framew,framew,framew);
  TGraph *grp1 = new TGraph;
  grp1->SetPoint(0,-framew ,-framew);
  grp1->SetPoint(1, framew ,-framew);
  grp1->SetPoint(2, framew , framew);
  grp1->SetPoint(3,-framew , framew);
  grp1->SetFillColor(11);
  grp1->Draw("F");
  TArc *arc1 = new TArc;
  arc1->SetFillColor(0);
  arc1->DrawArc(det[0],det[1],detr);
  // Vista Z-Y
  c->cd(2);
  gPad->Range(-framew,-framew,framew,framew);
  TGraph *grp2 = new TGraph;
  grp2->SetPoint(0,0,-framew);
  grp2->SetPoint(1,framew ,-framew);
  grp2->SetPoint(2, framew , framew);
  grp2->SetPoint(3,0, framew);
  grp2->SetFillColor(11);
  grp2->Draw("F");
  TArc *arc2 = new TArc;
  arc2->SetFillColor(0);
  arc2->DrawArc(det[2],det[1],detr);
 
}

bool Count(TVector3 x0, TVector3 d0, double& E,double& w, TGraph*& grxy, TGraph*& grzy, bool fillgraph){

  bool cnt=false;
  do {
    
    // Calcolo flusso
    double s_intersec = -d0*(x0-det);         //x0-det posizione rivelatore
    TVector3 intersec = x0 + s_intersec*d0; 
    if ((intersec-det).Mag()<detr){
      cnt=true;   //se il neutrone interseca il rivelatore
      //break;  
    }
    
    TVector3 x,d;
    double lambda     = 1/(n*sige);         //dovrei in teoria usare sigmaTot, ma qui no per semplicita'

    // calcolo libero cammino
    double rho0=0.25;                             
    double rho=rho0*(det-x0).Unit()*d0;           //allungo il libero cammino verso il rivelatore
    double s = -1/lambda*log(rnd.Rndm());       
    double w=w/(1-rho)*exp(-lambda*rho*s);

    x  = x0 + s*d0;
    x0 = x;

    // verifico che non sia uscito dal rivelatore
    if ( x.Z() < 0 || x.Z() > dz){
      break;
    }
    
    //Grafica
    if (fillgraph){
      grxy->SetPoint(grxy->GetN(),x.X(),x.Y());
      grzy->SetPoint(grzy->GetN(),x.Z(),x.Y());
    }

    double pabs= 0.30;  //sez.urto assorbim su totale

    /*if(rnd.Rndm()<pabs){
      E=0;
      break; //assorbo il neutrone
    }
    */
    
    //roulette russa
    w=w*(1-pabs);
    double pr=0.10;       //pr arbitrario, dipende dalla fisica
    if(rnd.Rndm()<pr){    //se il peso è troppo piccolo
      break;              //butto via
    } else{               //aumento il peso alle altre per non cambiare la normalizz
      w=w/(1-pr);
    }
    
    
    double thetaCM = acos(2*rnd.Rndm()-1);
    double theta   = acos( (1+A*cos(thetaCM))/sqrt(A*A+1+2*A*cos(thetaCM)) );
    double phi     = 2*TMath::Pi()*rnd.Rndm();
    E              = E/pow((A+1),2)*pow(cos(theta)+sqrt(A*A-sin(theta)*sin(theta)),2);
    if (d0.Z()==1){
      d.SetX(sin(theta)*cos(phi));
      d.SetY(sin(theta)*sin(phi));
      d.SetZ(cos(theta));
    } else if (d0.Z()==-1){
      d.SetX(-sin(theta)*cos(phi));
      d.SetY(-sin(theta)*sin(phi));
      d.SetZ(-cos(theta));
    } else {
      d.SetX(d0.X()*cos(theta)+sin(theta)/sqrt(1-pow(d0.Z(),2))*(d0.X()*d0.Z()*cos(phi)-d0.Y()*sin(phi)));
      d.SetY(d0.Y()*cos(theta)+sin(theta)/sqrt(1-pow(d0.Z(),2))*(d0.Y()*d0.Z()*cos(phi)+d0.X()*sin(phi)));
      d.SetZ(d0.Z()*cos(theta)-sqrt(1-pow(d0.Z(),2))*sin(theta)*cos(phi));
    }
    d0 = d;

    
  } while (E!=0);
		     
  return cnt;
}


void SimNeuVar(){

  rnd.SetSeed(time(NULL));

  // Disegno geometria e istogramma
  DrawGeom();
  TH1D *hE = new TH1D("Edep","",100,0,0.2);
  TH1D hc("hc","",1,0,1);
  c->cd(3);
  hE->Draw();
  TStopwatch t;

  double rflux = 0.05;    //raggio rivelatore
  int ndet     = 0;
  cout << "flux " << nneu/(rflux*rflux*TMath::Pi()) << endl;
  // Itero sui neutroni  
  for (int in=0; in<nneu; in++){
    
    TVector3 d0(0,0,1);                      //neutrone lungo z
    double  r02 = rnd.Rndm()*pow(rflux,2); 
    double  r0  = sqrt(r02);
    double  th0 = rnd.Rndm()*2*TMath::Pi();
    TVector3 x0(r0*cos(th0),r0*sin(th0),0);
    
    TGraph *grxy,*grzy;
    bool fillgraph=false;
    if (in<=ngr){
      
      grxy = new TGraph();
      grxy->SetMarkerColor(in+1);
      grxy->SetMarkerStyle(20);
      grxy->SetMarkerSize(0.40);

      grzy = new TGraph(*grxy);
      grzy->SetMarkerColor(in+1);
      grzy->SetMarkerStyle(20);
      grzy->SetMarkerSize(0.40);
      
      grxy->SetPoint(grzy->GetN(),x0.X(),x0.Y());
      grzy->SetPoint(grzy->GetN(),x0.Z(),x0.Y());
      fillgraph = true;
    }


    double Eterm = E;
    double w=1;					//viene aggiornato da count
    bool cnt= Count(x0,d0,Eterm,w,grxy,grzy,fillgraph);
    if (cnt){
      ndet++; //conteggio particelle rivelate
      hc.Fill(0.5,w);
    }
    
    if (fillgraph){
      c->cd(1);
      grxy->Draw("LP");
      c->cd(2);
      grzy->Draw("LP");
    }
    hE->Fill(Eterm,w);
  }

  ndet=hc.GetBinContent(1);
  double endet=hc.GetBinError(1);
  double p=ndet/nneu;
  double ep=endet/nneu;
  double var=endet*endet/(ndet*ndet);

  cout << ndet << endl;
  cout << "measured flux " << ndet/(detr*detr*TMath::Pi()) << endl;
  cout << "frazione rivelata " << p << " +/- " << ep << endl;
  double FOM=1/(var*t.CpuTime());
  cout << "FOM="<< FOM << endl;

}

#ifndef __CINT__
int main(){
  TApplication app("app",0,NULL);
  SimNeuVar();
  gPad->Update();
  app.Run(true);
}
#endif
