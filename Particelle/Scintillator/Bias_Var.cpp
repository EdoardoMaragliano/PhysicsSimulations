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

// g++ bias_var_var.cpp -o bias_var_var.x `root-config --cflags --glibs`

using namespace std;

double me   = 0.511;   
double n    = 3.67*1e6/(149.89)*6.022e23;

double r    = 0.02;
double dz   = 0.025;
double dist = 0.10;
double   E[3]={1.275,0.511,0.511};        //energie dei fotoni
int    nev  = 1e6;
int    ngr  = 100;

TStopwatch tstop;

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
    double PE_ov_tot = sigma_PE/sigma_tot;
    if(rnd.Rndm()<PE_ov_tot){ 					//c'era Comp ma secondo me ci va tot
      // Effetto fotoelettrico
    	Ed+=E;	
    	E=0;						//semplicemente assorbito -> esce dal do-while

    } else {
      // Effetto Compton
    	double cosTheta=1-(pow(1+2*E/me,rnd.Rndm())-1)*me/E;
    	double sinTheta=sin(TMath::ACos(cosTheta));
    	double phi=2*TMath::Pi()*rnd.Rndm();
   		//calcolo d
    	if(d0.Z()==1 || d0.Z()==-1){
    		d.SetX(sinTheta*cos(phi));
    		d.SetY(sinTheta*sin(phi));
    		d.SetZ(cosTheta); 
    	}
    	else{
    		d.SetX(d0.X()*cosTheta+sinTheta/sqrt(1-pow(d0.Z(),2))*(d0.X()*d0.Z()*cos(phi)-d0.Y()*sin(phi)));
    		d.SetY(d0.Y()*cosTheta+sinTheta/sqrt(1-pow(d0.Z(),2))*(d0.Y()*d0.Z()*cos(phi)-d0.X()*sin(phi)));
    		d.SetZ(d0.Z()*cosTheta-sqrt(1-pow(d0.Z(),2))*sinTheta*cos(phi));
    	}
      double Eout= E/(1+E/me*(1-cosTheta));	
      d0=d;
      Ed+=E-Eout;		//energia assoribita dall'elettrone deflesso
      E=Eout;
    }
  } while (E!=0);

  return Ed;
}


void bias_var(){

  rnd.SetSeed(time(NULL));

  // Disegno geometria e istogramma
  DrawGeom(r,dz);
  int naccept=0;

  TH1D *hE = new TH1D("hEdep","",100,0,2);
  c->cd(3);
  hE->Draw();
  double coslim = dist/sqrt(dist*dist+r*r);
  double theta_max = atan(r/dist);
  
  // Itero sui fotoni
  for (int iev=0;iev<nev;iev++){
    
    double   th[3];   //theta 
    double   ph[3];   //phi
    TVector3 u[3];    //direzioni
    TVector3 x0,d;
    int bias_varmask[3]={0,0,0};              //per riduzione varianza
    double w=1;                               //peso iniziale

    double pgeo=(1-coslim)/2;         //probabilità geometrica (sono nel giusto angolo solido?)
    double p1=pgeo*(1-pgeo);          //fotone 1 nel rilevatore
    double p2=2*pgeo*(1-pgeo);        //fotone 2 nel rilevatore
    double p12=2*pgeo*pgeo;           //entrambi nel rilevatore
    w=p1+p2+p12;
    double test=rnd.Rndm()*w;
    if (test<p1){
      bias_varmask[0]=1;                  //solo il primo
    } else if(test<p1+p2){
      bias_varmask[1]=1;                  //solo il secondo
    } else {
      bias_varmask[0]=1;
      bias_varmask[1]=1;                  //entrambi
    }
    
    double Ed=0;
    
    for (int k=0;k<3;k++){
      if(bias_varmask[k]==0) continue;  //se bias_varmask è 0 salto all'iter succ
      naccept++;
      double th=acos((1-coslim)*rnd.Rndm()+coslim);
      double ph=2*M_PI*rnd.Rndm();
      d=TVector3(sin(th)*cos(ph),sin(th)*sin(ph),cos(th));
      double R=dist/d.Z();
      x0=TVector3(R*d.X(),R*d.Y(),0);

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

      Ed += Edep(x0,d,E[k],grxy,grzy,fillgraph);
    
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

void isotropic(){

  rnd.SetSeed(time(NULL));

  // Disegno geometria e istogramma
  DrawGeom(r,dz);
  TH1D *hE = new TH1D("Edep","",100,0,2);
  c->cd(3);
  hE->Draw();
  int n = 0;
  double thMax=atan(r/dist);

	// Itero sui fotoni  
	for (int iph=0;iph<nev;iph++){
	 double Ed = 0;
	 double theta[3],phi[3];
	 TVector3 v[3], x0,d0;

	 theta[0] = acos(2*rnd.Rndm()-1);
	 theta[1] = acos(2*rnd.Rndm()-1);
	    
	 phi[0] = 2*M_PI*rnd.Rndm();
	 phi[1] = 2*M_PI*rnd.Rndm();
	    
		for (int i=0;i<2;i++){
			v[i].SetX(sin(theta[i])*cos(phi[i]));
			v[i].SetY(sin(theta[i])*sin(phi[i]));
			v[i].SetZ(cos(theta[i]));
		}
		v[2] = -v[1];
		    
		    
		//controllo che abbia preso il rilevatore
		for(int i = 0; i<3;i++){
	        if(acos(v[i].Z())<thMax){ 
			      n++;
		        x0 =TVector3(dist*v[i].X()/v[i].Z(),dist*v[i].Y()/v[i].Z(),0);
			      d0 = v[i];
		      
            TGraph *grxy,*grzy;
            bool fillgraph=false;
            if (n<=ngr){
              grxy = new TGraph();
              grxy->SetMarkerColor(n+1);
              grxy->SetMarkerStyle(20);
              grxy->SetMarkerSize(0.40);
              grzy = new TGraph(*grxy);
              grxy->SetPoint(grzy->GetN(),x0.X(),x0.Y());
              grzy->SetPoint(grzy->GetN(),x0.Z(),x0.Y());
              fillgraph = true;
			    }

			    Ed += Edep(x0,d0,E[i],grxy,grzy,fillgraph);
			    // cout << Ed << endl;
			    // Partendo dall'energia depositata estraggo l'energia misurata.
			    // si assuma un errore relativo sull'energia di
			    // (sigma(E)/E)^2 = somma in quad(0.02 + 0.04/sqrt(E))
		    
		      
			    if (fillgraph){
			      c->cd(1);
			      grxy->Draw("LP");
			      c->cd(2);
			      grzy->Draw("LP");
			    }
	        } 
	        if (Ed!=0){
			    Ed = rnd.Gaus(Ed, sqrt(0.02*0.02+0.04*0.04/Ed)*Ed);
			    hE->Fill(Ed);  
			}

		}
	}
  double sig=hE->GetBinContent(50);   //prendo un bin stabile (una spalla del grafico)
  double err=hE->GetBinError(50);     //prendo il suo errore
  double FOM=1/(tstop.CpuTime()*pow(err/sig,2));
  cout <<"FOM= "<<FOM<<endl;
}

void lin_beam(){

  rnd.SetSeed(time(NULL));

  // Disegno geometria e istogramma
  DrawGeom(r,dz);
  TH1D *hE = new TH1D("Edep","",100,0,2);
  c->cd(3);
  hE->Draw();
  
  // Itero sui fotoni  
  for (int iph=0;iph<nev;iph++){
    
    TVector3 d0(0,0,1);
    TVector3 x0(0,0,0);
    
    TGraph *grxy,*grzy;
    bool fillgraph=false;
    if (iph<=ngr){
      grxy = new TGraph();
      grxy->SetMarkerColor(iph+1);
      grxy->SetMarkerStyle(20);
      grxy->SetMarkerSize(0.40);
      grzy = new TGraph(*grxy);
      grxy->SetPoint(grzy->GetN(),x0.X(),x0.Y());
      grzy->SetPoint(grzy->GetN(),x0.Z(),x0.Y());
      fillgraph = true;
    }

    double Ed = Edep(x0,d0,E[0],grxy,grzy,fillgraph);
    
    // Partendo dall'energia depositata estraggo l'energia misurata.
    // si assuma un errore relativo sull'energia di
    // (sigma(E)/E)^2 = somma in quad(0.02 + 0.04/sqrt(E))
    if (Ed!=0){
      double sigmaE=Ed*sqrt(0.02*0.02+pow(0.04/sqrt(Ed),2));
      hE -> Fill(rnd.Gaus(Ed,sigmaE));
    }
    
    if (fillgraph){
      c->cd(1);
      grxy->Draw("LP");
      c->cd(2);
      grzy->Draw("LP");
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
    cout << "lin_beam"<< endl << "isotropic" << endl <<"bias_var" << endl;
    return 1;
  }
  TApplication app("app",0,NULL);
  if(string(argv[1])=="bias_var") bias_var();
  if(string(argv[1])=="isotropic") isotropic();
  if(string(argv[1])=="lin_beam") lin_beam();
  gPad->Update();
  app.Run(true);
}
#endif
