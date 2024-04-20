int main(){
  TApplication app("app",0,NULL);
  TLorentzVector *in=new TLorentzVector(0,0,pbeam,sqrt(pbeam*pbeam+pip->Mass()*pip->Mass()));
  Particle pion(pip,in);
  

  TH1D *Emeas = new TH1D("E(last measurment)","E(last measurment)",100,50,100);
  Emeas     ->  GetXaxis()->SetTitle("dEdx [MeV/cm]");
  //GRAFICI MEDIE TRONCATE
  TH1D *h1=new TH1D("Ecut 90%","Ecut 90%",100,65,75);
  TH1D *h2=new TH1D("Ecut 80%","Ecut 80%",100,65,75);
  TH1D *h3=new TH1D("Ecut 70%","Ecut 70%",100,65,75);
  TH1D *h4=new TH1D("Ecut 60%","Ecut 60%",100,65,75);
  h1    ->  GetXaxis()->SetTitle("dEdx [MeV/cm]");
  h2    ->  GetXaxis()->SetTitle("dEdx [MeV/cm]");
  h3    ->  GetXaxis()->SetTitle("dEdx [MeV/cm]");
  h4    ->  GetXaxis()->SetTitle("dEdx [MeV/cm]");

  TH1D hCut[4]={*h1,*h2,*h3,*h4};

  TGraph *grRMS=new TGraph();
  grRMS->SetTitle("RMS");
  TGraph *grDiff=new TGraph();
  grDiff->SetTitle("E_reco-E_gen");
    
  //NB in PDG le masse sono in GeV
  //NB in PDG le cariche sono in abs(e)/3
 
  double dEdx;
  double cut[4]={0.9,0.8,0.7,0.6};            //percentuali troncamento
  vector <double> data;
  double BB=BetheBloch(pion,Z,A);
  cout <<"BB="<<BB<<endl;

  for(int n=0; n<nEv;++n){
    data.clear();
    
    for(int i=0;i<nSamp;++i){                     //genero 20 campionamenti e ottengo la Landau
      dEdx=rnd.Landau(BB,1);
      data.push_back(dEdx); 
    }
    sort(data.begin(),data.end());               //ordino i dati

    for(int j=0;j<4;j++){                       //per ciascun troncamento
      for(int i=0;i<cut[j]*data.size();++i){    
        Emeas->Fill(data[i]);                       //fai il grafico troncato
      }
      double mean=Emeas->GetMean();          //estraggo la media dal gr troncato
      Emeas->Reset();                        //resetto il grafico di Emeas                   
      hCut[j].Fill(mean);                    //riempio hCut con la media
    }
  }
    for(int i=0;i<4;++i)
      hCut[i].Fit("gaus","L");
  for(auto mis:data)
    Emeas->Fill(mis);     //riempio con l'ultima misura
  for(int j=0;j<4;++j){
    grRMS->SetPoint(grRMS->GetN(),cut[j],hCut[j].GetStdDev());    //grafico rms vs troncamento
    grDiff->SetPoint(grDiff->GetN(),cut[j],abs(BB-hCut[j].GetMean()));
  }


  TCanvas *can1=new TCanvas("Emeas","Emeas",500,500);
  Emeas->Draw();
  TCanvas *can2=new TCanvas("RMS","RMS",800,600);
  can2->Divide(2);
  can2->cd(1);
  grDiff->Draw("A*");
  can2->cd(2);
  grRMS->Draw("A*");
  TCanvas *can3=new TCanvas("Medie troncate","Medie Troncate",800,600);
  can3->Divide(4);
  for(int i=0;i<4;i++){
    can3->cd(i+1);
    hCut[i].Draw();
  }

  cout << "quit ROOT to continue"<<endl;
  app.Run(true);

  cout <<"uso cut off al 60%" << endl;
  //distuggo un po' di oggetti per liberare memoria
  h3->~TH1D(); h2->~TH1D(); h1->~TH1D(); grDiff->~TGraph(); grRMS->~TGraph(); Emeas->~TH1D();
  can1->~TCanvas(); can2->~TCanvas(); can3->~TCanvas();

  //tengo solo una media troncata e la uso come E ricostruita
  //trascuro risoluzione sull'energia

  double TOF=L/c*sqrt(1+pow(pion.Mass()*c/pion.P(),2));
  cout << "TOF="<<TOF<<endl;
  double tres=100e-12; //sec
  rnd.SetSeed(1627479);
  TOF=rnd.Gaus(TOF,tres);
  cout << "TOF measurment="<<TOF<<endl;

  
  app.Run(true);






  return 0;
}