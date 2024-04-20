{
  
  double tau0 = 2; 
  double tau1 = 5;
  // eta > f/g > lambda/lambda_0 , f/g > tau_0/tau
  TRandom3 rnd;
  rnd.SetSeed(12345678);

  TH1D h("h","",20,0,40);
  for (int i=0;i<10000000;i++){
    double told=0,t=0;
    do {
      t    = tau0*(told/tau0-log(rnd.Rndm()));
      told = t;
    } while ( rnd.Rndm() > tau0/tau1 );
    h.Fill(t);
  }
  h.Draw();

}