{

  double mu1 = 1;
  double mu2 = 2;
  double s1  = 0.2;
  double s2  = 0.3;
  double rho = 0.5;

  TRandom3 rnd;

  TH2D h2("h2","",100,-3,5,100,-3,5);
  for (int i=0;i<1000;i++){
    double xi1 = rnd.Gaus(0,1);
    double xi2 = rnd.Gaus(0,1);
    double x1  = mu1 + s1*xi1;
    double x2  = mu2 + s2*rho*xi1 + sqrt(1-rho*rho)*s2*xi2;
    h2.Fill(x1,x2);
  }
  h2.Draw();
  cout << h2.GetCorrelationFactor() << endl;
  
}
