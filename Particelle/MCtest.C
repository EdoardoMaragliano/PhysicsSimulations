{
	double x1=1,x2=2;
	double s1=0.2,s2=0.3;
	double rho=0.5;

	TRandom3 rnd;
	rnd.SetSeed(12345678);

	TH2D* h = new TH2D("h2","",100,0,3,100,0,3);

	for (int i = 0; i < 10000; ++i)
	{
		double xi1= rnd.Gaus(0,1);
		double xi2= rnd.Gaus(0,1);
		double v1= x1+s1*xi1;
		double v2= x2+rho*s2*xi1 + sqrt(1-rho*rho)*s2*xi2; 
		h->Fill(v1,v2);
	}

	h->Draw();
	cout << h->GetCorrelationFactor() << endl;


}