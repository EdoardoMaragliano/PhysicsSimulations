{

    TGraphErrors gr;
	TCanvas can;
    gr.SetTitle("leggi di scala");
	gr.GetXaxis()->SetTitle("D/F");
	gr.GetYaxis()->SetTitle("nIsole");
	gr.SetMarkerStyle(7);
	TF1 f("f","[0]/(x**(1./3))",1E03,10E09);
	f.SetParameter(0,1000);
    can.SetLogx();
	
	ifstream filein("scala.dat");
	double x,y;
    vector<double> DF,N;
	while(filein>>x>>y){
        DF.push_back(x);
        N.push_back(y);
    }
    
    for(int i=0;i<DF.size();++i){
        double er=0.08*y+1;
    
		gr.SetPoint(i,DF[i],N[i]);
        gr.SetPointError(i,0.01,1);
	}
	gr.Draw("AP");
	gr.Fit("f","R");
}
