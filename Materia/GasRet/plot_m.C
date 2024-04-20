{
	TGraph gr;
	string str;
	cout << "nome file in" << endl;
	cin >> str;
	ifstream filein(str);
	int x, y, pieno;
	gr.SetMarkerStyle(7);
	while(filein>>x>>y>>pieno){
		if(pieno==1){
			gr.SetPoint(gr.GetN(),x,y);
		}
	}
	gr.Draw("AP");

}