/*
    //decad D*-, non salvo i prodotti perché non lo seguo più
    if(rnd.Rndm()<pD0){
      //cout << "D*m->D0+pim"<<endl;
      double masses[2]={D0->Mass(),pim->Mass()};  //D*m->D0+pim
      tgps.SetDecay(*v_Dm,2,masses);               
      w=w*tgps.Generate();                        //cout << "generate" << endl;
      Particle pion1(pim,tgps.GetDecay(1));
     // products.push_back(pion1);                   //aggiungo pione 1 ai prodotti
     // cout << "D0-> km+pip" << endl;
      TLorentzVector* vD0=tgps.GetDecay(0);       
      double m[2]={km->Mass(),pip->Mass()};       //decadimento D0 in k e pi
      tgps.SetDecay(*vD0,2,m);
      w=w*tgps.Generate();
    
      Particle pion2(pip,tgps.GetDecay(1));
      Particle kaon(km,tgps.GetDecay(0));
      //products.push_back(pion2);                  //aggiungo pione 2 ai prodotti
      //products.push_back(kaon);                   //aggiungo kaone ai prodotti

      //cout << pion1.P() << setw(10) << pion2.P() << setw(10) << kaon.P() << endl;

    } else {  

      //decadimento Dstar in 3 pi
      //cout << "D*m->pim+pip+pim" << endl;
      double masses[3]={pim->Mass(),pip->Mass(),pim->Mass()};
      tgps.SetDecay(*v_Dm,3,masses);
      w=w*tgps.Generate();
      Particle pion1(pim,tgps.GetDecay(0));
      Particle pion2(pip,tgps.GetDecay(1));
      Particle pion3(pim,tgps.GetDecay(2));
      //products.push_back(pion1);
      //products.push_back(pion2);
      //products.push_back(pion3);

      //cout << pion1.P() << setw(10) << pion2.P() << setw(10) << pion3.P() << endl;
    }*/
