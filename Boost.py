##DECADIMENTO JPSI IN SRMC -> SRLAB


import numpy as np
import matplotlib.pyplot as plt
from math import *
from ROOT import *

rnd = TRandom3();

mj   = 3.096916
me   = 0.511e-3
pj   = 100.

#Calcolo pstar
pstar=1/(2*mj)*sqrt((mj**2-(me+me)**2)*(mj**2-(me-me)**2))
Estar=sqrt(pstar**2+me**2)

hCM  = TH1D("hCM","",100,0,3.2)
hLab = TH1D("hLab","",100,0,3.2)

bp=TLorentzVector(0,0,pj,sqrt(pj**2+mj**2))   # Creo TLorentzVector per J/psi lungo z
bv=bp.BoostVector()

for i in range (1,1000):
    phi = 2*TMath.Pi()*rnd.Rndm()
    # Genero angoli phi e theta nel CM
    theta=acos(2*rnd.Rndm()-1)
    vp=TVector3 ()
    vp.SetMagThetaPhi(pstar,theta,phi)

# Creo TLorentzVector per e+ e-
    pv1=TLorentzVector(vp,Estar)
    pv2=TLorentzVector(-vp,Estar)
    # Riempo hCM
    hCM.Fill(pv1.Angle(pv2.Vect()))
  

    
    # Boosto el. e pos. nel sistema del laboratorio (metodo Boost a cui bisogna passare, 
    #come vettore di boost, il valore ritornato da BoostVector applicato alla J/psi)
    pv1CM=TLorentzVector(pv1)
    pv2CM=TLorentzVector(pv2)
    pv1.Boost(bv)	#modifico pv1
    pv2.Boost(bv)	#modifico pv2
    ang=pv1.Angle(pv2.Vect())
    hLab.Fill(ang)
    
    # Riempo hLab
can=TCanvas()
can.Divide(2)
can.cd(1)
hCM.Draw()
can.cd(2)
hLab.Draw()
gApplication.Run(True);
