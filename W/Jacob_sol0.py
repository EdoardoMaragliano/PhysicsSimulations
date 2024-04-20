import numpy as np
import matplotlib.pyplot as plt
from math import *
from ROOT import *

rnd = TRandom3();

# unita' in GeV
mw   = 80.4
pw   = 10


pstar = mw/2

file = TFile("WJacobianPeak.root")
hpTe = file.Get("pTe")
# TH1D *hpTe = (TH1D*)file.Get("pTe");

hpTe1 = TH1D(hpTe) # con pW = 0 e risoluzione perfetta
hpTe2 = TH1D(hpTe) # con pW = 0 e risol. sper.
hpTe3 = TH1D(hpTe) # con pW != 0 e risol. sper.

for i in range (1,10000):
    the  = acos(2*rnd.Rndm()-1)
    phi  = 2*TMath.Pi()*rnd.Rndm()

    vp   = TVector3(pstar*cos(phi)*sin(the),pstar*sin(phi)*sin(the),pstar*cos(the))
    pe   = TLorentzVector(vp ,sqrt(pow(vp.Mag(),2)))
    pmu  = TLorentzVector(-vp,sqrt(pow(vp.Mag(),2)))

    pt = pe.Pt()        ##.Pt() componente trasversa
    if pt>16:
        hpTe1.Fill(pt)
        
    pt = rnd.Gaus(pe.Pt(),0.15*sqrt(pe.Pt()))       ##con errore sperimentale
    if pt>16:
        hpTe2.Fill(pt)

    ##impulso non trasverso non nullo (W+e-)  -> non è Lorentz Inv 

    pWT  = -pw*log(rnd.Rndm())      ##p trasv W expo con media pw=1,5,10GeV
    phiW = 2*TMath.Pi()*rnd.Rndm()
    pWx  = pWT*cos(phiW)
    pWy  = pWT*sin(phiW)
    vW   = TLorentzVector(pWx,pWy,0,sqrt(pWx*pWx+pWy*pWy+mw*mw))
    ##posso mettere qualunque z che tanto non si trasforma

    pe.Boost(vW.BoostVector())
    pt = rnd.Gaus(pe.Pt(),0.15*sqrt(pe.Pt()))
    if pt>16:
        hpTe3.Fill(pt)

    
        
hpTe.SetMarkerStyle(20)
hpTe.SetMarkerSize(1.2)
hpTe.Draw("PE")
hpTe.SetMaximum(hpTe.GetMaximum()*3)
hpTe1.Draw("SAME")
hpTe1.SetLineColor(2)
hpTe1.Scale(hpTe.Integral()/hpTe1.Integral())
hpTe2.Draw("SAME")
hpTe2.SetLineColor(4)
hpTe2.Scale(hpTe.Integral()/hpTe2.Integral())
hpTe3.Draw("SAME")
hpTe3.SetLineColor(6)
hpTe3.Scale(hpTe.Integral()/hpTe3.Integral())
gApplication.Run(True)
