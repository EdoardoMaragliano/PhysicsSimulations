import numpy as np
import matplotlib.pyplot as plt
from math  import *
from ROOT  import *

rnd = TRandom3()
cx = TCanvas("cx","",5,5,400,400)

# unita' in GeV
mw_nom   = 80.4
gw       =  2
pw       =  10
me       =  0.511e-3
mnu      =  0

pstar = 1/(2*mw_nom)*sqrt((mw_nom**2-(me+mnu)**2)*(mw_nom**2-(me-mnu)**2))

file = TFile("WJacobianPeak.root")
hpTe = file.Get("pTe")
hpTW = file.Get("pTW")
print("pTW ",hpTW.GetMean())

hpTeMC0 = TH1D(hpTe) # 
#hpTeMC1 = TH1D(hpTe) # larghezza W
hpTeMC2 = TH1D(hpTe) # risol. sperimentale
hpTeMC3 = TH1D(hpTe) # con pTW != 0  
#hpTeMC4 = TH1D(hpTe) # e nu nu
hpTeMC0.Reset()
hpTeMC2.Reset()
hpTeMC3.Reset()

tgs    = TGenPhaseSpace()
masses = np.array([me,0,0]) ##elettrone e due neutrini
for i in range (1,10000):
    
    mw = mw_nom
    
    the  = acos(2*rnd.Rndm()-1)
    phi  = 2*TMath.Pi()*rnd.Rndm()
    pstar = 1/(2*mw)*sqrt((mw**2-(me+mnu)**2)*(mw**2-(me-mnu)**2))

    # Step 0
    vp   = TVector3()
    vp.SetMagThetaPhi(pstar,the,phi)
    pe = TLorentzVector(vp ,sqrt(pstar*pstar+me*me))
    pt = pe.Pt()
    if pt>16:
        hpTeMC0.Fill(pt)

    #Step 2
    pt = rnd.Gaus(pt,0.15*sqrt(pt))
    if pt>16:
        hpTeMC2.Fill(pt)

    #Step 3
    pWT  = -pw*log(rnd.Rndm())
    phiW = 2*TMath.Pi()*rnd.Rndm()
    pWx  = pWT*cos(phiW)
    pWy  = pWT*sin(phiW)
    vW   = TLorentzVector(pWx,pWy,0,sqrt(pWx*pWx+pWy*pWy+mw*mw))

    pe.Boost(vW.BoostVector())
    pt = pe.Pt()
    pt = rnd.Gaus(pt,0.15*sqrt(pt))
    if pt>16:
        hpTeMC3.Fill(pt)


hpTe.SetMarkerStyle(20)
hpTe.SetMarkerSize(0.8)
hpTeMC0.SetLineWidth(3)
hpTeMC2.SetLineWidth(3)
hpTeMC3.SetLineWidth(3)

hpTeMC0.SetLineColor(2)
hpTeMC2.SetLineColor(4)
hpTeMC3.SetLineColor(6)

hpTe.Draw("PE")
hpTe.SetMaximum(hpTe.GetMaximum()*2.5)
hpTe.GetXaxis().SetTitle("p^{T}_{e}")
hpTe.GetYaxis().SetTitle("N. events")

hpTeMC0.Draw("HSAME")
hpTeMC0.Scale(hpTe.Integral()/hpTeMC0.Integral())
hpTeMC2.Draw("HSAME")
hpTeMC2.Scale(hpTe.Integral()/hpTeMC2.Integral())
hpTeMC3.Draw("HSAME")
hpTeMC3.Scale(hpTe.Integral()/hpTeMC3.Integral())

print(hpTe.Chi2Test(hpTeMC0,"UW")) ##unwighted ##weighted
print(hpTe.Chi2Test(hpTeMC2,"UW")) ##dati      ##MC  
print(hpTe.Chi2Test(hpTeMC3,"UW")) ##MC è pesato perché ho normalizzato con Scale()

gApplication.Run(True)
