import numpy as np
import matplotlib.pyplot as plt
from math  import *
from ROOT  import *

rnd = TRandom3()
cx = TCanvas("cx","",5,5,400,400)

# unita' in GeV
mw_nom   = 80.4
gw       =  2       ##FWHM del W in GeV
pw       =  10
me       =  0.511e-3
mnu      =  0

pstar = 1/(2*mw_nom)*sqrt((mw_nom**2-(me+mnu)**2)*(mw_nom**2-(me-mnu)**2))

file = TFile("WJacobianPeak.root")
hpTe = file.Get("pTe")
hpTW = file.Get("pTW")
print("pTW ",hpTW.GetMean())

hpTeMC0 = TH1D(hpTe) # 
hpTeMC1 = TH1D(hpTe) # larghezza W
hpTeMC2 = TH1D(hpTe) # risol. sperimentale
hpTeMC3 = TH1D(hpTe) # con pTW != 0  
hpTeMC4 = TH1D(hpTe) # e nu nu

hpTeMC0.Reset()
hpTeMC1.Reset()
hpTeMC2.Reset()
hpTeMC3.Reset()
hpTeMC4.Reset()

tgs    = TGenPhaseSpace()
masses = np.array([me,0,0])
for i in range (1,10000):
    the  = acos(2*rnd.Rndm()-1)
    phi  = 2*TMath.Pi()*rnd.Rndm()
    pstar = 1/(2*mw_nom)*sqrt((mw_nom**2-(me+mnu)**2)*(mw_nom**2-(me-mnu)**2))

    # Step 0
    vp   = TVector3()
    vp.SetMagThetaPhi(pstar,the,phi)
    pe = TLorentzVector(vp ,sqrt(pstar*pstar+me*me))
    pt = pe.Pt()
    if pt>16:
        hpTeMC0.Fill(pt)

    # Step 1: genero le masse di W secondo la BreitWigner relativistica
    mw2    = mw_nom*mw_nom+ gw*mw_nom*tan(TMath.Pi()*(rnd.Rndm()-0.5))
    if mw2<0: ##escludo masse <0
        continue
    mw     = sqrt(mw2)
    pstar = 1/(2*mw)*sqrt((mw**2-(me+mnu)**2)*(mw**2-(me-mnu)**2))
    ##aggiorno pstar poiche' ho cambiato la massa
    
    vp.SetMagThetaPhi(pstar,the,phi)
    pe = TLorentzVector(vp ,sqrt(pstar*pstar+me*me))
    pt = pe.Pt()
    if pt>16:
        hpTeMC1.Fill(pt)

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

    #Step : genero decadimento a tre corpi
    tgs.SetDecay(vW,3,masses)               ##(P,nPart,masse)
    weight = tgs.Generate()                 ##genero il decadimento
    pe     = tgs.GetDecay(0)                ##quadrivetore elettrone
    pt     = pe.Pt()                        
    pt     = rnd.Gaus(pt,0.15*sqrt(pt))     ##aggiungo risoluzione in energia del rivelatore
    if pt>16:
        hpTeMC4.Fill(pt)

hpTe.SetMarkerStyle(20)
hpTe.SetMarkerSize(0.8)
hpTeMC0.SetLineWidth(2)
hpTeMC1.SetLineWidth(2)
hpTeMC2.SetLineWidth(2)
hpTeMC3.SetLineWidth(2)
hpTeMC4.SetLineWidth(2)

hpTeMC0.SetLineColor(2)     ##rosso
hpTeMC1.SetLineColor(3)     ##verde
hpTeMC2.SetLineColor(4)     ##blu
hpTeMC3.SetLineColor(6)     ##magenta
hpTeMC4.SetLineColor(8)

hpTe.Draw("PE")
hpTe.SetMaximum(hpTe.GetMaximum()*2.5)
hpTe.GetXaxis().SetTitle("p^{T}_{e}")
hpTe.GetYaxis().SetTitle("N. events")

hpTeMC0.Draw("HSAME")
hpTeMC0.Scale(hpTe.Integral()/hpTeMC1.Integral())
hpTeMC1.Draw("HSAME")
hpTeMC1.Scale(hpTe.Integral()/hpTeMC1.Integral())
hpTeMC2.Draw("HSAME")
hpTeMC2.Scale(hpTe.Integral()/hpTeMC2.Integral())
hpTeMC3.Draw("HSAME")
hpTeMC3.Scale(hpTe.Integral()/hpTeMC3.Integral())
hpTeMC4.Draw("HSAME")
hpTeMC4.Scale(hpTe.Integral()/hpTeMC4.Integral())
print(hpTe.Chi2Test(hpTeMC0,"UW"))
print(hpTe.Chi2Test(hpTeMC1,"UW"))
print(hpTe.Chi2Test(hpTeMC2,"UW"))
print(hpTe.Chi2Test(hpTeMC3,"UW"))
print(hpTe.Chi2Test(hpTeMC4,"UW"))
gApplication.Run(True)
