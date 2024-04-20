import numpy as np
import matplotlib.pyplot as plt
from math import *
from ROOT import *
from array import array
import time as t

#decadimento H3 -> He3+ + e- + v

uA  = 0.9315                ##massa atomo H1 (1amu=931.5MeV)
mv  = 1e-9
me  = 0.511e-3              ##massa elettrone
mH  = 3.0160492787*uA       ##massa trizio   
mHe = 3.016029321*uA        ##massa Elio

Q   = mH-mHe                ##energia disponibile per e- e v
mHe = mHe-me

P       = TLorentzVector(0.0,0.0,0.0,mH)        ##quadrimpulso iniziale H (fermo)
masses  = np.array([mHe,me,mv])                 ##array delle masse
massesn = np.array([mHe,me,mv*1000])             ##array delle masse per v troppo massivo(serve a vedere la deviazione nel plot di Kurie)

h      = TH1D("h1","",200,(Q-mv)*0.7,(Q-mv)*1.05)   ##distr energia e-
htotw  = TH1D("h3","",200,0,(Q-mv)*1.05)            ##distr energia tot
hn     = TH1D("h1n","",200,(Q-mv)*0.7,(Q-mv)*1.05)  
hntotw = TH1D("h3n","",200,0,(Q-mv)*1.05)


gRandom.SetSeed(int(t.time()))

event   = TGenPhaseSpace()                      ##genero spazio delle fasi
event.SetDecay(P, 3, masses)                    ##imposto spazio delle fasi a 3 corpi
eventn  = TGenPhaseSpace()
eventn.SetDecay(P, 3, massesn)

nev = 100000                                    ##numero neutrini
for i in range(0,nev):
    weight = event.Generate()
    el     = event.GetDecay(1)                      ##quadrivettore elettrone  
    nu     = event.GetDecay(2)                      ##quadrivettore neutrino
    Nu     = event.GetDecay(0)                      ##quadrivettore Elio?
    h.Fill(el.E()-me,weight*el.E()*nu.E()*Nu.E())
    htotw.Fill(el.E()-me,weight*el.E()*nu.E()*Nu.E())

h.SetMarkerStyle(20)
h.SetMarkerSize(0.5)
h.Draw("EP")
h.Scale(htotw.GetEntries()/htotw.Integral())

for i in range(0,nev):
    weight = eventn.Generate()
    el     = eventn.GetDecay(1)
    nu     = eventn.GetDecay(2)
    Nu     = eventn.GetDecay(0)
    hn.Fill(el.E()-me,weight*el.E()*nu.E()*Nu.E())
    hntotw.Fill(el.E()-me,weight*el.E()*nu.E()*Nu.E())

hn.Draw("HSAME")
hn.Scale(hntotw.GetEntries()/hntotw.Integral())
hn.SetFillColor(7)
hn.SetFillStyle(1001)
h.Draw("SAMEEP")
gApplication.Run(True)
