"""
Script Name: H3_Decay_Simulation.py
Description: This script simulates the decay process of tritium (H3) into helium-3, an electron, and a neutrino. 
It calculates the energy distribution of the emitted electron in both regular and high-mass neutrino scenarios and plots them for comparison.

H3 -> He3+ + e- + v

Dependencies:
    - numpy
    - ROOT (PyROOT)
    - array
    - time

Usage:
    - Make sure to have ROOT installed and properly configured.
    - Run the script in a Python environment with the necessary dependencies installed.

Note: Before running the script, ensure that you have the necessary ROOT environment properly set up and configured.
"""

import numpy as np
from math import *
from ROOT import TLorentzVector, gApplication, TH1D, TGenPhaseSpace, gRandom
from array import array
import time as t


uA  = 0.9315                ## atomic mass of H1 (1amu=931.5 MeV)
mv  = 1e-9
me  = 0.511e-3              ## electron mass
mH  = 3.0160492787*uA       ## tritium mass   
mHe = 3.016029321*uA        ## helium mass

Q   = mH-mHe                # Energy available for electron and neutrino
mHe = mHe-me

P       = TLorentzVector(0.0,0.0,0.0,mH)        ## inital 4-momentum of H, not moving
masses  = np.array([mHe,me,mv])                 ## array of masses
massesn = np.array([mHe,me,mv*1000])            ## array of masses assuming a very massive neutrino -> see Kurie Plot endpoint

h      = TH1D("h1","",200,(Q-mv)*0.7,(Q-mv)*1.05)   ## energy distribution e-
htotw  = TH1D("h3","",200,0,(Q-mv)*1.05)            ## total energy distribution
hn     = TH1D("h1n","",200,(Q-mv)*0.7,(Q-mv)*1.05)  
hntotw = TH1D("h3n","",200,0,(Q-mv)*1.05)


gRandom.SetSeed(int(t.time()))

event   = TGenPhaseSpace()                      ## generate phase space
event.SetDecay(P, 3, masses)                    ## setting 3 body phase space
eventn  = TGenPhaseSpace()
eventn.SetDecay(P, 3, massesn)

nev = 100000                                    ## number of neutrinos
for i in range(0,nev):
    weight = event.Generate()
    el     = event.GetDecay(1)                      ## electron 4-momentum   
    nu     = event.GetDecay(2)                      ## neutrino 4-momentum
    Nu     = event.GetDecay(0)                      
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
