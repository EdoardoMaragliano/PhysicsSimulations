"""
This script simulates the decay of the J/psi particle in the center-of-mass (CM) reference frame 
and then transforms the decay products into the laboratory (lab) reference frame. It generates random 
angles in the CM frame, computes the decay products' four-momenta, fills a histogram in the CM frame 
with the angle between the electron and positron, boosts these decay products to the lab frame, fills 
another histogram with the angle between the boosted particles, and then plots both histograms.

The decay process simulated is:
    J/psi -> e+ + e-

The script utilizes the ROOT framework for histogram plotting and random number generation.

Note: Proper setup of the ROOT environment is required to execute this script.
"""


import numpy as np
import matplotlib.pyplot as plt
from math import *
from ROOT import TRandom3, TMath, TH1D, TLorentzVector, TVector3, TCanvas, gApplication

# Initialize random number generator
rnd = TRandom3()

# J/psi mass and electron mass (in GeV)
mj = 3.096916
me = 0.511e-3

# J/psi momentum (in GeV)
pj = 100.

# Computing pstar
pstar = 1 / (2 * mj) * sqrt((mj**2 - (me + me)**2) * (mj**2 - (me - me)**2))
Estar = sqrt(pstar**2 + me**2)

# Histograms for angles in CM and lab reference frames
hCM = TH1D("hCM", "", 100, 0, 3.2)
hLab = TH1D("hLab", "", 100, 0, 3.2)

# Create TLorentzVector for J/psi along z-axis
bp = TLorentzVector(0, 0, pj, sqrt(pj**2 + mj**2))
bv = bp.BoostVector()  # Boost vector for J/psi

for i in range(1, 1000):
    phi = 2 * TMath.Pi() * rnd.Rndm()
    theta = acos(2 * rnd.Rndm() - 1)  # Generate angles phi and theta in CM
    vp = TVector3()
    vp.SetMagThetaPhi(pstar, theta, phi)

    # Create TLorentzVector for e+ and e- in CM
    pv1 = TLorentzVector(vp, Estar)
    pv2 = TLorentzVector(-vp, Estar)
    hCM.Fill(pv1.Angle(pv2.Vect()))  # Fill CM histogram

    # Boost electrons to the lab reference frame
    pv1.Boost(bv)  # Update pv1 to lab frame
    pv2.Boost(bv)  # Update pv2 to lab frame
    ang = pv1.Angle(pv2.Vect())
    hLab.Fill(ang)  # Fill lab histogram


# Plot histograms
can = TCanvas()
can.Divide(2)
can.cd(1)
hCM.Draw()
hCM.GetXaxis().SetTitle("Angle [radians] in CM frame")
hCM.GetYaxis().SetTitle("Entries")
can.cd(2)
hLab.GetXaxis().SetTitle("Angle [radians] in Lab frame")
hLab.GetYaxis().SetTitle("Entries")
hLab.Draw()
gApplication.Run(True)

