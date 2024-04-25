/*
File Name: SimNeu.cpp
Description: This program simulates the interaction of neutrons with a 
scintillator detector and calculates the deposited energy. 
It also visualizes the geometry and neutron trajectories.
Author: [Your Name]

Dependencies:
    - ROOT (CERN software framework for data analysis and visualization)
Compilation:
    - g++ SimNeu.cpp -o SimNeu.x `root-config --cflags --glibs`

Usage:
    - Compile the program using the provided compilation command.
    - Execute the compiled binary (SimNeu.x) in a terminal or command prompt.

Note: Make sure ROOT is properly installed and configured before compiling and running the program.
*/

#include <TH1D.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TArc.h>
#include <TApplication.h>
#include <iostream>
#include <TGraph.h>
#include <TF1.h>
#include <TStopwatch.h>
#include <TFile.h>

using namespace std;

// Global variables
double A = 10; // Atomic number
double n = 2460 / A * 6.022e23; // Density

double dz = 1.00; // Material thickness
double detr = 0.05; // Radius of the neutron detector sphere
TVector3 det(0.0, 0.0, 0.5); // Position vector of the detector center

double E = 0.1; // Neutron energy in MeV
int nneu = 1e5; // Number of neutrons
int ngr = 100; // Number of plotted points for visualization

TRandom3 rnd;
TCanvas *c;

double sige = 18e-28; // Cross-section for neutron interaction

// Function to draw the geometry
void DrawGeom() {
    // Draw the volume of the scintillator and the neutron detector
    c = new TCanvas("c", "", 10, 10, 900, 320);
    c->Divide(3, 1);
    // View X-Y
    c->cd(1);
    gPad->Range(-dz, -dz, dz, dz);
    TGraph *grp1 = new TGraph;
    grp1->SetPoint(0, -dz, -dz);
    grp1->SetPoint(1, dz, -dz);
    grp1->SetPoint(2, dz, dz);
    grp1->SetPoint(3, -dz, dz);
    grp1->SetFillColor(11);
    grp1->Draw("F");
    TArc *arc1 = new TArc;
    arc1->SetFillColor(0);
    arc1->DrawArc(det[0], det[1], detr);
    // View Z-Y
    c->cd(2);
    gPad->Range(-dz, -dz, dz, dz);
    TGraph *grp2 = new TGraph;
    grp2->SetPoint(0, 0, -dz);
    grp2->SetPoint(1, dz, -dz);
    grp2->SetPoint(2, dz, dz);
    grp2->SetPoint(3, 0, dz);
    grp2->SetFillColor(11);
    grp2->Draw("F");
    TArc *arc2 = new TArc;
    arc2->SetFillColor(0);
    arc2->DrawArc(det[2], det[1], detr);
}

/**
 * @brief Simulates neutron interaction with a scintillator detector and calculates deposited energy.
 *
 * This function simulates the trajectory of a neutron within a scintillator detector and calculates the deposited energy. 
 * It iterates through the trajectory until the neutron's energy is completely deposited or it exits the detector volume.
 * If requested, it fills trajectory graphs for visualization purposes.
 *
 * @param x0 Initial position vector of the neutron.
 * @param d0 Initial direction vector of the neutron.
 * @param E Reference to the neutron's energy (in MeV), which will be updated during the simulation.
 * @param grxy Pointer to a TGraph object representing the neutron trajectory in the X-Y plane.
 * @param grzy Pointer to a TGraph object representing the neutron trajectory in the Z-Y plane.
 * @param fillgraph Boolean flag indicating whether to fill the trajectory graphs (grxy and grzy) or not.
 * @return True if the neutron interacts with the detector during the simulation, indicating a successful detection event.
 *         False if the neutron does not interact with the detector.
 */
bool Count(TVector3 x0, TVector3 d0, double &E, TGraph *&grxy, TGraph *&grzy, bool fillgraph) {
    bool cnt = false;
    do {
        // Calculate intersection with detector
        double s_intersec = -d0 * (x0 - det);
        TVector3 intersec = x0 + s_intersec * d0;
        if ((intersec - det).Mag() < detr) {
            cnt = true; // Neutron interacts with detector
        }
        
        // Calculate next interaction point
        TVector3 x, d;
        double lambda = 1 / (n * sige);
        double s = -lambda * log(rnd.Rndm());
        x = x0 + s * d0;
        x0 = x;
        
        // Check if neutron exits detector volume
        if (x.Z() < 0 || x.Z() > dz) {
            break;
        }
        
        // Fill trajectory graphs if requested
        if (fillgraph) {
            grxy->SetPoint(grxy->GetN(), x.X(), x.Y());
            grzy->SetPoint(grzy->GetN(), x.Z(), x.Y());
        }
        
        // Calculate scattering angles
        double thetaCM = acos(2 * rnd.Rndm() - 1);
        double theta = acos((1 + A * cos(thetaCM)) / sqrt(A * A + 1 + 2 * A * cos(thetaCM)));
        double phi = 2 * TMath::Pi() * rnd.Rndm();
        
        // Update neutron's energy
        E = E / pow((A + 1), 2) * pow(cos(theta) + sqrt(A * A - sin(theta) * sin(theta)), 2);
        
        // Update neutron's direction
        if (d0.Z() == 1) {
            d.SetX(sin(theta) * cos(phi));
            d.SetY(sin(theta) * sin(phi));
            d.SetZ(cos(theta));
        } else if (d0.Z() == -1) {
            d.SetX(-sin(theta) * cos(phi));
            d.SetY(-sin(theta) * sin(phi));
            d.SetZ(-cos(theta));
        } else {
            d.SetX(d0.X() * cos(theta) + sin(theta) / sqrt(1 - pow(d0.Z(), 2)) * (d0.X() * d0.Z() * cos(phi) - d0.Y() * sin(phi)));
            d.SetY(d0.Y() * cos(theta) + sin(theta) / sqrt(1 - pow(d0.Z(), 2)) * (d0.Y() * d0.Z() * cos(phi) + d0.X() * sin(phi)));
            d.SetZ(d0.Z() * cos(theta) - sqrt(1 - pow(d0.Z(), 2)) * sin(theta) * cos(phi));
        }
        d0 = d;

    } while (E != 0); // Loop until neutron's energy is fully deposited
    return cnt;
}


/**
 * @brief Performs neutron simulation within a scintillator detector.
 *
 * This function initializes the random number generator, draws the geometry of the detector, and creates a histogram to 
 * store the deposited energy of neutrons. It iterates over a specified number of neutrons, simulating their trajectories 
 * within the detector. For each neutron, it generates a random initial position and direction, simulates its trajectory 
 * using the Count function, updates the histogram with the deposited energy, and optionally visualizes the neutron trajectory.
 *
 * @note This function requires the global variables `rnd`, `c`, `E`, `nneu`, `ngr`, `detr`, `det`, `dz`, `A`, `sige` to be properly initialized.
 */
void SimNeu() {
    // Initialize random number generator
    rnd.SetSeed(time(NULL));

    // Draw geometry and create histogram
    DrawGeom();
    TH1D *hE = new TH1D("Edep", "", 100, 0, 0.2);
    c->cd(3);
    hE->Draw();

    // Iterate over neutrons
    int ndet = 0;
    for (int in = 0; in < nneu; in++) {
        // Generate random initial position and direction
        TVector3 d0(0, 0, 1);
        double r02 = rnd.Rndm() * pow(detr, 2);
        double r0 = sqrt(r02);
        double th0 = rnd.Rndm() * 2 * TMath::Pi();
        TVector3 x0(r0 * cos(th0), r0 * sin(th0), 0);

        // Create trajectory graphs if requested
        TGraph *grxy, *grzy;
        bool fillgraph = false;
        if (in <= ngr) {
            grxy = new TGraph();
            grxy->SetMarkerColor(in + 1);
            grxy->SetMarkerStyle(20);
            grxy->SetMarkerSize(0.40);
            grzy = new TGraph(*grxy);
            grzy->SetMarkerColor(in + 1);
            grzy->SetMarkerStyle(20);
            grzy->SetMarkerSize(0.40);
            grxy->SetPoint(grzy->GetN(), x0.X(), x0.Y());
            grzy->SetPoint(grzy->GetN(), x0.Z(), x0.Y());
            fillgraph = true;
        }

        // Simulate neutron trajectory and update histogram
        double Eterm = E;
        bool cnt = Count(x0, d0, Eterm, grxy, grzy, fillgraph);
        if (cnt) {
            ndet++;
        }
        if (fillgraph) {
            c->cd(1);
            grxy->Draw("LP");
            c->cd(2);
            grzy->Draw("LP");
        }
        hE->Fill(Eterm);
    }

    // Print measured flux
    cout << ndet << endl;
    cout << "measured flux " << ndet / (detr * detr * TMath::Pi()) << endl;
}


#ifndef __CINT__
// Main function
int main() {
    TApplication app("app", 0, NULL);
    SimNeu();
    gPad->Update();
    app.Run(true);
}
#endif
