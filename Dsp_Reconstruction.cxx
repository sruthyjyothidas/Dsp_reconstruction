// Λc⁺ → p K⁻ π⁺ reconstruction with (6.28 ± 0.32)%
// Shyam Kumar; INFN Bari, Italy
// shyam.kumar@ba.infn.it; shyam055119@gmail.com 

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class PlotFile;
#endif

#ifndef __CINT__
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include "math.h"
#include "string.h"

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TSystem.h"
#include "TUnixSystem.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLatex.h"
#include "TMinuit.h"
#include "Math/Functor.h"
#include "Fit/Fitter.h"
#include "Math/Minimizer.h"
#endif
#include "StPhysicalHelix.h"
#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

using namespace std;

const double gPionMass = 0.13957;
const double gKaon1Mass = 0.493677;
//const double gProtonMass = 0.938272;
const double gKaon2Mass = 0.493677;

const double twoPi = 2.*3.1415927;
const double eMass = 0.000511;
//const int pdg_Lcp = 4122, pdg_p = 2212, pdg_k = 321, pdg_pi = 211;
const int pdg_Dsp = 431, pdg_kp = 321, pdg_km = 321, pdg_pi = 211;
bool debug = false;


const double bField = -1.7; // Tesla

void getDecayVertex_Chi2fit(const int index1, const int index2, const int index3, double &s1, double &s2, double &s3, TVector3 &vertex, float &chi2, double *parFitErr);
TVector3 GetDCAToPrimaryVertex(const int index, TVector3 vtx);

TLorentzVector GetDCADaughters(const int index1, const int index2, const int index3, TVector3 vtx,
  float *dcaDaughters, float &cosTheta, float &cosTheta_xy, float &decayLength, float &V0DcaToVtx, float &sigma_vtx, float & chi2_ndf, TVector3 &decayVertex, float &pT_kp, float &pT_km, float &pT_pi, double *parFitErr, float &px_kp, float &px_km, float &px_pi, float &py_kp, float &py_km, float &py_pi, float &pz_kp, float &pz_km, float &pz_pi);

TTreeReaderArray<float> *rcMomPx2;
TTreeReaderArray<float> *rcMomPy2;
TTreeReaderArray<float> *rcMomPz2;
TTreeReaderArray<float> *rcCharge2;

TTreeReaderArray<float> *rcTrkLoca2;
TTreeReaderArray<float> *rcTrkLocb2;
TTreeReaderArray<float> *rcTrkTheta2;
TTreeReaderArray<float> *rcTrkPhi2;
TTreeReaderArray<std::array<float, 21>> *rcTrkCov;

// Define function for Chi2 minimization between two helices	
struct Chi2Minimization {    
  StPhysicalHelix fhelix1, fhelix2, fhelix3;
  std::array<float, 21> fcov1, fcov2, fcov3; // full covariance matrix

  Chi2Minimization(StPhysicalHelix helix1, StPhysicalHelix helix2, StPhysicalHelix helix3, std::array<float, 21> cov1, std::array<float, 21> cov2, std::array<float, 21> cov3) : fhelix1(helix1),fhelix2(helix2), fhelix3(helix3), fcov1(cov1), fcov2(cov2), fcov3(cov3) {}
  
  // Implementation of the function to be minimized
  double operator() (const double *par) {
    double x = par[0];
    double y = par[1];
    double z = par[2];
    double s1 = par[3];
    double s2 = par[4];
    double s3 = par[5];    
    double f = 0;
    TVector3 vertex(x, y, z);
    TVector3 p1 = fhelix1.at(s1);
    TVector3 p2 = fhelix2.at(s2);
    TVector3 p3 = fhelix3.at(s3);
    TVector3 mom1 = fhelix1.momentumAt(s1,  bField * tesla);
    TVector3 mom2 = fhelix2.momentumAt(s2,  bField * tesla);
    TVector3 mom3 = fhelix3.momentumAt(s3,  bField * tesla);
    
    // x= −l0 sinϕ , y=l0 cosϕ , z=l
    // Recalculate l0 at PCA for error propagation
    float l0_track1 = p1.Pt(); float l1_track1 = p1.Z(); double phi_track1 = mom1.Phi();
    float l0_track2 = p2.Pt(); float l1_track2 = p2.Z(); double phi_track2 = mom2.Phi();
    float l0_track3 = p3.Pt(); float l1_track3 = p3.Z(); double phi_track3 = mom3.Phi();  
    // Track1: σx^2​=sin^2ϕ⋅σℓ0​^2​+ℓ0^2​cos^2ϕ⋅σϕ^2​+2⋅ℓ0​sinϕcosϕ⋅Cov(ℓ0​,ϕ)
    float sigx1_2 = sin(phi_track1)*sin(phi_track1)*fcov1[0] + l0_track1*l0_track1*cos(phi_track1)*cos(phi_track1)*fcov1[5]+ 2.0*l0_track1*sin(phi_track1)*cos(phi_track1)*fcov1[3]; 
    float sigx2_2 = sin(phi_track2)*sin(phi_track2)*fcov2[0] + l0_track2*l0_track2*cos(phi_track2)*cos(phi_track2)*fcov2[5]+ 2.0*l0_track2*sin(phi_track2)*cos(phi_track2)*fcov2[3];
    float sigx3_2 = sin(phi_track3)*sin(phi_track3)*fcov3[0] + l0_track3*l0_track3*cos(phi_track3)*cos(phi_track3)*fcov3[5]+ 2.0*l0_track3*sin(phi_track3)*cos(phi_track3)*fcov3[3];
  
    // σy^2​=cos^2ϕ⋅σℓ0^​2​+ℓ0^2​sin^2ϕ⋅σϕ^2​−2⋅ℓ0​sinϕcosϕ⋅Cov(ℓ0​,ϕ)
    float sigy1_2 = cos(phi_track1)*cos(phi_track1)*fcov1[0] + l0_track1*l0_track1*sin(phi_track1)*sin(phi_track1)*fcov1[5]-2.0*l0_track1*sin(phi_track1)*cos(phi_track1)*fcov1[3]; 
    float sigy2_2 = cos(phi_track2)*cos(phi_track2)*fcov2[0] + l0_track2*l0_track2*sin(phi_track2)*sin(phi_track2)*fcov2[5]-2.0*l0_track2*sin(phi_track2)*cos(phi_track2)*fcov2[3]; 
    float sigy3_2 = cos(phi_track3)*cos(phi_track3)*fcov3[0] + l0_track3*l0_track3*sin(phi_track3)*sin(phi_track3)*fcov3[5]-2.0*l0_track3*sin(phi_track3)*cos(phi_track3)*fcov3[3]; 
     
    // σz^2​
    float sigz1_2 = fcov1[2];
    float sigz2_2 = fcov2[2]; 
    float sigz3_2 = fcov3[2];   
    double d1_x = 10.*(vertex - p1).X(); double d2_x = 10.*(vertex - p2).X(); double d3_x = 10.*(vertex - p3).X();
    double d1_y = 10.*(vertex - p1).Y(); double d2_y = 10.*(vertex - p2).Y(); double d3_y = 10.*(vertex - p3).Y();
    double d1_z = 10.*(vertex - p1).Z(); double d2_z = 10.*(vertex - p2).Z(); double d3_z = 10.*(vertex - p3).Z();
    
    f = d1_x*d1_x/sigx1_2 + d2_x*d2_x/sigx2_2 + d3_x*d3_x/sigx3_2 + d1_y*d1_y/sigy1_2 + d2_y*d2_y/sigy2_2 + d3_y*d3_y/sigy3_2 + d1_z*d1_z/sigz1_2+ d2_z*d2_z/sigz2_2 + d3_z*d3_z/sigz3_2; // chi2 

    return f;
	}
};


// Define function for Distance2 minimization between three helices	
struct Distance2Minimization {    
     StPhysicalHelix fhelix1, fhelix2, fhelix3;
     Distance2Minimization(StPhysicalHelix helix1, StPhysicalHelix helix2, StPhysicalHelix helix3) : fhelix1(helix1),fhelix2(helix2), fhelix3(helix3) {}
    // Implementation of the function to be minimized
    double operator() (const double *par) {
    double x = par[0];
    double y = par[1];
    double z = par[2];
    double s1 = par[3];
    double s2 = par[4]; 
    double s3 = par[5];   
    double f = 0;
    
    TVector3 vertex(x, y, z);
    TVector3 p1 = fhelix1.at(s1);
    TVector3 p2 = fhelix2.at(s2);
    TVector3 p3 = fhelix3.at(s3);

    double d1 = (vertex - p1).Mag2();
    double d2 = (vertex - p2).Mag2();
    double d3 = (vertex - p3).Mag2();
    f = d1 + d2 + d3; // minimize total squared distance
    return f;	
	}
};

int main(int argc, char **argv)
{
  TString listname = "DISepInputFiles_q2_1to10_L501To1000.list";

  TChain *chain = new TChain("events");

  int nfiles = 0;
  char filename[512];
  ifstream *inputstream = new ifstream;
  inputstream->open(listname.Data());
  
  if(!inputstream){
    printf("[e] Cannot open file list: %s\n", listname.Data());
  }
  
  while(inputstream->good()){
    inputstream->getline(filename, 512);
    
    if (inputstream->good()){
      TFile *ftmp = TFile::Open(filename, "read");
      if (!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())) 
      {
        printf("[e] Could you open file: %s\n", filename);
      }else{
        cout<<"[i] Add "<<nfiles<<"th file: "<<filename<<endl;
        chain->Add(filename);
        nfiles++;
      }
    }
  }

  inputstream->close();
  printf("[i] Read in %d files with %lld events in total\n", nfiles, chain->GetEntries());

  TH1F *hEventStat = new TH1F("hEventStat", "Event statistics", 11, 0, 11);
  hEventStat->GetXaxis()->SetBinLabel(1, "MC events");
  hEventStat->GetXaxis()->SetBinLabel(2, "D_{s}^{+}");
  hEventStat->GetXaxis()->SetBinLabel(3, "D_{s}^{+} -> #Phi + #pi^{+}");
  hEventStat->GetXaxis()->SetBinLabel(4, "Reco Signal D_{s}^{+}/D_{s}^{-}");
  hEventStat->GetXaxis()->SetBinLabel(5, "Reco Signal D_{s}^{+}");
  hEventStat->GetXaxis()->SetBinLabel(6, "Reco Signal D_{s}^{-}");
  hEventStat->GetXaxis()->SetBinLabel(7, "D_{s}^{+} -> K^{+} + K^{-} + #pi^{+}");
  hEventStat->GetXaxis()->SetBinLabel(8, "Reco Signal D_{s}^{+}/D_{s}^{-}");
  hEventStat->GetXaxis()->SetBinLabel(9, "Reco Signal D_{s}^{+}");
  hEventStat->GetXaxis()->SetBinLabel(10, "Reco Signal D_{s}^{-}");
  hEventStat->GetXaxis()->SetBinLabel(11, "Reco Bkg D_{s}^{+}");

  TH1F *hMcMult = new TH1F("hMcMult", "MC multiplicity (|#eta| < 3.5);N_{MC}", 50, 0, 50);

  TH1F *hMcVtxX = new TH1F("hMcVtxX", "x position of MC vertex;x (mm)", 500, -5.0, 5.0);
  TH1F *hMcVtxY = new TH1F("hMcVtxY", "y position of MC vertex;y (mm)", 500, -5.0, 5.0);
  TH1F *hMcVtxZ = new TH1F("hMcVtxZ", "z position of MC vertex;z (mm)", 800, -200, 200);

  TH1F *hRecVtxX = new TH1F("hRecVtxX", "x position of Rec vertex;x (mm)", 500, -5.0, 5.0);
  TH1F *hRecVtxY = new TH1F("hRecVtxY", "y position of Rec vertex;y (mm)", 500, -5.0, 5.0);
  TH1F *hRecVtxZ = new TH1F("hRecVtxZ", "z position of Rec vertex;z (mm)", 800, -200, 200);

  TH1F *hPullVtxX = new TH1F("hPullVtxX", "Pull x position of MC vertex;(Vx_{rec}-Vx_{mc})/#sigma_{vx}; Entries (a.u.)", 1000, -10., 10.);
  TH1F *hPullVtxY = new TH1F("hPullVtxY", "Pull y position of MC vertex;(Vy_{rec}-Vy_{mc})/#sigma_{vy}; Entries (a.u.)", 1000, -10., 10.);
  TH1F *hPullVtxZ = new TH1F("hPullVtxZ", "Pull z position of MC vertex;(Vz_{rec}-Vz_{mc})/#sigma_{vz}; Entries (a.u.)", 2000, -20, 20);
 
  
  // Secondary vertex with chi2 fit method
  TH3F *hRes_SVx_Helixfit = new TH3F("hRes_SVx_Helixfit", "Fit method: Residual of SV (X); p_{T} (GeV/c); y ; SVx_{rec}-SVx_{mc} (mm)", 500, 0.0, 50.0, 80, -4.0, 4.0, 1000, -5.0, 5.0); 
  TH3F *hRes_SVy_Helixfit = new TH3F("hRes_SVy_Helixfit", "Fit method: Residual of SV (Y); p_{T} (GeV/c); y ; SVy_{rec}-SVy_{mc} (mm)", 500, 0.0, 50.0, 80, -4.0, 4.0, 1000, -5.0, 5.0); 
  TH3F *hRes_SVz_Helixfit = new TH3F("hRes_SVz_Helixfit", "Fit method: Residual of SV (Z); p_{T} (GeV/c); y ; SVz_{rec}-SVz_{mc} (mm)", 500, 0.0, 50.0, 80, -4.0, 4.0, 1000, -5.0, 5.0);   
 
  
  TH3F *hRes_SVx_Helixfit_pull = new TH3F("hRes_SVx_Helixfit_pull", "Fit method: Pull of SV (X); p_{T} (GeV/c); y ; SVx_{rec}-SVx_{mc} (mm)", 500, 0.0, 50.0, 80, -4.0, 4.0, 1000, -5.0, 5.0); 
  TH3F *hRes_SVy_Helixfit_pull = new TH3F("hRes_SVy_Helixfit_pull", "Fit method: Pull of SV (Y); p_{T} (GeV/c); y ; SVy_{rec}-SVy_{mc} (mm)", 500, 0.0, 50.0, 80, -4.0, 4.0, 1000, -5.0, 5.0); 
  TH3F *hRes_SVz_Helixfit_pull = new TH3F("hRes_SVz_Helixfit_pull", "Fit method: Pull of SV (Z); p_{T} (GeV/c); y ; SVz_{rec}-SVz_{mc} (mm)", 500, 0.0, 50.0, 80, -4.0, 4.0, 1000, -5.0, 5.0); 
  
  TH1F *hchi2_vtx = new TH1F("hchi2_vtx", "Helix Calculation: Chi2/ndf; #chi^{2}/ndf; Entries (a.u.)", 1000, 0.0, 50.0); 
  TH1F *hchi2_vtx_sig = new TH1F("hchi2_vtx_sig", "Helix Calculation: Chi2/ndf; #chi^{2}/ndf; Entries (a.u.)", 1000, 0.0, 50.0); 
  TH1F *hchi2_vtx_bkg = new TH1F("hchi2_vtx_bkg", "Helix Calculation: Chi2/ndf; #chi^{2}/ndf; Entries (a.u.)", 1000, 0.0, 50.0); 

  TH2F *hDspDecayVxVy = new TH2F("hDspDecayVxVy", "D_{s}^{+} decay vertex to primary vertex;#Deltav_{x} (mm);#Deltav_{y} (mm)", 400, -1-0.0025, 1-0.0025, 400, -1-0.0025, 1-0.0025);
  TH2F *hDspDecayVrVz = new TH2F("hDspDecayVrVz", "D_{s}^{+} decay vertex to primary vertex;#Deltav_{z} (mm);#Deltav_{r} (mm)", 100, -2, 2, 100, -0.2, 1.8);

  TH2F *hMCDspPtRap = new TH2F("hMCDspPtRap", "MC D_{s}^{+};y;p_{T} (GeV/c)", 20, -5, 5, 100, 0, 10);

  TH2F *hMcKpPtEta = new TH2F("hMcKpPtEta", "MC P from D_{s}^{+} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);
  TH2F *hMcKpPtEtaReco = new TH2F("hMcKpPtEtaReco", "RC P from D_{s}^{+} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);

  TH2F *hMcKmPtEta = new TH2F("hMcKmPtEta", "MC K from D_{s}^{+} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);
  TH2F *hMcKmPtEtaReco = new TH2F("hMcKmPtEtaReco", "RC K from D_{s}^{+} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);

  TH2F *hMcPiPtEta = new TH2F("hMcPiPtEta", "MC #pi from D_{s}^{+} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);
  TH2F *hMcPiPtEtaReco = new TH2F("hMcPiPtEtaReco", "RC #pi from D_{s}^{+} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);

  TH1F *hNRecoVtx = new TH1F("hNRecoVtx", "Number of reconstructed vertices;N", 10, 0, 10);

  const char* part_name[3] = {"P", "K", "Pi"};
  const char* part_title[3] = {"P", "K", "#pi"};
  TH3F *hRcSecPartLocaToRCVtx[3];
  TH3F *hRcSecPartLocbToRCVtx[3];
  TH3F *hRcPrimPartLocaToRCVtx[3];
  TH3F *hRcPrimPartLocbToRCVtx[3];
 
  for (int i=0; i<3; i++){
    hRcSecPartLocaToRCVtx[i] = new TH3F(Form("hRcSec%sLocaToRCVtx",part_name[i]), Form( "DCA_{xy} distribution for D_{s}^{+} decayed %s;p_{T} (GeV/c);#eta;DCA_{xy} (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
    hRcSecPartLocbToRCVtx[i] = new TH3F(Form("hRcSec%sLocbToRCVtx",part_name[i]), Form( "DCA_{z} distribution for D_{s}^{+} decayed %s;p_{T} (GeV/c);#eta;DCA_{z} (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, -0.5, 0.5);
    hRcPrimPartLocaToRCVtx[i] = new TH3F(Form("hRcPrim%sLocaToRCVtx",part_name[i]), Form( "DCA_{xy} distribution for primary %s;p_{T} (GeV/c);#eta;DCA_{xy} (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
    hRcPrimPartLocbToRCVtx[i] = new TH3F(Form("hRcPrim%sLocbToRCVtx",part_name[i]), Form( "DCA_{z} distribution for primary %s;p_{T} (GeV/c);#eta;DCA_{z} (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, -0.5, 0.5);
  }

  const char* axis_name[3] = {"x", "y", "z"};
  const int nDimDca = 4;
  const int nBinsDca[nDimDca] = {50, 20, 500, 50};
  const double minBinDca[nDimDca] = {0, -5, -1+0.002, 0};
  const double maxBinDca[nDimDca] = {5, 5, 1+0.002, 50};
  THnSparseF *hPrimTrkDcaToRCVtx[3][3];

  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      hPrimTrkDcaToRCVtx[i][j] = new THnSparseF(Form("hPrim%sDca%sToRCVtx",part_name[i],axis_name[j]), Form("DCA_{%s} distribution for primary %s;p_{T} (GeV/c);#eta;DCA_{%s} (mm);N_{MC}",axis_name[j],part_title[i],axis_name[j]), nDimDca, nBinsDca, minBinDca, maxBinDca);
    }
  }

  TH3F *h3PairDca12[2], *h3PairDca23[2], *h3PairDca13[2];
  TH3F *h3PairCosTheta[2];
  TH3F *h3PairDca[2];
  TH3F *h3PairDecayLength[2];
  const char* pair_name[2] = {"signal", "bkg"};
  const char* pair_title[2] = {"Signal", "Background"};
  
  for (int i=0; i<2; i++){
    h3PairDca12[i] = new TH3F(Form("h3PairDca12_%s", pair_name[i]), Form("%s pair DCA_{12};p_{T} (GeV/c);#eta;DCA_{12} (mm)", pair_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
    h3PairDca23[i] = new TH3F(Form("h3PairDca23_%s", pair_name[i]), Form("%s pair DCA_{23};p_{T} (GeV/c);#eta;DCA_{23} (mm)", pair_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
    h3PairDca13[i] = new TH3F(Form("h3PairDca13_%s", pair_name[i]), Form("%s pair DCA_{13};p_{T} (GeV/c);#eta;DCA_{13} (mm)", pair_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
    h3PairCosTheta[i] = new TH3F(Form("h3PairCosTheta_%s", pair_name[i]), Form("%s pair cos(#theta);p_{T} (GeV/c);#eta;cos(#theta)", pair_title[i]), 100, 0, 10, 20, -5, 5, 100, -1, 1);
  
    h3PairDca[i] = new TH3F(Form("h3PairDca_%s", pair_name[i]), Form("%s pair DCA;p_{T} (GeV/c);#eta;DCA_{pair} (mm)", pair_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);

    h3PairDecayLength[i] = new TH3F(Form("h3PairDecayLength_%s", pair_name[i]), Form("%s pair decay length;p_{T} (GeV/c);#eta;L (mm)", pair_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
  } 

  // Invariant mass
  const char* cut_name[2] = {"all", "DCA"};
  TH3F *h3InvMass[2][2];

  for(int i=0; i<2; i++){
    for(int j=0; j<2; j++){
      //h3InvMass[i][j] = new TH3F(Form("h3InvMass_%s_%s", pair_name[i], cut_name[j]), "Invariant mass of unlike-sign #piK pairs;p_{T} (GeV/c);y;M_{#piK} (GeV/c^{2})", 100, 0, 10, 20, -5, 5, 100, 2.0, 3.0);
      h3InvMass[i][j] = new TH3F(Form("h3InvMass_%s_%s", pair_name[i], cut_name[j]), "Invariant mass of unlike-sign #piK pairs;p_{T} (GeV/c);y;M_{#piK} (GeV/c^{2})", 100, 0, 10, 20, -5, 5, 100, 1.5, 2.5);
    }
  }

  TTreeReader treereader(chain);
  // MC  
  TTreeReaderArray<int> mcPartGenStatus = {treereader, "MCParticles.generatorStatus"};
  TTreeReaderArray<int> mcPartPdg = {treereader, "MCParticles.PDG"};
  TTreeReaderArray<float> mcPartCharge = {treereader, "MCParticles.charge"};
  TTreeReaderArray<unsigned int> mcPartParent_begin = {treereader, "MCParticles.parents_begin"};
  TTreeReaderArray<unsigned int> mcPartParent_end = {treereader, "MCParticles.parents_end"};
  TTreeReaderArray<int> mcPartParent_index = {treereader, "_MCParticles_parents.index"};
  TTreeReaderArray<unsigned int> mcPartDaughter_begin = {treereader, "MCParticles.daughters_begin"};
  TTreeReaderArray<unsigned int> mcPartDaughter_end = {treereader, "MCParticles.daughters_end"};
  TTreeReaderArray<int> mcPartDaughter_index = {treereader, "_MCParticles_daughters.index"};
  TTreeReaderArray<double> mcPartMass = {treereader, "MCParticles.mass"};
  TTreeReaderArray<double> mcPartVx = {treereader, "MCParticles.vertex.x"};
  TTreeReaderArray<double> mcPartVy = {treereader, "MCParticles.vertex.y"};
  TTreeReaderArray<double> mcPartVz = {treereader, "MCParticles.vertex.z"};
  TTreeReaderArray<double> mcMomPx = {treereader, "MCParticles.momentum.x"};
  TTreeReaderArray<double> mcMomPy = {treereader, "MCParticles.momentum.y"};
  TTreeReaderArray<double> mcMomPz = {treereader, "MCParticles.momentum.z"};
  TTreeReaderArray<double> mcEndPointX = {treereader, "MCParticles.endpoint.x"};
  TTreeReaderArray<double> mcEndPointY = {treereader, "MCParticles.endpoint.y"};
  TTreeReaderArray<double> mcEndPointZ = {treereader, "MCParticles.endpoint.z"};

  TTreeReaderArray<int> assocChSimID = {treereader, "_ReconstructedChargedParticleAssociations_sim.index"};
  TTreeReaderArray<int> assocChRecID = {treereader, "_ReconstructedChargedParticleAssociations_rec.index"};
  TTreeReaderArray<float> assocWeight = {treereader, "ReconstructedChargedParticleAssociations.weight"};

  TTreeReaderArray<float> rcMomPx = {treereader, "ReconstructedChargedParticles.momentum.x"};
  TTreeReaderArray<float> rcMomPy = {treereader, "ReconstructedChargedParticles.momentum.y"};
  TTreeReaderArray<float> rcMomPz = {treereader, "ReconstructedChargedParticles.momentum.z"};
  TTreeReaderArray<float> rcPosx = {treereader, "ReconstructedChargedParticles.referencePoint.x"};
  TTreeReaderArray<float> rcPosy = {treereader, "ReconstructedChargedParticles.referencePoint.y"};
  TTreeReaderArray<float> rcPosz = {treereader, "ReconstructedChargedParticles.referencePoint.z"};
  TTreeReaderArray<float> rcCharge = {treereader, "ReconstructedChargedParticles.charge"};
  TTreeReaderArray<int>   rcPdg = {treereader, "ReconstructedChargedParticles.PDG"};

  TTreeReaderArray<float> rcTrkLoca = {treereader, "CentralCKFTrackParameters.loc.a"};
  TTreeReaderArray<float> rcTrkLocb = {treereader, "CentralCKFTrackParameters.loc.b"};
  TTreeReaderArray<float> rcTrkqOverP = {treereader, "CentralCKFTrackParameters.qOverP"};
  TTreeReaderArray<float> rcTrkTheta = {treereader, "CentralCKFTrackParameters.theta"};
  TTreeReaderArray<float> rcTrkPhi = {treereader, "CentralCKFTrackParameters.phi"};

  rcMomPx2 = new TTreeReaderArray<float>{treereader, "ReconstructedChargedParticles.momentum.x"};
  rcMomPy2 = new TTreeReaderArray<float>{treereader, "ReconstructedChargedParticles.momentum.y"};
  rcMomPz2 = new TTreeReaderArray<float>{treereader, "ReconstructedChargedParticles.momentum.z"};
  rcCharge2 = new TTreeReaderArray<float>{treereader, "ReconstructedChargedParticles.charge"};

  rcTrkLoca2 = new TTreeReaderArray<float>{treereader, "CentralCKFTrackParameters.loc.a"};
  rcTrkLocb2 = new TTreeReaderArray<float>{treereader, "CentralCKFTrackParameters.loc.b"};
  rcTrkTheta2 = new TTreeReaderArray<float>{treereader, "CentralCKFTrackParameters.theta"};
  rcTrkPhi2 = new TTreeReaderArray<float>{treereader, "CentralCKFTrackParameters.phi"};
  rcTrkCov = new TTreeReaderArray<std::array<float, 21>>{treereader, "CentralCKFTrackParameters.covariance.covariance[21]"};

  TTreeReaderArray<float> CTVx = {treereader, "CentralTrackVertices.position.x"};
  TTreeReaderArray<float> CTVy = {treereader, "CentralTrackVertices.position.y"};
  TTreeReaderArray<float> CTVz = {treereader, "CentralTrackVertices.position.z"};
  TTreeReaderArray<int> CTVndf = {treereader, "CentralTrackVertices.ndf"};
  TTreeReaderArray<float> CTVchi2 = {treereader, "CentralTrackVertices.chi2"};
  TTreeReaderArray<float> CTVerr_xx = {treereader, "CentralTrackVertices.positionError.xx"};
  TTreeReaderArray<float> CTVerr_yy = {treereader, "CentralTrackVertices.positionError.yy"};
  TTreeReaderArray<float> CTVerr_zz = {treereader, "CentralTrackVertices.positionError.zz"};
  TTreeReaderArray<float> incQ2 = {treereader, "InclusiveKinematicsTruth.Q2"};
  TTreeReaderArray<float> incxB = {treereader, "InclusiveKinematicsTruth.x"};

  TTreeReaderArray<int> prim_vtx_index = {treereader, "PrimaryVertices_objIdx.index"};

  TTreeReaderArray<unsigned int> vtxAssocPart_begin = {treereader, "CentralTrackVertices.associatedParticles_begin"};
  TTreeReaderArray<unsigned int> vtxAssocPart_end = {treereader, "CentralTrackVertices.associatedParticles_end"};
  TTreeReaderArray<int> vtxAssocPart_index = {treereader, "_CentralTrackVertices_associatedParticles.index"};

  //TFile *outfile = new TFile("DsNtuple_eAu.root", "RECREATE");
  //TFile *outfile = new TFile("DsNtuple_ep_2000To2399.root", "RECREATE");
  TFile *outfile = new TFile("DISNtuple_ep_q2_1to10_L501To1000.root", "RECREATE");
  outfile->cd();

  // Create a ROOT file to store the Ntuple
  //TFile *file_signal = new TFile(signal.Data(), "RECREATE");
  TTree *tree_sig = new TTree("treeMLSig", "treeMLSig"); 

  // Define variables to store in the Ntuple
  float d0_kp_sig, d0_km_sig, d0_pi_sig, pT_kp_sig, pT_km_sig, pT_pi_sig, px_kp_sig, px_km_sig, px_pi_sig, py_kp_sig, py_km_sig, py_pi_sig, pz_kp_sig, pz_km_sig, pz_pi_sig, d0xy_kp_sig, d0xy_km_sig, d0xy_pi_sig, sum_d0xy_sig, dca_12_sig, dca_Dsp_sig, decay_length_sig, xB_sig, Q2_sig;
  float costheta_sig, costhetaxy_sig, pt_Dsp_sig, y_Dsp_sig, mass_Dsp_sig, sigma_vtx_sig, mult_sig, chi2_sig;
  int decay;

  // Link the variables to the TTree branches
  tree_sig->Branch("d0_kp", &d0_kp_sig, "d0_kp/F");
  tree_sig->Branch("d0_km", &d0_km_sig, "d0_km/F");
  tree_sig->Branch("d0_pi", &d0_pi_sig, "d0_pi/F"); 
  tree_sig->Branch("pT_kp", &pT_kp_sig, "pT_kp/F");
  tree_sig->Branch("pT_km", &pT_km_sig, "pT_km/F");
  tree_sig->Branch("pT_pi", &pT_pi_sig, "pT_pi/F"); 
  tree_sig->Branch("px_kp", &px_kp_sig, "px_kp/F");
  tree_sig->Branch("px_km", &px_km_sig, "px_km/F");
  tree_sig->Branch("px_pi", &px_pi_sig, "px_pi/F"); 
  tree_sig->Branch("py_kp", &py_kp_sig, "py_kp/F");
  tree_sig->Branch("py_km", &py_km_sig, "py_km/F");
  tree_sig->Branch("py_pi", &py_pi_sig, "py_pi/F"); 
  tree_sig->Branch("pz_kp", &pz_kp_sig, "pz_kp/F");
  tree_sig->Branch("pz_km", &pz_km_sig, "pz_km/F");
  tree_sig->Branch("pz_pi", &pz_pi_sig, "pz_pi/F"); 
  tree_sig->Branch("d0xy_kp", &d0xy_kp_sig, "d0xy_kp/F");
  tree_sig->Branch("d0xy_km", &d0xy_km_sig, "d0xy_km/F");
  tree_sig->Branch("d0xy_pi", &d0xy_pi_sig, "d0xy_pi/F");  
  tree_sig->Branch("sum_d0xy", &sum_d0xy_sig, "sum_d0xy/F");        
  tree_sig->Branch("dca_12", &dca_12_sig, "dca_12/F");
  tree_sig->Branch("dca_Dsp", &dca_Dsp_sig, "dca_Dsp/F");
  tree_sig->Branch("pt_Dsp", &pt_Dsp_sig, "pt_Dsp/F");  
  tree_sig->Branch("y_Dsp", &y_Dsp_sig, "y_Dsp/F");  
  tree_sig->Branch("mass_Dsp", &mass_Dsp_sig, "mass_Dsp/F");              
  tree_sig->Branch("decay_length", &decay_length_sig, "decay_length/F");   
  tree_sig->Branch("costheta", &costheta_sig, "costheta/F"); 
  tree_sig->Branch("costheta_xy", &costhetaxy_sig, "costheta_xy/F"); 
  tree_sig->Branch("sigma_vtx", &sigma_vtx_sig, "sigma_vtx/F"); 
  tree_sig->Branch("mult", &mult_sig, "mult/F");
  tree_sig->Branch("chi2", &chi2_sig, "chi2/F");
  tree_sig->Branch("xB", &xB_sig, "xB/F");
  tree_sig->Branch("Q2", &Q2_sig, "Q2/F");                          
  tree_sig->Branch("decay", &decay, "decay/I");                          

  //TFile *file_bkg = new TFile(bkg.Data(), "RECREATE");
  TTree *tree_bkg = new TTree("treeMLBkg", "treeMLBkg"); 

  // Define variables to store in the Ntuple
  float d0_kp_bkg, d0_km_bkg, d0_pi_bkg, pT_kp_bkg, pT_km_bkg, pT_pi_bkg, px_kp_bkg, px_km_bkg, px_pi_bkg, py_kp_bkg, py_km_bkg, py_pi_bkg, pz_kp_bkg, pz_km_bkg, pz_pi_bkg, d0xy_kp_bkg, d0xy_km_bkg, d0xy_pi_bkg, sum_d0xy_bkg, dca_12_bkg, dca_Dsp_bkg, decay_length_bkg, xB_bkg, Q2_bkg;
  float costheta_bkg, costhetaxy_bkg, pt_Dsp_bkg, y_Dsp_bkg, mass_Dsp_bkg, sigma_vtx_bkg, mult_bkg, chi2_bkg;
  // Link the variables to the TTree branches
  tree_bkg->Branch("d0_kp", &d0_kp_bkg, "d0_kp/F");
  tree_bkg->Branch("d0_km", &d0_km_bkg, "d0_km/F");
  tree_bkg->Branch("d0_pi", &d0_pi_bkg, "d0_pi/F");
  tree_bkg->Branch("pT_kp", &pT_kp_bkg, "pT_kp/F");
  tree_bkg->Branch("pT_km", &pT_km_bkg, "pT_km/F");
  tree_bkg->Branch("pT_pi", &pT_pi_bkg, "pT_pi/F"); 
  tree_bkg->Branch("px_kp", &px_kp_bkg, "px_kp/F");
  tree_bkg->Branch("px_km", &px_km_bkg, "px_km/F");
  tree_bkg->Branch("px_pi", &px_pi_bkg, "px_pi/F"); 
  tree_bkg->Branch("py_kp", &py_kp_bkg, "py_kp/F");
  tree_bkg->Branch("py_km", &py_km_bkg, "py_km/F");
  tree_bkg->Branch("py_pi", &py_pi_bkg, "py_pi/F"); 
  tree_bkg->Branch("pz_kp", &pz_kp_bkg, "pz_kp/F");
  tree_bkg->Branch("pz_km", &pz_km_bkg, "pz_km/F");
  tree_bkg->Branch("pz_pi", &pz_pi_bkg, "pz_pi/F"); 
  tree_bkg->Branch("d0xy_kp", &d0xy_kp_bkg, "d0xy_kp/F");
  tree_bkg->Branch("d0xy_km", &d0xy_km_bkg, "d0xy_km/F");
  tree_bkg->Branch("d0xy_pi", &d0xy_pi_bkg, "d0xy_pi/F");  
  tree_bkg->Branch("sum_d0xy", &sum_d0xy_bkg, "sum_d0xy/F");            
  tree_bkg->Branch("dca_12", &dca_12_bkg, "dca_12/F");
  tree_bkg->Branch("dca_Dsp", &dca_Dsp_bkg, "dca_Dsp/F");
  tree_bkg->Branch("pt_Dsp", &pt_Dsp_bkg, "pt_Dsp/F");  
  tree_bkg->Branch("y_Dsp", &y_Dsp_bkg, "y_Dsp/F");  
  tree_bkg->Branch("mass_Dsp", &mass_Dsp_bkg, "mass_Dsp/F");              
  tree_bkg->Branch("decay_length", &decay_length_bkg, "decay_length/F");   
  tree_bkg->Branch("costheta", &costheta_bkg, "costheta/F");  
  tree_bkg->Branch("costheta_xy", &costhetaxy_bkg, "costheta_xy/F"); 
  tree_bkg->Branch("sigma_vtx", &sigma_vtx_bkg, "sigma_vtx/F"); 
  tree_bkg->Branch("mult", &mult_bkg, "mult/F");
  tree_bkg->Branch("chi2", &chi2_bkg, "chi2/F"); 
  tree_bkg->Branch("xB", &xB_bkg, "xB/F"); 
  tree_bkg->Branch("Q2", &Q2_bkg, "Q2/F"); 

  //TFile *fout_mcgen = new TFile(signalmc.Data(),"RECREATE");
  TTree *tree_Dsp = new TTree("DsTree","Ds Baryon Kinematics");

  float Dsp_px, Dsp_py, Dsp_pz;
  float Dsp_pt, Dsp_eta, Dsp_y, Dsp_mass;
  float Phi_px, Phi_py, Phi_pz;
  float Phi_pt, Phi_eta, Phi_y, Phi_mass;

  int Dsp_decay;

  tree_Dsp->Branch("px",&Dsp_px,"px/F");
  tree_Dsp->Branch("py",&Dsp_py,"py/F");
  tree_Dsp->Branch("pz",&Dsp_pz,"pz/F");
  tree_Dsp->Branch("pt",&Dsp_pt,"pt/F");
  tree_Dsp->Branch("eta",&Dsp_eta,"eta/F");
  tree_Dsp->Branch("rapidity",&Dsp_y,"rapidity/F");
  tree_Dsp->Branch("mass",&Dsp_mass,"mass/F");                                     
  tree_Dsp->Branch("decay",&Dsp_decay,"decay/I");                                     

  tree_Dsp->Branch("px_phi",&Phi_px,"px_phi/F");
  tree_Dsp->Branch("py_phi",&Phi_py,"py_phi/F");
  tree_Dsp->Branch("pz_phi",&Phi_pz,"pz_phi/F");
  tree_Dsp->Branch("pt_phi",&Phi_pt,"pt_phi/F");
  tree_Dsp->Branch("eta_phi",&Phi_eta,"eta_phi/F");
  tree_Dsp->Branch("rapidity_phi",&Phi_y,"rapidity_phi/F");
  tree_Dsp->Branch("mass_phi",&Phi_mass,"mass_phi/F");                                     

  int nevents = 0;
  //int count = 0;

  while(treereader.Next()){

    if(nevents%1000==0) printf("\nEvent No.-----> %d\n",nevents);
    
    // find MC primary vertex
    float Q2_mc = -1000.;
    float xB_mc = -1000.;
    if (incQ2.GetSize()>0){
      Q2_mc = incQ2[0];
      xB_mc = incxB[0];
    }
    int nMCPart = mcPartMass.GetSize();

    TVector3 vertex_mc(-999., -999., -999.);
    for(int imc=0; imc<nMCPart; imc++){
      // Status 4 describes the incoming electron/proton and the end point of electron or proton is MC Vertex  
      if(mcPartGenStatus[imc] == 4 && mcPartPdg[imc] == 11){
        vertex_mc.SetXYZ(mcEndPointX[imc], mcEndPointY[imc], mcEndPointZ[imc]);
        break;
      }  
    }
    hEventStat->Fill(0.5);
    hMcVtxX->Fill(vertex_mc.x());
    hMcVtxY->Fill(vertex_mc.y());
    hMcVtxZ->Fill(vertex_mc.z());


    // get RC primary vertex
    TVector3 vertex_rc(-999., -999., -999.);
    TVector3 err_vertex_rc(-999., -999., -999.);

    if (prim_vtx_index.GetSize()>0){
      int rc_vtx_index = prim_vtx_index[0];
      vertex_rc.SetXYZ(CTVx[rc_vtx_index], CTVy[rc_vtx_index], CTVz[rc_vtx_index]);
      err_vertex_rc.SetXYZ(sqrt(CTVerr_xx[rc_vtx_index]), sqrt(CTVerr_yy[rc_vtx_index]), sqrt(CTVerr_zz[rc_vtx_index]));
    }
    
    hRecVtxX->Fill(vertex_rc.x());
    hRecVtxY->Fill(vertex_rc.y());
    hRecVtxZ->Fill(vertex_rc.z());

    hPullVtxX->Fill((vertex_rc.x()-vertex_mc.x())/err_vertex_rc.x()); 
    hPullVtxY->Fill((vertex_rc.y()-vertex_mc.y())/err_vertex_rc.y()); 
    hPullVtxZ->Fill((vertex_rc.z()-vertex_mc.z())/err_vertex_rc.z());

    // map MC and RC particles
    int nAssoc = assocChRecID.GetSize();
    map<int, int> assoc_map_mc_to_rc;  // MC --> RC Associations
    map<int, int> assoc_map_rc_to_mc;  // RC --> MC Associations

    for (unsigned int rc_index=0; rc_index<rcMomPx.GetSize(); rc_index++){
	    // loop over the association to find the matched MC particle with largest weight
      double max_weight = 0;
      int matched_mc_index = -1;
      for (int j=0; j<nAssoc; j++){
        if (assocChRecID[j] != rc_index) continue;
        if (assocWeight[j] > max_weight){
          max_weight = assocWeight[j];
          matched_mc_index = assocChSimID[j];
        }
      }

	    // build the map
      assoc_map_mc_to_rc[matched_mc_index] = rc_index;
      assoc_map_rc_to_mc[rc_index] = matched_mc_index;
    }
    
    // MC --> RC associations (kv; key,value for map)
    // std::cout << "MC to RC map:\n";
    // for (const auto& kv : assoc_map_mc_to_rc) {
    //     std::cout << "MC index " << kv.first << " ---> RC index " << kv.second << "\n";
    // }

    // RC --> MC associations
    // std::cout << "\nRC to MC map:\n";
    // for (const auto& kv : assoc_map_rc_to_mc) {
    //    std::cout << "RC index " << kv.first << " --->  MC index " << kv.second << "\n";
    // }


    // Loop over primary particles
    int nMcPart = 0;
    for (int imc=0; imc<nMCPart; imc++){
      if (mcPartGenStatus[imc] == 1 && mcPartCharge[imc] != 0){
        double dist = sqrt( pow(mcPartVx[imc]-vertex_mc.x(),2) + pow(mcPartVy[imc]-vertex_mc.y(),2) + pow(mcPartVz[imc]-vertex_mc.z(),2));      
        if (dist < 1e-4){
		      // count primary charged particles within |eta| < 3.5
          TVector3 mc_mom(mcMomPx[imc], mcMomPy[imc], mcMomPz[imc]);
          double mcEta = mc_mom.PseudoRapidity();
          if(fabs(mcEta) < 3.5) nMcPart++;
        }
      }
    }
    // Primary charged multiplicity at generated level with acceptance cut
    hMcMult->Fill(nMcPart);

    for (int imc=0; imc<nMCPart; imc++){
      if (mcPartGenStatus[imc] == 1 && mcPartCharge[imc] != 0){
        double dist = sqrt( pow(mcPartVx[imc]-vertex_mc.x(),2) + pow(mcPartVy[imc]-vertex_mc.y(),2) + pow(mcPartVz[imc]-vertex_mc.z(),2));      
        if (dist < 1e-4){		  
          // check if the MC particle is reconstructed
          int rc_index = -1;
          if (assoc_map_mc_to_rc.find(imc) != assoc_map_mc_to_rc.end()) rc_index = assoc_map_mc_to_rc[imc];

          if (rc_index>=0){
            TVector3 dcaToVtx = GetDCAToPrimaryVertex(rc_index, vertex_rc);

            int ip = -1;
            if (fabs(mcPartPdg[imc]) == pdg_kp) ip = 0;
            if (fabs(mcPartPdg[imc]) == pdg_km) ip = 1;
            if (fabs(mcPartPdg[imc]) == pdg_pi) ip = 2;
            if (ip>=0){
              TVector3 mom(rcMomPx[rc_index], rcMomPy[rc_index], rcMomPz[rc_index]);
              if(ip<3){
                hRcPrimPartLocaToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.Perp());
                hRcPrimPartLocbToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.z());
              }

              double fill1[] = {mom.Pt(), mom.Eta(), dcaToVtx.x(), nMcPart*1.};
              double fill2[] = {mom.Pt(), mom.Eta(), dcaToVtx.y(), nMcPart*1.};
              double fill3[] = {mom.Pt(), mom.Eta(), dcaToVtx.z(), nMcPart*1.};
              hPrimTrkDcaToRCVtx[ip][0]->Fill(fill1);
              hPrimTrkDcaToRCVtx[ip][1]->Fill(fill2);
              hPrimTrkDcaToRCVtx[ip][2]->Fill(fill3);
            }
          }
        }
      }
    }

    // Look for Ds⁺ → K+ K⁻ π⁺
    //bool hasDs = false;
    vector<int> mc_index_Dsp_kp;
    vector<int> mc_index_Dsp_km;
    vector<int> mc_index_Dsp_pi;
    vector<int> mc_index_Dsp_nr_kp;
    vector<int> mc_index_Dsp_nr_km;
    vector<int> mc_index_Dsp_nr_pi;
    mc_index_Dsp_kp.clear();      
    mc_index_Dsp_km.clear();
    mc_index_Dsp_pi.clear();
    mc_index_Dsp_nr_kp.clear();      
    mc_index_Dsp_nr_km.clear();
    mc_index_Dsp_nr_pi.clear();

    for (int imc=0; imc<nMCPart; imc++){
      if (fabs(mcPartPdg[imc]) != pdg_Dsp) continue;
      hEventStat->Fill(1.5);
      
      // Fill Dsp gen tree with properties
      TLorentzVector Dsp_vec;

      Dsp_px   = mcMomPx[imc];
      Dsp_py   = mcMomPy[imc];
      Dsp_pz   = mcMomPz[imc];
      Dsp_mass = mcPartMass[imc];

      Dsp_vec.SetXYZM(Dsp_px, Dsp_py, Dsp_pz, Dsp_mass);

      Dsp_pt  = Dsp_vec.Pt();
      Dsp_eta = Dsp_vec.Eta();
      Dsp_y   = Dsp_vec.Rapidity();

      tree_Dsp->Fill();

      // printf("MC Particle (Momentum) = (%f, %f, %f) \n",mcMomPx[imc], mcMomPy[imc], mcMomPz[imc]);		  	
      int nDaughters = mcPartDaughter_end[imc]-mcPartDaughter_begin[imc];
      //if (nDaughters!=2) continue;

      if (nDaughters!=2 && nDaughters!=3) continue;
      //cout << mcPartPdg[imc] << " " << nDaughters << endl; 
      bool is_phipi_decay = false;	  
      int phi_dau1_index = -1;
      int phi_dau2_index = -1;
      int phi_index = -1;
      int pion_index = -1;

      if (nDaughters==2){
        int daug_index_1 = mcPartDaughter_index[mcPartDaughter_begin[imc]];
        int daug_index_2 = mcPartDaughter_index[mcPartDaughter_begin[imc]+1];
        int daug_pdg_1 = mcPartPdg[daug_index_1];
        int daug_pdg_2 = mcPartPdg[daug_index_2];
        
        // Check Ds -> phi + pi
        if (fabs(daug_pdg_1) == 333 && fabs(daug_pdg_2) == 211) {
          phi_index  = daug_index_1;
          pion_index = daug_index_2;
        }
        if (fabs(daug_pdg_1) == 211 && fabs(daug_pdg_2) == 333) {
          pion_index  = daug_index_1;
          phi_index  = daug_index_2;
        }
        //cout << daug_pdg_1 << " in index " << phi_index << " and " << daug_pdg_2 << " in index " << pion_index << endl;

        // If Ds -> phi pi, check phi -> K+ K-
        if (phi_index >= 0) {
          int nPhiDaughters = mcPartDaughter_end[phi_index] - mcPartDaughter_begin[phi_index];
      
          if (nPhiDaughters == 2) {
            int phi_dau1 = mcPartDaughter_index[mcPartDaughter_begin[phi_index]];
            int phi_dau2 = mcPartDaughter_index[mcPartDaughter_begin[phi_index] + 1];
            int pdg1 = mcPartPdg[phi_dau1];
            int pdg2 = mcPartPdg[phi_dau2];
            
            bool true_phi_to_KK =
              ((pdg1 == 321 && pdg2 == -321) ||
               (pdg1 == -321 && pdg2 == 321));

            if (true_phi_to_KK) {
              is_phipi_decay = true;
              //cout << "TRUE Ds -> phi pi, phi -> K+ K-" << endl;
            
              phi_dau1_index = phi_dau1;
              phi_dau2_index = phi_dau2;
              //cout << pdg1 << " in index " << phi_dau1_index << " and " << pdg2 << " in index " << phi_dau2_index << endl;
            }
          }
        }
      }
      
      // find Ds+ that decay into KKPi
      bool is_kkpi_decay = false;	 
      int daug_index_1 = -1;
      int daug_index_2 = -1;
      int daug_index_3 = -1;
      int daug_pdg_1 = -1;
      int daug_pdg_2 = -1;
      int daug_pdg_3 = -1;
      if (nDaughters==3){
        daug_index_1 = mcPartDaughter_index[mcPartDaughter_begin[imc]];
        daug_index_2 = mcPartDaughter_index[mcPartDaughter_begin[imc]+1];
        daug_index_3 = mcPartDaughter_index[mcPartDaughter_begin[imc]+2];	  
        daug_pdg_1 = mcPartPdg[daug_index_1];
        daug_pdg_2 = mcPartPdg[daug_index_2];
        daug_pdg_3 = mcPartPdg[daug_index_3];
        
        if (((fabs(daug_pdg_1)==321 && fabs(daug_pdg_2)==321 && fabs(daug_pdg_3)==211)
          || (fabs(daug_pdg_1)==321 && fabs(daug_pdg_2)==211 && fabs(daug_pdg_3)==321)
          || (fabs(daug_pdg_1)==211 && fabs(daug_pdg_2)==321 && fabs(daug_pdg_3)==321))){
          is_kkpi_decay = true;
        }
      }

      //cout << phi_index << " " << pion_index << endl;
      //cout << phi_dau1_index << " " << phi_dau2_index << endl;
      //cout << is_phipi_decay << endl;

      //int pdg1 = mcPartPdg[phi_dau1];
      //int pdg2 = mcPartPdg[phi_dau2];


      /*if (nDaughters == 3 &&
         ((fabs(daug_pdg_1)==pdg_kp && fabs(daug_pdg_2)==pdg_km && fabs(daug_pdg_3)==pdg_pi) 
         || (fabs(daug_pdg_1)==pdg_kp && fabs(daug_pdg_2)==pdg_pi && fabs(daug_pdg_3)==pdg_km)
        //|| (fabs(daug_pdg_1)==pdg_km && fabs(daug_pdg_2)==pdg_kp && fabs(daug_pdg_3)==pdg_pi) 
        //|| (fabs(daug_pdg_1)==pdg_km && fabs(daug_pdg_2)==pdg_pi && fabs(daug_pdg_3)==pdg_kp)
        //|| (fabs(daug_pdg_1)==pdg_pi && fabs(daug_pdg_2)==pdg_km && fabs(daug_pdg_3)==pdg_kp)
        || (fabs(daug_pdg_1)==pdg_pi && fabs(daug_pdg_2)==pdg_kp && fabs(daug_pdg_3)==pdg_km))){
        is_kkpi_decay = true;
      }
      */
      //if (!is_phipi_decay) continue;
      
      
      if (is_phipi_decay || is_kkpi_decay){
        //if (!is_phipi_decay) continue;

        //if (!is_phipi_decay && !is_kkpi_decay) continue;
        //if (is_phipi_decay && is_kkpi_decay) cout << "both decay" << endl;
        // Efficiency for Ds+
        //TLorentzVector mc_part_gen;
        //mc_part_gen.SetXYZM(mcMomPx[imc], mcMomPy[imc], mcMomPz[imc], mcPartMass[imc]);
        
        // Fill the momentum at Gen level Aug3, 2025
        //float mcRap_gen = mc_part_gen.Rapidity();
        //float mcPt_gen = mc_part_gen.Pt();  

        if (is_phipi_decay){
          hEventStat->Fill(2.5);
          int Ds_dau1_pdg = mcPartPdg[phi_dau1_index]; 
          int Ds_dau2_pdg = mcPartPdg[phi_dau2_index]; 
          int Ds_dau3_pdg = mcPartPdg[pion_index]; 
          
          //cout << " look for " << Ds_dau1_pdg << " " << Ds_dau2_pdg << " " << Ds_dau3_pdg << endl;
          
          // Store true K+ and K- MC indices
          if (Ds_dau1_pdg == 321 && Ds_dau2_pdg == -321) {
              mc_index_Dsp_kp.push_back(phi_dau1_index);
              mc_index_Dsp_km.push_back(phi_dau2_index);
          }else{
              mc_index_Dsp_kp.push_back(phi_dau2_index);
              mc_index_Dsp_km.push_back(phi_dau1_index);
          }

          //if (fabs(Ds_dau1_pdg)==321 && fabs(Ds_dau2_pdg)==-321) {
          //  mc_index_Dsp_kp.push_back(phi_dau1_index);
          //  mc_index_Dsp_km.push_back(phi_dau2_index);
          //}
          // Store true pi+ MC index
          mc_index_Dsp_pi.push_back(pion_index);
          //cout << "stored index: kp " << phi_dau1_index << ", km " << phi_dau2_index << ", pi " << pion_index << endl;
        }

        if (is_kkpi_decay){

          hEventStat->Fill(6.5);
          int kp_index = -1;
          int km_index = -1;
          int pi_index = -1;

          int dau_index[3] = {
            daug_index_1,
            daug_index_2,
            daug_index_3
          };
          int dau_pdg[3] = {
            daug_pdg_1,
            daug_pdg_2,
            daug_pdg_3
          };

          for (int i = 0; i < 3; i++) {

            if (dau_pdg[i] == 321) {
              kp_index = dau_index[i];
            }else if (dau_pdg[i] == -321) {
              km_index = dau_index[i];
            }else if (dau_pdg[i] == 211 || dau_pdg[i] == -211) {
              pi_index = dau_index[i];
            }
          }
          if (kp_index >= 0 && km_index >= 0 && pi_index >= 0) {

            mc_index_Dsp_nr_kp.push_back(kp_index);
            mc_index_Dsp_nr_km.push_back(km_index);
            mc_index_Dsp_nr_pi.push_back(pi_index);
          }
          /*if((fabs(daug_pdg_1)==321 && fabs(daug_pdg_2)==321 && fabs(daug_pdg_3)==211)){
            mc_index_Dsp_nr_kp.push_back(daug_index_1);
            mc_index_Dsp_nr_km.push_back(daug_index_2);
            mc_index_Dsp_nr_pi.push_back(daug_index_3);
          }else if((fabs(daug_pdg_1)==321 && fabs(daug_pdg_2)==211 && fabs(daug_pdg_3)==321)){
            mc_index_Dsp_nr_kp.push_back(daug_index_1);
            mc_index_Dsp_nr_km.push_back(daug_index_3);
            mc_index_Dsp_nr_pi.push_back(daug_index_2);
          }else{
            mc_index_Dsp_nr_kp.push_back(daug_index_2);
            mc_index_Dsp_nr_km.push_back(daug_index_3);
            mc_index_Dsp_nr_pi.push_back(daug_index_1);
          }*/
        }
        //hasDs = true;

        // Dsp kinematics
        TLorentzVector mc_mom_vec;
        mc_mom_vec.SetXYZM(mcMomPx[imc], mcMomPy[imc], mcMomPz[imc], mcPartMass[imc]);
        //printf("Mass  = %f \n",mcPartMass[imc]);

        double mcRap = mc_mom_vec.Rapidity();
        double mcPt = mc_mom_vec.Pt();
        hMCDspPtRap->Fill(mcRap, mcPt); 

        // decay daughter kinematics
        if (is_phipi_decay){
          
          for(int ip = 0; ip<3; ip++){
            int mc_part_index;

            if(ip==0) mc_part_index = mc_index_Dsp_kp[mc_index_Dsp_kp.size()-1];
            if(ip==1) mc_part_index = mc_index_Dsp_km[mc_index_Dsp_km.size()-1];
            if(ip==2) mc_part_index = mc_index_Dsp_pi[mc_index_Dsp_pi.size()-1];
            
            //cout << ip << " " << mc_part_index << endl;

            TLorentzVector mc_part_vec;
            mc_part_vec.SetXYZM(mcMomPx[mc_part_index], mcMomPy[mc_part_index], mcMomPz[mc_part_index], mcPartMass[mc_part_index]);
            if(ip==0) hMcKpPtEta->Fill(mc_part_vec.Eta(), mc_part_vec.Pt());
            if(ip==1) hMcKmPtEta->Fill(mc_part_vec.Eta(), mc_part_vec.Pt());
            if(ip==2) hMcPiPtEta->Fill(mc_part_vec.Eta(), mc_part_vec.Pt());	      

            int rc_part_index = -1;
            if(assoc_map_mc_to_rc.find(mc_part_index) != assoc_map_mc_to_rc.end()) rc_part_index = assoc_map_mc_to_rc[mc_part_index];
            if (debug) 
              printf("Rec: ProngNo., MC Index, Reco Index = (%d, %d, %d) \n",ip,mc_part_index,rc_part_index);
           
            if (rc_part_index>=0){
              TVector3 dcaToVtx = GetDCAToPrimaryVertex(rc_part_index, vertex_rc);
              TVector3 mom(rcMomPx[rc_part_index], rcMomPy[rc_part_index], rcMomPz[rc_part_index]);
              hRcSecPartLocaToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.Pt());
              hRcSecPartLocbToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.z());

             // printf("Sec %d: (%2.4f, %2.4f, %2.4f), mcStartPoint = (%2.4f, %2.4f, %2.4f)\n", rc_part_index, pos.x(), pos.y(), pos.z(), mcPartVx[mc_part_index], mcPartVy[mc_part_index], mcPartVz[mc_part_index]);
            }
          }
        }
        // NR decay daughter kinematics
        if (is_kkpi_decay){
          for(int ip = 0; ip<3; ip++){
            int mc_part_index;
            if(ip==0) mc_part_index = mc_index_Dsp_nr_kp[mc_index_Dsp_nr_kp.size()-1];
            if(ip==1) mc_part_index = mc_index_Dsp_nr_km[mc_index_Dsp_nr_km.size()-1];
            if(ip==2) mc_part_index = mc_index_Dsp_nr_pi[mc_index_Dsp_nr_pi.size()-1];

            TLorentzVector mc_part_vec;
            mc_part_vec.SetXYZM(mcMomPx[mc_part_index], mcMomPy[mc_part_index], mcMomPz[mc_part_index], mcPartMass[mc_part_index]);
            if(ip==0) hMcKpPtEta->Fill(mc_part_vec.Eta(), mc_part_vec.Pt());
            if(ip==1) hMcKmPtEta->Fill(mc_part_vec.Eta(), mc_part_vec.Pt());
            if(ip==2) hMcPiPtEta->Fill(mc_part_vec.Eta(), mc_part_vec.Pt());	      

            int rc_part_index = -1;
            if(assoc_map_mc_to_rc.find(mc_part_index) != assoc_map_mc_to_rc.end()) rc_part_index = assoc_map_mc_to_rc[mc_part_index];
            if (debug) 
              printf("Rec: ProngNo., MC Index, Reco Index = (%d, %d, %d) \n",ip,mc_part_index,rc_part_index);
           
            if (rc_part_index>=0){
              TVector3 dcaToVtx = GetDCAToPrimaryVertex(rc_part_index, vertex_rc);
              TVector3 mom(rcMomPx[rc_part_index], rcMomPy[rc_part_index], rcMomPz[rc_part_index]);
              hRcSecPartLocaToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.Pt());
              hRcSecPartLocbToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.z());

             // printf("Sec %d: (%2.4f, %2.4f, %2.4f), mcStartPoint = (%2.4f, %2.4f, %2.4f)\n", rc_part_index, pos.x(), pos.y(), pos.z(), mcPartVx[mc_part_index], mcPartVy[mc_part_index], mcPartVz[mc_part_index]);
            }
          }
        }
      }
    }

    // Get reconstructed protons, kaons and pions
    hNRecoVtx->Fill(CTVx.GetSize());
    const int pid_mode = 1; // 0 - truth; 1 - realistic
    vector<unsigned int> kp_index;
    vector<unsigned int> km_index;      
    vector<unsigned int> pi_index;
    kp_index.clear();
    km_index.clear();      
    pi_index.clear();

    for (unsigned int rc_index=0; rc_index<rcMomPx.GetSize(); rc_index++){	 
      //cout << rc_index << " " << rcPdg[rc_index] << endl;
      if (pid_mode==0){
        int iSimPartID = -1;
        if(assoc_map_rc_to_mc.find(rc_index) != assoc_map_rc_to_mc.end()) iSimPartID = assoc_map_rc_to_mc[rc_index];
        if(iSimPartID>=0){
          if(fabs(mcPartPdg[iSimPartID]) == pdg_kp) kp_index.push_back(rc_index);	
          if(fabs(mcPartPdg[iSimPartID]) == pdg_km) km_index.push_back(rc_index);
          if(fabs(mcPartPdg[iSimPartID]) == pdg_pi) pi_index.push_back(rc_index);
        }
      }else if(pid_mode==1){
        if(fabs(rcPdg[rc_index]) == pdg_kp) kp_index.push_back(rc_index);
        if(fabs(rcPdg[rc_index]) == pdg_km) km_index.push_back(rc_index);
        if(fabs(rcPdg[rc_index]) == pdg_pi) pi_index.push_back(rc_index);	      
        //if (rcPdg[rc_index] == -321) cout << "where k- is " <<  rc_index << endl;
        //if (rcPdg[rc_index] == -321) km_index.push_back(rc_index);
        //if (rcPdg[rc_index] == 211) pi_index.push_back(rc_index);
      }
    }

    // printf("KaonP, KaonM, and Pion vector sizes: = (%d, %d, %d) \n",kp_index.size(), km_index.size(),pi_index.size()); 
    // Combinatorics proton, kaon, and pion 
    for(unsigned int i=0; i<kp_index.size(); i++){
      TVector3 dcaToVtx_kp = GetDCAToPrimaryVertex(kp_index[i], vertex_rc);
      int q_kaonp = rcCharge[kp_index[i]];

      for(unsigned int j=0; j<km_index.size(); j++){
        TVector3 dcaToVtx_km = GetDCAToPrimaryVertex(km_index[j], vertex_rc);
        int q_kaonm = rcCharge[km_index[j]];
        //cout << "k- charge is " << q_kaonm << endl;
        for(unsigned int k=0; k<pi_index.size(); k++){
          TVector3 dcaToVtx_pi = GetDCAToPrimaryVertex(pi_index[k], vertex_rc);
          int q_pion = rcCharge[pi_index[k]];

          if ((q_kaonp == +1 && q_kaonm == -1 && q_pion == +1) || (q_kaonp == +1 && q_kaonm == -1 && q_pion == -1)){
          //if ((q_kaonp == -1 && q_kaonm == +1 && q_pion == -1)){
            bool is_Dsp_kkpi = false;
            bool is_Dsp_phipi = false;

            int  mc_index_kp = -1, mc_index_km = -1, mc_index_pi = -1;
            if(assoc_map_rc_to_mc.find(kp_index[i]) != assoc_map_rc_to_mc.end()) mc_index_kp = assoc_map_rc_to_mc[kp_index[i]];
            if(assoc_map_rc_to_mc.find(km_index[j])  != assoc_map_rc_to_mc.end()) mc_index_km  = assoc_map_rc_to_mc[km_index[j]];
            if(assoc_map_rc_to_mc.find(pi_index[k])  != assoc_map_rc_to_mc.end()) mc_index_pi  = assoc_map_rc_to_mc[pi_index[k]];

            for(unsigned int idaugh=0; idaugh<mc_index_Dsp_pi.size(); idaugh++){
              if(mc_index_kp==mc_index_Dsp_kp[idaugh] && mc_index_km==mc_index_Dsp_km[idaugh] && mc_index_pi==mc_index_Dsp_pi[idaugh]){
                is_Dsp_phipi = true;
                //cout << "Reco to mc match " << is_Dsp_phipi << " " << q_kaonp << " " << q_kaonm << " " << q_pion << endl;
                //cout << "                      " << endl;
                break;
              }
            }
            for(unsigned int idaugh=0; idaugh<mc_index_Dsp_nr_pi.size(); idaugh++){
              if(mc_index_kp==mc_index_Dsp_nr_kp[idaugh] && mc_index_km==mc_index_Dsp_nr_km[idaugh] && mc_index_pi==mc_index_Dsp_nr_pi[idaugh]){
                is_Dsp_kkpi = true;
                break;
              }
            }

            float dcaDaughters[3], cosTheta, decayLength, V0DcaToVtx, cosTheta_xy, sigma_vtx, pT_kp, pT_km, pT_pi, px_kp, px_km, px_pi, py_kp, py_km, py_pi, pz_kp, pz_km, pz_pi;
            float chi2_ndf = 0.;
            double err_Par[6];  // or whatever size is appropriate

            TVector3 decayVertex;
            TLorentzVector parent = GetDCADaughters(kp_index[i], km_index[j], pi_index[k], vertex_rc, dcaDaughters, cosTheta, cosTheta_xy, decayLength, V0DcaToVtx, sigma_vtx, chi2_ndf, decayVertex, pT_kp, pT_km, pT_pi, err_Par, px_kp, px_km, px_pi, py_kp, py_km, py_pi, pz_kp, pz_km, pz_pi);
            hchi2_vtx->Fill(chi2_ndf);
            
            int decay_channel = -1;
            if (is_Dsp_kkpi){
              hEventStat->Fill(7.5);
              if (q_kaonp == 1 && q_kaonm == -1 && q_pion == 1)
                hEventStat->Fill(8.5);   // Ds⁺
              //else if (q_kaonp == -1 && q_kaonm == 1 && q_pion == -1)
              else if (q_kaonp == +1 && q_kaonm == -1 && q_pion == -1)
                hEventStat->Fill(9.5);   // Ds⁻
              decay_channel = 0;
            }else if (is_Dsp_phipi){
              hEventStat->Fill(3.5);
              if (q_kaonp == 1 && q_kaonm == -1 && q_pion == 1)
                hEventStat->Fill(4.5);   // Ds⁺
              //else if (q_kaonp == -1 && q_kaonm == 1 && q_pion == -1)
              else if (q_kaonp == +1 && q_kaonm == -1 && q_pion == -1)
                hEventStat->Fill(5.5);   // Ds⁻
              decay_channel = 1;
            }else{
              hEventStat->Fill(10.5);
            }

            if (is_Dsp_kkpi || is_Dsp_phipi){
              //hEventStat->Fill(7.5);
              hchi2_vtx_sig->Fill(chi2_ndf);
              TVector3 MCVertex_KaonP(mcPartVx[mc_index_kp], mcPartVy[mc_index_kp], mcPartVz[mc_index_kp]);
              //TVector3 MCVertex_Pion(mcPartVx[mc_index_pi], mcPartVy[mc_index_pi], mcPartVz[mc_index_pi]);
         
              // printf("Signal MC Vertex Kaon = (%f, %f, %f)\n",MCVertex_Kaon.X(), MCVertex_Kaon.Y(), MCVertex_Kaon.Z());	 	
              // printf("Signal MC Vertex Pion = (%f, %f, %f)\n",MCVertex_Pion.X(), MCVertex_Pion.Y(), MCVertex_Pion.Z());
         
              // cout<<"Signal MC Vertex Kaon (cm): "<<MCVertex_Kaon.X()*0.1<<"\t"<<MCVertex_Kaon.Y()*0.1<<"\t"<<MCVertex_Kaon.Z()*0.1<<endl;
              // cout<<"Signal MC Vertex Pion (cm): "<<MCVertex_Pion.X()*0.1<<"\t"<<MCVertex_Pion.Y()*0.1<<"\t"<<MCVertex_Pion.Z()*0.1<<endl; 
         
              hRes_SVx_Helixfit->Fill(parent.Pt(), parent.Rapidity(), (decayVertex.X()-MCVertex_KaonP.X()*0.1)*10);
              hRes_SVy_Helixfit->Fill(parent.Pt(), parent.Rapidity(), (decayVertex.Y()-MCVertex_KaonP.Y()*0.1)*10);
              hRes_SVz_Helixfit->Fill(parent.Pt(), parent.Rapidity(), (decayVertex.Z()-MCVertex_KaonP.Z()*0.1)*10);  

              hRes_SVx_Helixfit_pull->Fill(parent.Pt(), parent.Rapidity(), ((decayVertex.X()-MCVertex_KaonP.X()*0.1))/(err_Par[0]));
              hRes_SVy_Helixfit_pull->Fill(parent.Pt(), parent.Rapidity(), ((decayVertex.Y()-MCVertex_KaonP.Y()*0.1))/(err_Par[1]));
              hRes_SVz_Helixfit_pull->Fill(parent.Pt(), parent.Rapidity(), ((decayVertex.Z()-MCVertex_KaonP.Z()*0.1))/(err_Par[2])); 
        		     
              h3PairDca12[0]->Fill(parent.Pt(), parent.Rapidity(), dcaDaughters[0]);
              h3PairDca23[0]->Fill(parent.Pt(), parent.Rapidity(), dcaDaughters[1]); 
              h3PairDca13[0]->Fill(parent.Pt(), parent.Rapidity(), dcaDaughters[2]);		      
              h3PairCosTheta[0]->Fill(parent.Pt(), parent.Rapidity(), cosTheta);
              h3PairDca[0]->Fill(parent.Pt(), parent.Rapidity(), V0DcaToVtx);
              h3PairDecayLength[0]->Fill(parent.Pt(), parent.Rapidity(), decayLength);
              //printf("Signal: dca12 = %2.4f, cosTheta = %2.4f, D0dca = %2.4f, decay = %2.4f\n", dcaDaughters, cosTheta, V0DcaToVtx, decayLength);
              h3InvMass[0][0]->Fill(parent.Pt(), parent.Rapidity(), parent.M());

              // Toplogical Variables for Signal
              d0_kp_sig = dcaToVtx_kp.Mag();
              d0_km_sig = dcaToVtx_km.Mag();
              d0_pi_sig = dcaToVtx_pi.Mag();
              pT_kp_sig = pT_kp;	
              pT_km_sig = pT_km; 
              pT_pi_sig = pT_pi; 
              px_kp_sig = px_kp;	
              px_km_sig = px_km; 
              px_pi_sig = px_pi; 
              py_kp_sig = py_kp;	
              py_km_sig = py_km; 
              py_pi_sig = py_pi; 
              pz_kp_sig = pz_kp;	
              pz_km_sig = pz_km; 
              pz_pi_sig = pz_pi; 
              d0xy_kp_sig = dcaToVtx_kp.Perp();
              d0xy_km_sig = dcaToVtx_km.Perp();
              d0xy_pi_sig = dcaToVtx_pi.Perp();		      
              sum_d0xy_sig = sqrt(d0xy_kp_sig*d0xy_kp_sig+d0xy_km_sig*d0xy_km_sig+d0xy_pi_sig*d0xy_pi_sig);		      
              dca_12_sig = *min_element(dcaDaughters, dcaDaughters + 3);
              dca_Dsp_sig = V0DcaToVtx;
              decay_length_sig = decayLength;
              costheta_sig = cosTheta;
              costhetaxy_sig = cosTheta_xy;
              pt_Dsp_sig = parent.Pt();
              y_Dsp_sig = parent.Rapidity();
              mass_Dsp_sig = parent.M();
              sigma_vtx_sig = sigma_vtx;
              mult_sig = nMcPart;
              chi2_sig = chi2_ndf;
              xB_sig = xB_mc;
              Q2_sig = Q2_mc;
              decay = decay_channel; 
              tree_sig->Fill();  
            
            }else{
            
              hEventStat->Fill(10.5);
              hchi2_vtx_bkg->Fill(chi2_ndf);		      
              h3PairDca12[1]->Fill(parent.Pt(), parent.Rapidity(), dcaDaughters[0]);
              h3PairDca23[1]->Fill(parent.Pt(), parent.Rapidity(), dcaDaughters[1]); 
              h3PairDca13[1]->Fill(parent.Pt(), parent.Rapidity(), dcaDaughters[2]);			      
              h3PairCosTheta[1]->Fill(parent.Pt(), parent.Rapidity(), cosTheta);
              h3PairDca[1]->Fill(parent.Pt(), parent.Rapidity(), V0DcaToVtx);
              h3PairDecayLength[1]->Fill(parent.Pt(), parent.Rapidity(), decayLength);

              //printf("Bkg: dca12 = %2.4f, cosTheta = %2.4f, D0dca = %2.4f, decay = %2.4f\n", dcaDaughters, cosTheta, V0DcaToVtx, decayLength);
              h3InvMass[1][0]->Fill(parent.Pt(), parent.Rapidity(), parent.M());

              // Toplogical Variables for Bkg
              d0_kp_bkg = dcaToVtx_kp.Mag();
              d0_km_bkg = dcaToVtx_km.Mag();
              d0_pi_bkg = dcaToVtx_pi.Mag();	
              pT_kp_bkg =  pT_kp; 
              pT_km_bkg =  pT_km; 
              pT_pi_bkg =  pT_pi; 	      
              px_kp_bkg = px_kp;	
              px_km_bkg = px_km; 
              px_pi_bkg = px_pi; 
              py_kp_bkg = py_kp;	
              py_km_bkg = py_km; 
              py_pi_bkg = py_pi; 
              pz_kp_bkg = pz_kp;	
              pz_km_bkg = pz_km; 
              pz_pi_bkg = pz_pi; 
              d0xy_kp_bkg = dcaToVtx_kp.Perp();
              d0xy_km_bkg = dcaToVtx_km.Perp();
              d0xy_pi_bkg = dcaToVtx_pi.Perp();		      
              sum_d0xy_bkg = sqrt(d0xy_kp_bkg*d0xy_kp_bkg+d0xy_km_bkg*d0xy_km_bkg+d0xy_pi_bkg*d0xy_pi_bkg);	      
              dca_12_bkg = *min_element(dcaDaughters, dcaDaughters + 3);
              dca_Dsp_bkg = V0DcaToVtx;
              decay_length_bkg = decayLength;
              costheta_bkg = cosTheta;
              costhetaxy_bkg = cosTheta_xy;
              pt_Dsp_bkg = parent.Pt();
              y_Dsp_bkg = parent.Rapidity();
              mass_Dsp_bkg = parent.M();
              sigma_vtx_bkg = sigma_vtx;
              mult_bkg = nMcPart;
              chi2_bkg = chi2_ndf;
              xB_bkg = xB_mc;          
              Q2_bkg = Q2_mc;
            
              tree_bkg->Fill();
            } 
          }
        }
      }
    }

    nevents++;
  }


  tree_sig->Write();
  tree_bkg->Write();
  tree_Dsp->Write();
  hEventStat->Write();
  /*hMcMult->Write();
  hMcVtxX->Write();
  hMcVtxY->Write();
  hMcVtxZ->Write();
  hRecVtxX->Write();
  hRecVtxY->Write();
  hRecVtxZ->Write();
  hPullVtxX->Write();
  hPullVtxY->Write();
  hPullVtxZ->Write();
  hchi2_vtx->Write();
  hchi2_vtx_sig->Write();		      
  hchi2_vtx_bkg->Write();
  hDspDecayVxVy->Write();
  hDspDecayVrVz->Write();
  hRes_SVx_Helixfit->Write();
  hRes_SVy_Helixfit->Write();
  hRes_SVz_Helixfit->Write();
  hRes_SVx_Helixfit_pull->Write();
  hRes_SVy_Helixfit_pull->Write();
  hRes_SVz_Helixfit_pull->Write();

  hMCDspPtRap->Write();
  hMcKpPtEta->Write();
  hMcKpPtEtaReco->Write();
  hMcPiPtEta->Write();
  hMcPiPtEtaReco->Write();
  hMcKmPtEta->Write();
  hMcKmPtEtaReco->Write();

  hNRecoVtx->Write();

  for(int ip=0; ip<3; ip++){
    hRcSecPartLocaToRCVtx[ip]->Write();
    hRcSecPartLocbToRCVtx[ip]->Write();
    hRcPrimPartLocaToRCVtx[ip]->Write();
    hRcPrimPartLocbToRCVtx[ip]->Write();
  }

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      hPrimTrkDcaToRCVtx[i][j]->Write();
    }
  }

  for(int i=0; i<2; i++)
  {
    h3PairDca12[i]->Write();
    h3PairDca23[i]->Write();
    h3PairDca13[i]->Write();            
    h3PairCosTheta[i]->Write();
    h3PairDca[i]->Write();
    h3PairDecayLength[i]->Write();
  } 

  for(int i=0; i<2; i++){
    for(int j=0; j<2; j++){
      h3InvMass[i][j]->Write();
    }
  }
  */
  outfile->Close();

}

//======================================
TVector3 GetDCAToPrimaryVertex(const int index, TVector3 vtx)
{
  //printf("check %d: (%2.4f, %2.4f, %2.4f) =? (%2.4f, %2.4f, %2.4f)\n", index, mom.x(), mom.y(), mom.z(), rcMomPx2->At(index), rcMomPy2->At(index), rcMomPz2->At(index));

  TVector3 pos(rcTrkLoca2->At(index) * sin(rcTrkPhi2->At(index)) * -1 * millimeter, rcTrkLoca2->At(index) * cos(rcTrkPhi2->At(index)) * millimeter, rcTrkLocb2->At(index) * millimeter);
  TVector3 mom(rcMomPx2->At(index), rcMomPy2->At(index), rcMomPz2->At(index));

  StPhysicalHelix pHelix(mom, pos, bField * tesla, rcCharge2->At(index));

  TVector3 vtx_tmp;
  vtx_tmp.SetXYZ(vtx.x()*millimeter, vtx.y()*millimeter, vtx.z()*millimeter);
  
  pHelix.moveOrigin(pHelix.pathLength(vtx_tmp));
  TVector3 dcaToVtx = pHelix.origin() - vtx_tmp;

  dcaToVtx.SetXYZ(dcaToVtx.x()/millimeter, dcaToVtx.y()/millimeter, dcaToVtx.z()/millimeter);
  
  return dcaToVtx;
}

//======================================
// Local to global conversion 
// x = - l0 Sin (phi), y = l0 Cos (phi), z = l1 
TLorentzVector GetDCADaughters(const int index1, const int index2, const int index3, TVector3 vtx,
  float *dcaDaughters, float &cosTheta, float &cosTheta_xy, float &decayLength, float &V0DcaToVtx, float &sigma_vtx, float & chi2_ndf, TVector3 &decayVertex, float &pT_kp, float &pT_km, float &pT_pi, double *parFitErr, float &px_kp, float &px_km, float &px_pi, float &py_kp, float &py_km, float &py_pi, float &pz_kp, float &pz_km, float &pz_pi)
{
  // -- get helix
  TVector3 pos1(rcTrkLoca2->At(index1) * sin(rcTrkPhi2->At(index1)) * -1 * millimeter, rcTrkLoca2->At(index1) * cos(rcTrkPhi2->At(index1)) * millimeter, rcTrkLocb2->At(index1) * millimeter);
  TVector3 pos2(rcTrkLoca2->At(index2) * sin(rcTrkPhi2->At(index2)) * -1 * millimeter, rcTrkLoca2->At(index2) * cos(rcTrkPhi2->At(index2)) * millimeter, rcTrkLocb2->At(index2) * millimeter);
  TVector3 pos3(rcTrkLoca2->At(index3) * sin(rcTrkPhi2->At(index3)) * -1 * millimeter, rcTrkLoca2->At(index3) * cos(rcTrkPhi2->At(index3)) * millimeter, rcTrkLocb2->At(index3) * millimeter);
  
  TVector3 mom1(rcMomPx2->At(index1), rcMomPy2->At(index1), rcMomPz2->At(index1));
  TVector3 mom2(rcMomPx2->At(index2), rcMomPy2->At(index2), rcMomPz2->At(index2));
  TVector3 mom3(rcMomPx2->At(index3), rcMomPy2->At(index3), rcMomPz2->At(index3));

  float charge1 = rcCharge2->At(index1);
  float charge2 = rcCharge2->At(index2);
  float charge3 = rcCharge2->At(index3);  
  
  StPhysicalHelix p1Helix(mom1, pos1, bField * tesla, charge1);
  StPhysicalHelix p2Helix(mom2, pos2, bField * tesla, charge2);
  StPhysicalHelix p3Helix(mom3, pos3, bField * tesla, charge3);

  TVector3 vtx_tmp;
  vtx_tmp.SetXYZ(vtx.x()*millimeter, vtx.y()*millimeter, vtx.z()*millimeter);
  
  // Decay vertex with distance minimization
  double s1, s2, s3;
  
  getDecayVertex_Chi2fit(index1,index2,index3,s1,s2,s3,decayVertex,chi2_ndf, parFitErr);
  
  // Points of closest approach among three helices
  TVector3 const p1 = p1Helix.at(s1);
  TVector3 const p2 = p2Helix.at(s2);
  TVector3 const p3 = p3Helix.at(s3);
  // printf("Chi2/ndf = %f \n",chi2_ndf);
  // Three dcaDaughters; dcaDaughters_12, dcaDaughters_23, dcaDaughters_13
  dcaDaughters[0] = (p1 - p2).Mag()/millimeter;
  dcaDaughters[1] = (p2 - p3).Mag()/millimeter;
  dcaDaughters[2] = (p3 - p1).Mag()/millimeter; 

  // -- calculate Momentum of each daughters at the DCA point
  TVector3 const p1MomAtDca = p1Helix.momentumAt(s1,  bField * tesla);  
  TVector3 const p2MomAtDca = p2Helix.momentumAt(s2,  bField * tesla); 
  TVector3 const p3MomAtDca = p3Helix.momentumAt(s3,  bField * tesla);
  pT_kp = p1MomAtDca.Pt();   pT_km = p2MomAtDca.Pt();   pT_pi = p3MomAtDca.Pt();
  px_kp = p1MomAtDca.X();   px_km = p2MomAtDca.X();   px_pi = p3MomAtDca.X();
  py_kp = p1MomAtDca.Y();   py_km = p2MomAtDca.Y();   py_pi = p3MomAtDca.Y();
  pz_kp = p1MomAtDca.Z();   pz_km = p2MomAtDca.Z();   pz_pi = p3MomAtDca.Z();
  
  TLorentzVector p1FourMom(p1MomAtDca, sqrt(p1MomAtDca.Mag2()+gKaon1Mass*gKaon1Mass));
  TLorentzVector p2FourMom(p2MomAtDca, sqrt(p2MomAtDca.Mag2()+gKaon2Mass*gKaon2Mass));
  TLorentzVector p3FourMom(p3MomAtDca, sqrt(p3MomAtDca.Mag2()+gPionMass*gPionMass));  
  
  TLorentzVector parent = p1FourMom + p2FourMom + p3FourMom;

  // Confirm by evaluating the distance with vertices
  sigma_vtx = sqrt((p1-decayVertex).Mag2()+(p2-decayVertex).Mag2()+(p3-decayVertex).Mag2())/millimeter;

  // -- calculate pointing angle and decay length with respect to primary vertex
  //    if decay vertex is a tertiary vertex
  //    -> only rough estimate -> needs to be updated after secondary vertex is found
  TVector3 vtxToV0 = decayVertex - vtx_tmp;
  TVector3 vtxToV0_xy(vtxToV0.x(), vtxToV0.y(), 0.);
  TVector3 parent_xy(parent.Vect().x(),parent.Vect().y(),0.);
  float pointingAngle = vtxToV0.Angle(parent.Vect());
  float pointingAngle_xy = vtxToV0_xy.Angle(parent_xy);
  cosTheta = std::cos(pointingAngle);
  cosTheta_xy = std::cos(pointingAngle_xy);  
  decayLength = vtxToV0.Mag()/millimeter;

  // -- calculate V0 DCA to primary vertex
  V0DcaToVtx = decayLength * std::sin(pointingAngle);

  //TVector3 dcaToVtx = GetDCAToPrimaryVertex(parent.Vect(), decayVertex, 0, vtx);
  //V0DcaToVtx = dcaToVtx.Mag();

  return parent;
}

void getDecayVertex_Chi2fit(const int index1, const int index2, const int index3,  double &s1, double &s2, double &s3, TVector3 &vertex, float &chi2, double *parFitErr)
{
    TVector3 pos1(rcTrkLoca2->At(index1) * sin(rcTrkPhi2->At(index1)) * -1 * millimeter,
                  rcTrkLoca2->At(index1) * cos(rcTrkPhi2->At(index1)) * millimeter,
                  rcTrkLocb2->At(index1) * millimeter);

    TVector3 mom1(rcMomPx2->At(index1), rcMomPy2->At(index1), rcMomPz2->At(index1));

    TVector3 pos2(rcTrkLoca2->At(index2) * sin(rcTrkPhi2->At(index2)) * -1 * millimeter,
                  rcTrkLoca2->At(index2) * cos(rcTrkPhi2->At(index2)) * millimeter,
                  rcTrkLocb2->At(index2) * millimeter);

    TVector3 mom2(rcMomPx2->At(index2), rcMomPy2->At(index2), rcMomPz2->At(index2));
    
    TVector3 pos3(rcTrkLoca2->At(index3) * sin(rcTrkPhi2->At(index3)) * -1 * millimeter,
                  rcTrkLoca2->At(index3) * cos(rcTrkPhi2->At(index3)) * millimeter,
                  rcTrkLocb2->At(index3) * millimeter);

    TVector3 mom3(rcMomPx2->At(index3), rcMomPy2->At(index3), rcMomPz2->At(index3));


    float charge1 = rcCharge2->At(index1);
    float charge2 = rcCharge2->At(index2);
    float charge3 = rcCharge2->At(index3);


    StPhysicalHelix helix1(mom1, pos1, bField * tesla, charge1);
    StPhysicalHelix helix2(mom2, pos2, bField * tesla, charge2);
    StPhysicalHelix helix3(mom3, pos3, bField * tesla, charge3); 
    
    pair<double, double> const ss12 = helix1.pathLengths(helix2);
    pair<double, double> const ss13 = helix1.pathLengths(helix3);
    TVector3 const p1_init = helix1.at(ss12.first);
    TVector3 const p2_init = helix2.at(ss12.second); 
    TVector3 const p3_init = helix3.at(ss13.second);
    std::array<float, 21>& fcov1 = rcTrkCov->At(index1);
    std::array<float, 21>& fcov2 = rcTrkCov->At(index2); 
    std::array<float, 21>& fcov3 = rcTrkCov->At(index3);  

      // Perform Minimization
    const Int_t nPar = 6;
    Chi2Minimization d2Function(helix1,helix2,helix3,fcov1,fcov2,fcov3);
    ROOT::Math::Functor fcn(d2Function,nPar); // 5 parameters
    ROOT::Fit::Fitter fitter;

    double x0 =1./3*(p1_init.X() + p2_init.X()+ p3_init.X()); double y0 =1./3*(p1_init.Y() + p2_init.Y() + p3_init.Y()); double z0 =1./3*(p1_init.Z() + p2_init.Z() + p3_init.Z());  
    double pStart[nPar] = {x0,y0,z0,ss12.first,ss12.second,ss13.second};
    fitter.SetFCN(fcn, pStart,nPar,1);
    fitter.Config().ParSettings(0).SetName("x0");
    fitter.Config().ParSettings(0).SetStepSize(0.01);
   // fitter.Config().ParSettings(0).SetLimits(-1., 1.);    
    // No limits for x, y, z

    fitter.Config().ParSettings(1).SetName("y0");
    fitter.Config().ParSettings(1).SetStepSize(0.01);
   // fitter.Config().ParSettings(1).SetLimits(-1., 1.);

    fitter.Config().ParSettings(2).SetName("z0");
    fitter.Config().ParSettings(2).SetStepSize(0.01);
   // fitter.Config().ParSettings(2).SetLimits(-10., 10.);    

    fitter.Config().ParSettings(3).SetName("s1");
    fitter.Config().ParSettings(3).SetValue(0.0);
    fitter.Config().ParSettings(3).SetStepSize(0.01);
   // fitter.Config().ParSettings(3).SetLimits(-1., 1.);

    fitter.Config().ParSettings(4).SetName("s2");
    fitter.Config().ParSettings(4).SetValue(0.0);
    fitter.Config().ParSettings(4).SetStepSize(0.01);
   // fitter.Config().ParSettings(4).SetLimits(-1., 1.);  

    fitter.Config().ParSettings(5).SetName("s3");
    fitter.Config().ParSettings(5).SetValue(0.0);
    fitter.Config().ParSettings(5).SetStepSize(0.01);
   // fitter.Config().ParSettings(5).SetLimits(-1., 1.);          
    
    fitter.Config().MinimizerOptions().SetMaxIterations(10000);	
    // do the fit 

    Bool_t ok = fitter.FitFCN();
    //if (!ok) Error("Fitting","Fitting failed");
    const ROOT::Fit::FitResult & result = fitter.Result();
    //chi2 = fitter.Result().Chi2()/3.0;
    chi2 = fitter.Result().MinFcnValue()/3.0;  // Minimum value of your function
    
    //int status = fitter.Result().Status();
    //if (status>0 ) {printf("Fit Failed!!!!\n");}
    // if (status>0 || chi2_ndf>10. ) return;
    // cout <<"\033[1;31m Fit Result Chi2:\033[0m"<<chi2<<endl;
    //result.Print(std::cout);
   
    // Get the covariance matrix
    //  TMatrixDSym covMatrix(5);
    // result.GetCovarianceMatrix(covMatrix); // Matrix for the parameter errors
    // covMatrix.Print();
   
    const double *parFit = result.GetParams();
    const double *FitErr = result.GetErrors();
    
    for (int i = 0; i < nPar; ++i) parFitErr[i] = FitErr[i];
   
    vertex.SetXYZ(parFit[0], parFit[1], parFit[2]);
    s1 = parFit[3]; s2 = parFit[4]; s3 = parFit[5];

}

