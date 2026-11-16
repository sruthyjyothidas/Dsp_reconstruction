//This code produces a root file for signal and background invariant mass distribution in ep and eAu sample for D0 and Ds mesons.
void Dmeson_invariant_mass(){
   
    bool D0 = 0;
    bool Ds = 1;

    bool ep = 0;
    bool eAu = 1;
    
    const char* spec = "";
    if (D0) spec = "D0";
    if (Ds) spec = "Ds";
    const char* sys = "";
    if (ep) sys = "ep";
    if (eAu) sys = "eAu";

    TFile* fSig = new TFile(Form("../root/inputFiles/%s/%sNtuple_%s.root",spec,spec,sys),"READ"); 
    TFile* fBkg = new TFile(Form("../root/inputFiles/%s/%sNtuple_%sDIS.root",spec,spec,sys),"READ"); 

    TTree *treeMLSig = (TTree*)fSig->Get("treeMLSig");
    TTree *treeMLBkg = (TTree*)fBkg->Get("treeMLBkg");
    
    Long64_t nSig = treeMLSig->GetEntries();
    Long64_t nBkg = treeMLBkg->GetEntries();

    std::cout << "Signal entries: " << nSig << std::endl;
    std::cout << "Background entries: " << nBkg << std::endl;
    
    Float_t mass_sig;
    Float_t pt_sig;
    Float_t costheta_sig, decay_length_sig, d0_k_sig, d0_pi_sig, ds_kp_sig, ds_km_sig, ds_pi_sig, dca_Dsp_sig;
    Int_t decay;
    Float_t px_kp_sig, px_km_sig, px_pi_sig;
    Float_t py_kp_sig, py_km_sig, py_pi_sig;
    Float_t pz_kp_sig, pz_km_sig, pz_pi_sig;

    treeMLSig->SetBranchAddress("costheta", &costheta_sig);
    treeMLSig->SetBranchAddress("decay_length", &decay_length_sig);
    if (D0){
        treeMLSig->SetBranchAddress("mass_D0", &mass_sig);
        treeMLSig->SetBranchAddress("pt_D0", &pt_sig);
        treeMLSig->SetBranchAddress("d0_k", &d0_k_sig);
        treeMLSig->SetBranchAddress("d0_pi", &d0_pi_sig);
    }
    if (Ds){
        treeMLSig->SetBranchAddress("mass_Dsp", &mass_sig);
        treeMLSig->SetBranchAddress("pt_Dsp", &pt_sig);
        treeMLSig->SetBranchAddress("decay", &decay);
        treeMLSig->SetBranchAddress("dca_Dsp", &dca_Dsp_sig);
        treeMLSig->SetBranchAddress("d0_kp", &ds_kp_sig);
        treeMLSig->SetBranchAddress("d0_km", &ds_km_sig);
        treeMLSig->SetBranchAddress("d0_pi", &ds_pi_sig);
        treeMLSig->SetBranchAddress("px_kp", &px_kp_sig);
        treeMLSig->SetBranchAddress("py_kp", &py_kp_sig);
        treeMLSig->SetBranchAddress("pz_kp", &pz_kp_sig);
        treeMLSig->SetBranchAddress("px_km", &px_km_sig);
        treeMLSig->SetBranchAddress("py_km", &py_km_sig);
        treeMLSig->SetBranchAddress("pz_km", &pz_km_sig);
        treeMLSig->SetBranchAddress("px_pi", &px_pi_sig);
        treeMLSig->SetBranchAddress("py_pi", &py_pi_sig);
        treeMLSig->SetBranchAddress("pz_pi", &pz_pi_sig);
    }

    Float_t mass_bkg;
    Float_t pt_bkg;
    Float_t costheta_bkg, decay_length_bkg, d0_k_bkg, d0_pi_bkg, ds_kp_bkg, ds_km_bkg, ds_pi_bkg, dca_Dsp_bkg;
    Float_t Q2;
    Float_t px_kp_bkg, px_km_bkg, px_pi_bkg;
    Float_t py_kp_bkg, py_km_bkg, py_pi_bkg;
    Float_t pz_kp_bkg, pz_km_bkg, pz_pi_bkg;
    
    treeMLBkg->SetBranchAddress("Q2", &Q2);
    treeMLBkg->SetBranchAddress("costheta", &costheta_bkg);
    treeMLBkg->SetBranchAddress("decay_length", &decay_length_bkg);
    
    if (D0){
        treeMLBkg->SetBranchAddress("mass_D0", &mass_bkg);
        treeMLBkg->SetBranchAddress("pt_D0", &pt_bkg);
        treeMLBkg->SetBranchAddress("d0_k", &d0_k_bkg);
        treeMLBkg->SetBranchAddress("d0_pi", &d0_pi_bkg);
    }
    if (Ds){
        treeMLBkg->SetBranchAddress("mass_Dsp", &mass_bkg);
        treeMLBkg->SetBranchAddress("pt_Dsp", &pt_bkg);
        treeMLBkg->SetBranchAddress("dca_Dsp", &dca_Dsp_bkg);
        treeMLBkg->SetBranchAddress("d0_kp", &ds_kp_bkg);
        treeMLBkg->SetBranchAddress("d0_km", &ds_km_bkg);
        treeMLBkg->SetBranchAddress("d0_pi", &ds_pi_bkg);
        treeMLBkg->SetBranchAddress("px_kp", &px_kp_bkg);
        treeMLBkg->SetBranchAddress("py_kp", &py_kp_bkg);
        treeMLBkg->SetBranchAddress("pz_kp", &pz_kp_bkg);
        treeMLBkg->SetBranchAddress("px_km", &px_km_bkg);
        treeMLBkg->SetBranchAddress("py_km", &py_km_bkg);
        treeMLBkg->SetBranchAddress("pz_km", &pz_km_bkg);
        treeMLBkg->SetBranchAddress("px_pi", &px_pi_bkg);
        treeMLBkg->SetBranchAddress("py_pi", &py_pi_bkg);
        treeMLBkg->SetBranchAddress("pz_pi", &pz_pi_bkg);
    }

    float ml; float mu;
    if (D0){
        ml = 1.73;
        mu = 2; 
    }
    if (Ds){
        ml = 1.8;
        mu = 2.12; 
    }


    TH1D* hMassSig = new TH1D("hMassSig","",100, ml, mu);
    TH1D* hMassBkg = new TH1D("hMassBkg","",100, ml, mu);
    TH1D* hEventQ2 = (TH1D*)fBkg->Get("hEventQ2");
    TH1D* hEventStat = (TH1D*)fBkg->Get("hEventStat");

    const int nbins = 5; 
    TH1D* hMassSigPt[nbins];
    TH1D* hMassBkgPt[nbins];
    for (int i = 0; i < nbins; i++) {
        hMassSigPt[i] = new TH1D(Form("hMassSigPt%d", i),"",100, ml, mu);
        hMassBkgPt[i] = new TH1D(Form("hMassBkgPt%d", i),"",100, ml, mu);
    }

    float wt_sig;
    if (Ds && ep) wt_sig = 0.137;
    if (D0 && ep) wt_sig = 1.6666667;
    if (Ds && eAu) wt_sig = 0.094;
    if (D0 && eAu) wt_sig = 0.94339623;
    
    float wt_bkg_ep[4] = {299.87, 25.85, 1.159, 0.094};
    float wt_bkg_eAu[3] = {209.19, 21.82, 0.609};

    for (Long64_t i = 0; i < treeMLSig->GetEntries(); i++) {
        treeMLSig->GetEntry(i);
        
        if (Ds && decay!=1) continue; //1-Phi+pi (2.32), 0-KKpi(5.5)
        
        //SV cuts 
        if (costheta_sig < 0.9) continue;
        if (decay_length_sig < 0.05) continue;

        //displaced tracks
        if (D0){ 
            if (d0_k_sig < 0.005) continue;
            if (d0_pi_sig < 0.005) continue;
        }
        if (Ds){
            if (ds_kp_sig < 0.005) continue;
            if (ds_km_sig < 0.005) continue;
            if (ds_pi_sig < 0.005) continue;
        }
        //Phi mass cut for Ds
        if (Ds){
            TLorentzVector kp, km, pi;
            kp.SetXYZM(px_kp_sig, py_kp_sig, pz_kp_sig, 0.493677);
            km.SetXYZM(px_km_sig, py_km_sig, pz_km_sig, 0.493677);
            pi.SetXYZM(px_pi_sig, py_pi_sig, pz_pi_sig, 0.139);

            TLorentzVector phi = kp + km;
            double phiMass = phi.M();
            if (fabs(phiMass - 1.0195) > 0.02) continue;
        }

        hMassSig->Fill(mass_sig, wt_sig);
        if (pt_sig > 0.5 && pt_sig < 1.0) hMassSigPt[0]->Fill(mass_sig, wt_sig);
        if (pt_sig > 1.0 && pt_sig < 1.5) hMassSigPt[1]->Fill(mass_sig, wt_sig);
        if (pt_sig > 1.5 && pt_sig < 3.0) hMassSigPt[2]->Fill(mass_sig, wt_sig);
        if (pt_sig > 3.0 && pt_sig < 5.0) hMassSigPt[3]->Fill(mass_sig, wt_sig);
        if (pt_sig > 5.0 && pt_sig < 20.0) hMassSigPt[4]->Fill(mass_sig, wt_sig);
    }

    for (Long64_t i = 0; i < treeMLBkg->GetEntries(); i++) {
        treeMLBkg->GetEntry(i);
        float wt_q2;
        if (ep){
            if (Q2 > 1 && Q2 <= 10) wt_q2 = wt_bkg_ep[0];
            if (Q2 > 10 && Q2 <= 100) wt_q2 = wt_bkg_ep[1];
            if (Q2 > 100 && Q2 <= 1000) wt_q2 = wt_bkg_ep[2];
            if (Q2 > 1000) wt_q2 = wt_bkg_ep[3];
        }
        if (eAu){
            if (Q2 > 1 && Q2 <= 10) wt_q2 = wt_bkg_eAu[0];
            if (Q2 > 10 && Q2 <= 100) wt_q2 = wt_bkg_eAu[1];
            if (Q2 > 100) wt_q2 = wt_bkg_eAu[2];
        }
        
        //SV cuts
        if (costheta_bkg < 0.9) continue;
        if (decay_length_bkg < 0.05) continue;

        //Displaced tracks
        if (D0){
            if (d0_k_bkg < 0.005) continue;
            if (d0_pi_bkg < 0.005) continue;
        }
        if (Ds){
            if (ds_kp_bkg < 0.005) continue;
            if (ds_km_bkg < 0.005) continue;
            if (ds_pi_bkg < 0.005) continue;
        }
        
        //Phi mass cut for Ds
        if (Ds){
            TLorentzVector kp, km, pi;
            kp.SetXYZM(px_kp_bkg, py_kp_bkg, pz_kp_bkg, 0.493677);
            km.SetXYZM(px_km_bkg, py_km_bkg, pz_km_bkg, 0.493677);
            pi.SetXYZM(px_pi_bkg, py_pi_bkg, pz_pi_bkg, 0.139);

            TLorentzVector phi = kp + km;
            double phiMass = phi.M();
            if (fabs(phiMass - 1.0195) > 0.02) continue;
        }

        hMassBkg->Fill(mass_bkg, wt_q2);
        if (pt_bkg > 0.5 && pt_bkg < 1.0) hMassBkgPt[0]->Fill(mass_bkg, wt_q2);
        if (pt_bkg > 1.0 && pt_bkg < 1.5) hMassBkgPt[1]->Fill(mass_bkg, wt_q2);
        if (pt_bkg > 1.5 && pt_bkg < 3.0) hMassBkgPt[2]->Fill(mass_bkg, wt_q2);
        if (pt_bkg > 3.0 && pt_bkg < 5.0) hMassBkgPt[3]->Fill(mass_bkg, wt_q2);
        if (pt_bkg > 5.0 && pt_bkg < 20.0) hMassBkgPt[4]->Fill(mass_bkg, wt_q2);
    }
    
    TFile* fout = new TFile(Form("../root/outputFiles/invariant_mass_%s_%s.root",sys,spec),"RECREATE");
    fout->cd();
    hMassSig->Write();
    hMassBkg->Write();
    for (int i = 0; i < nbins; i++) {
        hMassSigPt[i]->Write();
        hMassBkgPt[i]->Write();
    }
    fout->Close();
    
    /*TFile* fout = new TFile(Form("../root/output_%s_srdas.root",tag),"RECREATE");
    fout->cd();   
    hEventStat->Write(); 
    hPhiMassSig->Write(); 
    hPhiMassBkg->Write(); 
    hDsMassSig->Write();
    hDsNRMass->Write();
    hDsMassBkg->Write();
    for (int i = 0; i < nbins; i++) {
        hDsMasspT[i]->Write();
    }
    fout->Close();
    */
}
