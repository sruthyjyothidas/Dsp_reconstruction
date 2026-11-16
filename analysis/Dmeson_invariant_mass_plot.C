//This macro plots signal and background distribution and their optimizations
//Manually needs to set optimized bin
Double_t pol2Fit(Double_t *x, Double_t *par)
{
    Double_t a = par[0];
    Double_t b = par[1];
    Double_t c = par[2];

    return a + b*x[0] + c*x[0]*x[0];
}
//For pretty plot
Double_t background(Double_t *x, Double_t *par) {
   return (par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0]);
}
Double_t doubleGaussian(Double_t *x, Double_t *par) {
  return par[0] * (exp(-0.5* TMath::Power((x[0]-par[1])/par[2],2)) + (par[3]*exp(-0.5* TMath::Power((x[0]-par[1])/(par[4]*par[2]),2))));
}
Double_t fitFunction(Double_t *x, Double_t *par) {
  return background(x,par) + doubleGaussian(x,&par[4]);
}
void Dmeson_invariant_mass_plot(){
    
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
     
    float ml; float mu;
    const char* dau; 
    if (D0){
        ml = 1.73;
        mu = 2;
        dau = "m_{K^{+} #pi^{-}}";
    }
    if (Ds){
        ml = 1.8;
        mu = 2.12;
        dau = "m_{K^{+} K^{-} #pi^{+}}";
    }
    const float PtValues[6] = {0.5, 1.0, 1.5, 3.0, 5.0, 20.0}; 

    int Nbins = 100;
    float binw = (mu-ml)*1000./Nbins;
    
    TFile* f = new TFile(Form("../root/outputFiles/invariant_mass_%s_%s.root",sys,spec),"READ");
    TH1D* hMassSig = (TH1D*)f->Get("hMassSig");
    TH1D* hMassBkg = (TH1D*)f->Get("hMassBkg");
    const int nPtbins = 5;
    TH1D* hMassSigPt[nPtbins];
    TH1D* hMassBkgPt[nPtbins];
    for (int i = 0; i < nPtbins; i++) {
        hMassSigPt[i] = (TH1D*)f->Get(Form("hMassSigPt%d",i));
        hMassBkgPt[i] = (TH1D*)f->Get(Form("hMassBkgPt%d",i));
    }

    //Correct uncertainties
    for (int ig = 0; ig < hMassSig->GetNbinsX(); ig++){
        hMassSig->SetBinError(ig+1, sqrt(hMassSig->GetBinContent(ig+1)));
        for (int i = 0; i < nPtbins; i++) {
            hMassSigPt[i]->SetBinError(ig+1, sqrt(hMassSigPt[i]->GetBinContent(ig+1)));
        }
    }

    TF1* fit = new TF1(Form("fit"),pol2Fit, ml, mu, 3);
    fit->SetParameters(1.0, 1.0, 0.0);
    fit->FixParameter(2, 0.0);
    hMassBkg->Fit(fit,"RQN");
    double nSamples = hMassBkg->Integral();
    int nRandom = (int)nSamples;
    //cout << nRandom << endl;
    TH1D* hMassBkgSam = (TH1D*)hMassBkg->Clone("hMassBkgSam");
    hMassBkgSam->Reset();
    hMassBkgSam->Sumw2();

    for (int j = 0; j < nRandom; ++j) {
        double mass = fit->GetRandom();
        hMassBkgSam->Fill(mass);
    }

    //Fits for Pt-differential cases
    TF1* fitPt[nPtbins];
    TH1D* hMassBkgSamPt[nPtbins];
    for (int i = 0; i < nPtbins; i++) {
        fitPt[i] = new TF1(Form("fit%d",i),pol2Fit, ml, mu, 3);
        fitPt[i]->SetParameters(1.0, 1.0, 0.0);
        fitPt[i]->FixParameter(2, 0.0);
        hMassBkgPt[i]->Fit(fitPt[i],"RQN");
        double nSamples = hMassBkgPt[i]->Integral();
        int nRandom = (int)nSamples;
        hMassBkgSamPt[i] = (TH1D*)hMassBkgPt[i]->Clone(Form("hMassBkgSamPt%d",i));
        hMassBkgSamPt[i]->Reset();
        hMassBkgSamPt[i]->Sumw2();

        for (int j = 0; j < nRandom; ++j) {
            double mass = fitPt[i]->GetRandom();
             hMassBkgSamPt[i]->Fill(mass);
        }
    }

    TH1D* hMassCombined = (TH1D*)hMassSig->Clone("hMassCombined");
    hMassCombined->Add(hMassBkgSam);
    
    TH1D* hMassCombinedPt[nPtbins];
    for (int i = 0; i < nPtbins; i++) {
        hMassCombinedPt[i] = (TH1D*)hMassSigPt[i]->Clone(Form("hMassCombinedPt%d",i));
        hMassCombinedPt[i]->Add(hMassBkgSamPt[i]);
    }

    auto leg_sys = new TLegend(0.35,0.95,0.41,0.96);
    if (eAu) leg_sys = new TLegend(0.32,0.95,0.38,0.96);
    leg_sys->SetBorderSize(0);
    leg_sys->SetTextSize(0.045);
    leg_sys->SetTextFont(42);
    leg_sys->SetMargin(0.0);
    if (ep) leg_sys->AddEntry((TObject*)0,"PYTHIA ep @ 9x275 GeV, Q^{2} > 1","");
    if (eAu) leg_sys->AddEntry((TObject*)0,"BeAGLE eAu @ 9x100 GeV, Q^{2} > 1","");
    
    auto leg_epic = new TLegend(0.14,0.85,0.15,0.88);
    leg_epic->SetBorderSize(0);
    leg_epic->SetTextSize(0.05);
    leg_epic->SetTextFont(42);
    leg_epic->SetMargin(0.5);
    leg_epic->AddEntry((TObject*)0,"#bf{ePIC} #it{Simulation}","");

    TCanvas* c1 = new TCanvas("c1", "New Canvas", 50, 50, 600, 600);
    TPad* pad1 = new TPad("pad1","",0.,0.,1.0,1.0);

    pad1->Draw();
    pad1->cd();
    pad1->SetTopMargin(0.06);
    pad1->SetLeftMargin(0.11);
    
    TH1F* hframe1 = pad1->DrawFrame(ml,0,mu,1.3*hMassSig->GetMaximum());
    hframe1->GetXaxis()->SetTitle(Form("%s (GeV)",dau));
    hframe1->GetYaxis()->SetTitle(Form("Candidates / %.1f MeV", (binw)));
    hframe1->GetXaxis()->CenterTitle(true);
    hframe1->GetYaxis()->CenterTitle(true);
    hframe1->GetXaxis()->SetTitleOffset(1.1);
    hframe1->GetYaxis()->SetTitleOffset(1.0);
    hframe1->GetXaxis()->SetTitleSize(0.052);
    hframe1->GetXaxis()->SetLabelSize(0.045);
    hframe1->GetYaxis()->SetTitleSize(0.05);
    hframe1->GetYaxis()->SetLabelSize(0.045);
    hframe1->GetXaxis()->SetNdivisions(509,kTRUE);
    hframe1->GetYaxis()->SetNdivisions(506,kTRUE);
    hframe1->GetXaxis()->SetTickLength(0.025);
    hframe1->GetYaxis()->SetMaxDigits(3);
    hframe1->Draw("AXIS");
    gStyle->SetErrorX(0);

    hMassSig->Draw("pesame");
    hMassSig->SetMarkerSize(1.0);
    
    auto leg1 = new TLegend(0.15,0.7,0.25,0.8);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.04);
    leg1->SetTextFont(42);
    leg1->SetMargin(0.1);
    if (D0) leg1->AddEntry("",Form("#it{N}_{D_{0}^{#pm}}^{Signal}: %d",int(hMassSig->Integral())),"");
    if (Ds) leg1->AddEntry("",Form("#it{N}_{D_{s}^{#pm}}^{Signal}: %d",int(hMassSig->Integral())),"");
    leg1->Draw();

    leg_sys->Draw();
    leg_epic->Draw();

    c1->SaveAs(Form("../plots/DMeson/Invariant_mass_%s_%s_signal.pdf",sys, spec));
    c1->Close();
    delete c1;
    
    TCanvas* c2 = new TCanvas("c2", "New Canvas", 50, 50, 600, 600);
    TPad* pad2 = new TPad("pad2","",0.,0.,1.0,1.0);

    pad2->Draw();
    pad2->cd();
    pad2->SetTopMargin(0.06);
    pad2->SetLeftMargin(0.11);
    
    TH1F* hframe2 = pad2->DrawFrame(ml,0,mu,1.6*hMassBkg->GetMaximum());
    hframe2->GetXaxis()->SetTitle(Form("%s (GeV)",dau));
    hframe2->GetYaxis()->SetTitle(Form("Candidates / %.1f MeV", (binw)));
    hframe2->GetXaxis()->CenterTitle(true);
    hframe2->GetYaxis()->CenterTitle(true);
    hframe2->GetXaxis()->SetTitleOffset(1.1);
    hframe2->GetYaxis()->SetTitleOffset(1.0);
    hframe2->GetXaxis()->SetTitleSize(0.052);
    hframe2->GetXaxis()->SetLabelSize(0.045);
    hframe2->GetYaxis()->SetTitleSize(0.05);
    hframe2->GetYaxis()->SetLabelSize(0.045);
    hframe2->GetXaxis()->SetNdivisions(509,kTRUE);
    hframe2->GetYaxis()->SetNdivisions(506,kTRUE);
    hframe2->GetXaxis()->SetTickLength(0.025);
    hframe2->GetYaxis()->SetMaxDigits(3);
    hframe2->Draw("AXIS");
    gStyle->SetErrorX(0);

    hMassBkg->Draw("pesame");
    hMassBkg->SetMarkerSize(1.0);
    fit->Draw("lsame");

    auto leg2 = new TLegend(0.15,0.7,0.25,0.8);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.04);
    leg2->SetTextFont(42);
    leg2->SetMargin(0.1);
    if (D0) leg2->AddEntry("",Form("#it{N}_{D_{0}^{#pm}}^{Bkg}: %d",int(hMassBkg->Integral())),"");
    if (Ds) leg2->AddEntry("",Form("#it{N}_{D_{s}^{#pm}}^{Bkg}: %d",int(hMassBkg->Integral())),"");
    leg2->Draw();
    
    auto leg2a = new TLegend(0.6,0.7,0.7,0.8);
    leg2a->SetBorderSize(0);
    leg2a->SetTextSize(0.04);
    leg2a->SetTextFont(42);
    leg2a->SetMargin(0.6);
    leg2a->AddEntry(hMassBkg,"Background", "pe");
    leg2a->AddEntry(fit, "Pol. 1 fit", "l");
    leg2a->Draw();

    leg_sys->Draw();
    leg_epic->Draw();

    c2->SaveAs(Form("../plots/DMeson/Invariant_mass_%s_%s_bkg.pdf",sys, spec));
    c2->Close();
    delete c2;
    
    TCanvas* c3 = new TCanvas("c3", "New Canvas", 50, 50, 600, 600);
    TPad* pad3 = new TPad("pad3","",0.,0.,1.0,1.0);

    pad3->Draw();
    pad3->cd();
    pad3->SetTopMargin(0.06);
    pad3->SetLeftMargin(0.11);
    
    TH1F* hframe3 = pad3->DrawFrame(ml,0,mu,1.6*hMassBkg->GetMaximum());
    hframe3->GetXaxis()->SetTitle(Form("%s (GeV)",dau));
    hframe3->GetYaxis()->SetTitle(Form("Candidates / %.1f MeV", (binw)));
    hframe3->GetXaxis()->CenterTitle(true);
    hframe3->GetYaxis()->CenterTitle(true);
    hframe3->GetXaxis()->SetTitleOffset(1.1);
    hframe3->GetYaxis()->SetTitleOffset(1.0);
    hframe3->GetXaxis()->SetTitleSize(0.052);
    hframe3->GetXaxis()->SetLabelSize(0.045);
    hframe3->GetYaxis()->SetTitleSize(0.05);
    hframe3->GetYaxis()->SetLabelSize(0.045);
    hframe3->GetXaxis()->SetNdivisions(509,kTRUE);
    hframe3->GetYaxis()->SetNdivisions(506,kTRUE);
    hframe3->GetXaxis()->SetTickLength(0.025);
    hframe3->GetYaxis()->SetMaxDigits(3);
    hframe3->Draw("AXIS");
    gStyle->SetErrorX(0);

    hMassBkgSam->Draw("pesame");
    hMassBkgSam->SetMarkerSize(1.0);
    
    auto leg3 = new TLegend(0.15,0.7,0.25,0.8);
    leg3->SetBorderSize(0);
    leg3->SetTextSize(0.04);
    leg3->SetTextFont(42);
    leg3->SetMargin(0.1);
    if (D0) leg3->AddEntry("",Form("#it{N}_{D_{0}^{#pm}}^{Bkg}: %d",int(hMassBkg->Integral())),"");
    if (Ds) leg3->AddEntry("",Form("#it{N}_{D_{s}^{#pm}}^{Bkg}: %d",int(hMassBkg->Integral())),"");
    leg3->Draw();

    leg_sys->Draw();
    leg_epic->Draw();

    c3->SaveAs(Form("../plots/DMeson/Invariant_mass_%s_%s_bkgsampled.pdf",sys, spec));
    c3->Close();
    delete c3;

    //Signal optimization
    double mass;
    if (D0) mass = 1.86484;
    if (Ds) mass = 1.96834;
    const int nbins = 47;
    int begin; 
    int end;
    if (D0){
        begin = 1;
        end = 100;
    }
    if (Ds){
        begin = 6;
        end = 100;
    }
    float x[nbins] = {};
    float y[nbins] = {};
    TGraphAsymmErrors* hSigSigEff = new TGraphAsymmErrors(nbins,x,y);
    TGraphAsymmErrors* hSigEff = new TGraphAsymmErrors(nbins,x,y);
    TGraphAsymmErrors* hBkgEff = new TGraphAsymmErrors(nbins,x,y);
    
    float opt_eff = 0.95;
    int best_begin = 1;
    int best_end = 100; 
    for (int ig = 0; ig < nbins; ig++){
        double eff = hMassSig->Integral(begin+ig,end-ig)/hMassSig->Integral();
        double sig = hMassSig->Integral(begin+ig,end-ig);
        double bkg = hMassBkgSam->Integral(begin+ig,end-ig);
        double purity = sig/bkg;
        double sigsig = sig/sqrt(sig+bkg);
        cout << begin+ig << " " << hMassSig->GetBinCenter(begin+ig) << " " << end-ig << " " << hMassSig->GetBinCenter(end-ig) << " eff: " << eff << " sig: " << sig << " bkg: " << bkg << " sigsig: " << sigsig << " unc: " << 1.0/sigsig << endl;
        if (Ds && abs(eff-opt_eff) <= 3.56e-3){ 
            cout << "My best bins " << begin+ig << " " << end-ig << " " << abs(eff-opt_eff) << endl;
            best_begin = begin+ig;
            best_end = end-ig;
        }
        if (D0 && abs(eff-opt_eff) <= 2e-3){
            cout << "My best bins " << begin+ig << " " << end-ig << " " << abs(eff-opt_eff) << endl;
            best_begin = begin+ig;
            best_end = end-ig;
        }
        hSigSigEff->SetPoint(ig, eff, sigsig);
        hSigEff->SetPoint(ig, eff, sig);
        hBkgEff->SetPoint(ig, eff, bkg);
    }
    
    TCanvas* c4 = new TCanvas("c4", "New Canvas", 50, 50, 600, 600);
    TPad* pad4 = new TPad("pad4","",0.,0.,1.0,1.0);

    pad4->Draw();
    pad4->cd();
    pad4->SetTopMargin(0.06);
    pad4->SetLeftMargin(0.11);
    
    TH1F* hframe4 = pad4->DrawFrame(ml,0,mu,1.4*hMassCombined->GetMaximum());
    hframe4->GetXaxis()->SetTitle(Form("%s (GeV)",dau));
    hframe4->GetYaxis()->SetTitle(Form("Candidates / %.1f MeV", (binw)));
    hframe4->GetXaxis()->CenterTitle(true);
    hframe4->GetYaxis()->CenterTitle(true);
    hframe4->GetXaxis()->SetTitleOffset(1.1);
    hframe4->GetYaxis()->SetTitleOffset(1.0);
    hframe4->GetXaxis()->SetTitleSize(0.052);
    hframe4->GetXaxis()->SetLabelSize(0.045);
    hframe4->GetYaxis()->SetTitleSize(0.05);
    hframe4->GetYaxis()->SetLabelSize(0.045);
    hframe4->GetXaxis()->SetNdivisions(509,kTRUE);
    hframe4->GetYaxis()->SetNdivisions(506,kTRUE);
    hframe4->GetXaxis()->SetTickLength(0.025);
    hframe4->GetYaxis()->SetMaxDigits(3);
    hframe4->Draw("AXIS");
    gStyle->SetErrorX(0);

    hMassCombined->Draw("pesame");
    hMassCombined->SetMarkerSize(0.5);
    hMassCombined->SetLineWidth(2);
    
    double signal_significance = hMassSig->Integral(best_begin,best_end)/sqrt(hMassSig->Integral(best_begin,best_end)+hMassBkgSam->Integral(best_begin,best_end));

    auto leg4 = new TLegend(0.15,0.7,0.25,0.8);
    leg4->SetBorderSize(0);
    leg4->SetTextSize(0.04);
    leg4->SetTextFont(42);
    leg4->SetMargin(0.1);
    leg4->AddEntry("",Form("S/#sqrt{S+B} = %.1f",signal_significance),"");
    leg4->Draw();

    leg_sys->Draw();
    leg_epic->Draw();

    c4->SaveAs(Form("../plots/DMeson/Invariant_mass_%s_%s_sigbkgcombined.pdf",sys, spec));
    c4->Close();
    delete c4;
    
    //Pretty Ds plots
    if (Ds){
        TCanvas* c4 = new TCanvas("c4", "New Canvas", 50, 50, 600, 600);
        TPad* pad4 = new TPad("pad4","",0.,0.,1.0,1.0);

        pad4->Draw();
        pad4->cd();
        pad4->SetTopMargin(0.06);
        pad4->SetLeftMargin(0.11);

        TH1F* hframe4 = pad4->DrawFrame(ml,0,mu,1.4*hMassCombined->GetMaximum());
        hframe4->GetXaxis()->SetTitle(Form("%s (GeV)",dau));
        hframe4->GetYaxis()->SetTitle(Form("Candidates / %.1f MeV", hMassCombined->GetBinWidth(1)*1000.0));
        hframe4->GetXaxis()->CenterTitle(true);
        hframe4->GetYaxis()->CenterTitle(true);
        hframe4->GetXaxis()->SetTitleOffset(1.1);
        hframe4->GetYaxis()->SetTitleOffset(1.0);
        hframe4->GetXaxis()->SetTitleSize(0.052);
        hframe4->GetXaxis()->SetLabelSize(0.045);
        hframe4->GetYaxis()->SetTitleSize(0.05);
        hframe4->GetYaxis()->SetLabelSize(0.045);
        hframe4->GetXaxis()->SetNdivisions(509,kTRUE);
        hframe4->GetYaxis()->SetNdivisions(506,kTRUE);
        hframe4->GetXaxis()->SetTickLength(0.025);
        hframe4->GetYaxis()->SetMaxDigits(3);
        hframe4->Draw("AXIS");
        gStyle->SetErrorX(0);

        hMassCombined->Draw("pesame");
        hMassCombined->SetMarkerSize(0.5);
        hMassCombined->SetLineWidth(2);
        
        //Fit
        const int npar = 9;
        Double_t par[npar];
        TF1 *fitFcn = new TF1("fitFcn",fitFunction,ml,mu,npar);
        // initial values
        fitFcn->SetParameters(0,0,0,0,  // bkg
                          1,        // A1
                          mass,    // mean
                          0.0015,    // sigma1
                          0.3,      // A2/A1
                          2.0);     // sigma2/sigma1
        
        if (eAu) fitFcn->SetParLimits(4, 0, 2e3);
        fitFcn->SetParLimits(5, mass-0.1, mass+0.1);
        //fitFcn->SetParLimits(6, 0.001, 0.02); // sigma1
        //fitFcn->SetParLimits(7, 0.0, 1.0);    // A2/A1
        //fitFcn->SetParLimits(8, 1.0, 3.0);    // sigma2/sigma1

        fitFcn->SetNpx(500);
        fitFcn->SetLineWidth(2);
        fitFcn->SetLineColor(kRed);

        TFitResultPtr r;
        for (int i = 0; i < 20; i++){
            if (i!=19){
                r = hMassCombined->Fit("fitFcn","SRQN","");
            }else{
                r = hMassCombined->Fit("fitFcn","SRN","");
            }
        }

        TF1 *backFcn = new TF1("backFcn",background,ml,mu,4);
        backFcn->SetLineColor(kGreen+2);
        backFcn->SetLineWidth(2);
        backFcn->SetLineStyle(3);

        TF1 *signalFcn = new TF1("signalFcn",doubleGaussian,ml,mu,5);
        signalFcn->SetLineColor(kAzure+1);
        signalFcn->SetLineWidth(2);
        signalFcn->SetNpx(500);

        fitFcn->GetParameters(par);
        backFcn->SetParameters(par);
        signalFcn->SetParameters(&par[4]);

        backFcn->Draw("same");
        signalFcn->Draw("same");
        fitFcn->Draw("same");
        //hMassSig->Draw("pesame");
        
        float sig = signalFcn->Integral(ml,mu)/hMassSig->GetBinWidth(1);
        float bkg = backFcn->Integral(ml,mu)/hMassSig->GetBinWidth(1);
        
        TLegend *legend = new TLegend(0.65,0.54,0.72,0.77);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.04);
        legend->SetTextFont(42);
        legend->SetMargin(0.8);
        legend->AddEntry(hMassCombined,"Data ","pe");
        legend->AddEntry(fitFcn,"Fit","l");
        legend->AddEntry(signalFcn,"Signal","l");
        legend->AddEntry(backFcn,"Background","l");
        legend->Draw();

        auto leg1 = new TLegend(0.15,0.76,0.2,0.82);
        leg1->SetBorderSize(0);
        leg1->SetTextSize(0.04);
        leg1->SetTextFont(42);
        leg1->SetMargin(0.1);
        leg1->AddEntry("",Form("D_{s}^{+} #rightarrow #phi + #pi^{+} #rightarrow K^{+} + K^{-} + #pi^{+} & c.c."), "");
        leg1->Draw();

        //double signal_significance = hMassSig->Integral(best_begin,best_end)/sqrt(hMassSig->Integral(best_begin,best_end)+hMassBkgSam->Integral(best_begin,best_end));

        //auto leg4 = new TLegend(0.15,0.7,0.25,0.8);
        //leg4->SetBorderSize(0);
        //leg4->SetTextSize(0.04);
        //leg4->SetTextFont(42);
        //leg4->SetMargin(0.1);
        //leg4->AddEntry("",Form("S/#sqrt{S+B} = %.1f, %.1f",sig/sqrt(sig+bkg), hMassSig->Integral()/sqrt(hMassSig->Integral()+hMassBkgSam->Integral())),"");
        //leg4->AddEntry("",Form("S = %.1f, %.1f", signalFcn->Integral(ml,mu)/hMassSig->GetBinWidth(1), hMassSig->Integral( hMassSig->FindBin(ml),  hMassSig->FindBin(mu))), "");
        //leg4->AddEntry("",Form("B = %.1f, %.1f", backFcn->Integral(ml,mu)/hMassSig->GetBinWidth(1), hMassBkgSam->Integral()), "");
        //leg4->Draw();

        leg_sys->Draw();
        leg_epic->Draw();

        c4->SaveAs(Form("../plots/DMeson/Invariant_mass_%s_%s_sigbkgcombined_pretty.pdf",sys, spec));
        c4->Close();
        delete c4;
    }

    TCanvas* c5 = new TCanvas("c5", "New Canvas", 50, 50, 600, 600);
    TPad* pad5 = new TPad("pad5","",0.,0.,1.0,1.0);

    pad5->Draw();
    pad5->cd();
    pad5->SetTopMargin(0.06);
    pad5->SetLeftMargin(0.11);

    TH1F* hframe5 = pad5->DrawFrame(0.5,0,1.0,hBkgEff->GetPointY(0));
    hframe5->GetXaxis()->SetTitle(Form("Signal efficiency"));
    hframe5->GetYaxis()->SetTitle(Form("Counts"));
    hframe5->GetXaxis()->CenterTitle(true);
    hframe5->GetYaxis()->CenterTitle(true);
    hframe5->GetXaxis()->SetTitleOffset(1.1);
    hframe5->GetYaxis()->SetTitleOffset(1.0);
    hframe5->GetXaxis()->SetTitleSize(0.052);
    hframe5->GetXaxis()->SetLabelSize(0.045);
    hframe5->GetYaxis()->SetTitleSize(0.05);
    hframe5->GetYaxis()->SetLabelSize(0.045);
    hframe5->GetXaxis()->SetNdivisions(509,kTRUE);
    hframe5->GetYaxis()->SetNdivisions(506,kTRUE);
    hframe5->GetXaxis()->SetTickLength(0.025);
    hframe5->GetYaxis()->SetMaxDigits(3);
    hframe5->Draw("AXIS");
    gStyle->SetErrorX(0);
    
    hSigEff->Draw("pesame");
    hBkgEff->Draw("pesame");
    
    hSigEff->SetMarkerColor(2);
    hSigEff->SetLineColor(2);

    leg_sys->Draw();
    leg_epic->Draw();
    
    auto leg5 = new TLegend(0.15,0.7,0.25,0.8);
    leg5->SetBorderSize(0);
    leg5->SetTextSize(0.04);
    leg5->SetTextFont(42);
    leg5->SetMargin(0.5);
    leg5->AddEntry(hSigEff,Form("Signal"),"pe");
    leg5->AddEntry(hBkgEff,Form("Background"),"pe");
    leg5->Draw();

    c5->SaveAs(Form("../plots/DMeson/Invariant_mass_%s_%s_signaloptimization.pdf",sys, spec));
    c5->Close();
    delete c5;

    TCanvas* c6 = new TCanvas("c6", "New Canvas", 50, 50, 600, 600);
    TPad* pad6 = new TPad("pad6","",0.,0.,1.0,1.0);

    pad6->Draw();
    pad6->cd();
    pad6->SetTopMargin(0.06);
    pad6->SetLeftMargin(0.13);

    TH1F* hframe6 = pad6->DrawFrame(0.5,0,1.0,D0?350:100);
    hframe6->GetXaxis()->SetTitle(Form("Signal efficiency"));
    hframe6->GetYaxis()->SetTitle(Form("S/#sqrt{S+B}"));
    hframe6->GetXaxis()->CenterTitle(true);
    hframe6->GetYaxis()->CenterTitle(true);
    hframe6->GetXaxis()->SetTitleOffset(1.1);
    hframe6->GetYaxis()->SetTitleOffset(1.3);
    hframe6->GetXaxis()->SetTitleSize(0.052);
    hframe6->GetXaxis()->SetLabelSize(0.045);
    hframe6->GetYaxis()->SetTitleSize(0.05);
    hframe6->GetYaxis()->SetLabelSize(0.045);
    hframe6->GetXaxis()->SetNdivisions(509,kTRUE);
    hframe6->GetYaxis()->SetNdivisions(506,kTRUE);
    hframe6->GetXaxis()->SetTickLength(0.025);
    hframe6->GetYaxis()->SetMaxDigits(3);
    hframe6->Draw("AXIS");
    gStyle->SetErrorX(0);
    
    hSigSigEff->Draw("pesame");

    leg_sys->Draw();
    leg_epic->Draw();

    TLine *line = new TLine(0.95, 0, 0.95, 400);
    line->Draw();
    
    c6->SaveAs(Form("../plots/DMeson/Invariant_mass_%s_%s_signaloptimization_signalsignificance.pdf",sys, spec));
    c6->Close();
    delete c6;
    
    //Plot pT-dependent mass histograms 
    for (int i = 0; i < nPtbins; i++) {
        TCanvas* c4 = new TCanvas("c4", "New Canvas", 50, 50, 600, 600);
        TPad* pad4 = new TPad("pad4","",0.,0.,1.0,1.0);

        pad4->Draw();
        pad4->cd();
        pad4->SetTopMargin(0.06);
        pad4->SetLeftMargin(0.11);
        
        TH1F* hframe4 = pad4->DrawFrame(ml,0,mu,1.4*hMassCombinedPt[i]->GetMaximum());
        hframe4->GetXaxis()->SetTitle(Form("%s (GeV)",dau));
        hframe4->GetYaxis()->SetTitle(Form("Candidates / %.1f MeV", (binw)));
        hframe4->GetXaxis()->CenterTitle(true);
        hframe4->GetYaxis()->CenterTitle(true);
        hframe4->GetXaxis()->SetTitleOffset(1.1);
        hframe4->GetYaxis()->SetTitleOffset(1.0);
        hframe4->GetXaxis()->SetTitleSize(0.052);
        hframe4->GetXaxis()->SetLabelSize(0.045);
        hframe4->GetYaxis()->SetTitleSize(0.05);
        hframe4->GetYaxis()->SetLabelSize(0.045);
        hframe4->GetXaxis()->SetNdivisions(509,kTRUE);
        hframe4->GetYaxis()->SetNdivisions(506,kTRUE);
        hframe4->GetXaxis()->SetTickLength(0.025);
        hframe4->GetYaxis()->SetMaxDigits(3);
        hframe4->Draw("AXIS");
        gStyle->SetErrorX(0);

        hMassCombinedPt[i]->Draw("pesame");
        hMassCombinedPt[i]->SetMarkerSize(0.5);
        hMassCombinedPt[i]->SetLineWidth(2);
        
        double signal_significance = hMassSigPt[i]->Integral(best_begin,best_end)/sqrt(hMassSigPt[i]->Integral(best_begin,best_end)+hMassBkgSamPt[i]->Integral(best_begin,best_end));

        auto leg4 = new TLegend(0.15,0.7,0.25,0.85);
        leg4->SetBorderSize(0);
        leg4->SetTextSize(0.04);
        leg4->SetTextFont(42);
        leg4->SetMargin(0.1);
        leg4->AddEntry("",Form("%.1f < p_{T} < %.1f",PtValues[i],PtValues[i+1]),"");
        leg4->AddEntry("",Form("S/#sqrt{S+B} = %.1f",signal_significance),"");
        leg4->Draw();
    
        leg_sys->Draw();
        leg_epic->Draw();

        c4->SaveAs(Form("../plots/DMeson/Invariant_mass_%s_%s_sigbkgcombined_Pt%d.pdf",sys, spec, i));
        c4->Close();
        delete c4;
    }
    
    //Root file for statistical projection with optimization at 95%
    float xS[nPtbins] = {0.707, 1.225, 2.121, 3.873, 10.0};
    float yS[nPtbins] = {};
    double xlow[nPtbins]  = {0.207, 0.225, 0.621, 0.873, 5.0};
    double xhigh[nPtbins] = {0.293, 0.275, 0.879, 1.127, 10.0};
    TGraphAsymmErrors* hStat = new TGraphAsymmErrors(nPtbins,xS,yS);
    for (int i = 0; i < nPtbins; i++) {
        double sig = hMassSigPt[i]->Integral(best_begin,best_end);
        double bkg = hMassBkgSamPt[i]->Integral(best_begin,best_end);
        double statUnc = sqrt(sig+bkg)/sig;
        hStat->SetPoint(i, xS[i], 1);
        hStat->SetPointError(i, xlow[i], xhigh[i], statUnc, statUnc);
    }

    TFile *fout = new TFile(Form("../root/outputFiles/Statistical_projection_%s_%s.root", sys, spec), "RECREATE");
    fout->cd();
    hStat->Write("hStatUnc");
    fout->Close();
}
