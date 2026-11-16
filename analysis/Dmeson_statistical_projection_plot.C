void Dmeson_statistical_projection_plot(){

    TFile* fepD0 = new TFile("../root/outputFiles/Statistical_projection_ep_D0.root", "READ");
    TFile* fepDs = new TFile("../root/outputFiles/Statistical_projection_ep_Ds.root", "READ");
    TFile* feAuD0 = new TFile("../root/outputFiles/Statistical_projection_eAu_D0.root", "READ");
    TFile* feAuDs = new TFile("../root/outputFiles/Statistical_projection_eAu_Ds.root", "READ");
    
    TGraphAsymmErrors* gStatUncepD0 = (TGraphAsymmErrors*)fepD0->Get("hStatUnc");
    TGraphAsymmErrors* gStatUncepDs = (TGraphAsymmErrors*)fepDs->Get("hStatUnc");
    TGraphAsymmErrors* gStatUnceAuD0 = (TGraphAsymmErrors*)feAuD0->Get("hStatUnc");
    TGraphAsymmErrors* gStatUnceAuDs = (TGraphAsymmErrors*)feAuDs->Get("hStatUnc");
    
    const int nPtbins = 5;
    float xS[nPtbins] = {0.707, 1.225, 2.121, 3.873, 10.0};
    float yS[nPtbins] = {};
    double xlow[nPtbins]  = {0.207, 0.225, 0.621, 0.873, 5.0};
    double xhigh[nPtbins] = {0.293, 0.275, 0.879, 1.127, 10.0};

    TGraphAsymmErrors* gStatUncepDsD0Ratio = new TGraphAsymmErrors(nPtbins,xS,yS);
    TGraphAsymmErrors* gStatUnceAuDsD0Ratio = new TGraphAsymmErrors(nPtbins,xS,yS);
   
    for (int i = 0; i < nPtbins; i++) {
        double xepDs, yepDs, xepD0, yepD0;
        double xeAuDs, yeAuDs, xeAuD0, yeAuD0;

        gStatUncepDs->GetPoint(i, xepDs, yepDs);
        gStatUncepD0->GetPoint(i, xepD0, yepD0);
        gStatUnceAuDs->GetPoint(i, xeAuDs, yeAuDs);
        gStatUnceAuD0->GetPoint(i, xeAuD0, yeAuD0);

        double epRatio = yepDs / yepD0;
        double eAuRatio = yeAuDs / yeAuD0;
        
        // Statistical uncertainties
        double epErrLow = epRatio * sqrt(
            pow(gStatUncepDs->GetErrorYlow(i) / yepDs, 2) +
            pow(gStatUncepD0->GetErrorYhigh(i) / yepD0, 2)
        );

        double epErrHigh = epRatio * sqrt(
            pow(gStatUncepDs->GetErrorYhigh(i) / yepDs, 2) +
            pow(gStatUncepD0->GetErrorYlow(i) / yepD0, 2)
        );
        
        double eAuErrLow = eAuRatio * sqrt(
            pow(gStatUnceAuDs->GetErrorYlow(i) / yeAuDs, 2) +
            pow(gStatUnceAuD0->GetErrorYhigh(i) / yeAuD0, 2)
        );

        double eAuErrHigh = eAuRatio * sqrt(
            pow(gStatUnceAuDs->GetErrorYhigh(i) / yeAuDs, 2) +
            pow(gStatUnceAuD0->GetErrorYlow(i) / yeAuD0, 2)
        );

        gStatUncepDsD0Ratio->SetPoint(i, xepDs, epRatio);
        gStatUncepDsD0Ratio->SetPointError(i, gStatUncepDs->GetErrorXlow(i), gStatUncepDs->GetErrorXhigh(i), epErrLow, epErrHigh);
        
        gStatUnceAuDsD0Ratio->SetPoint(i, xeAuDs, eAuRatio);
        gStatUnceAuDsD0Ratio->SetPointError(i, gStatUnceAuDs->GetErrorXlow(i), gStatUnceAuDs->GetErrorXhigh(i), eAuErrLow, eAuErrHigh);
        
        cout << epErrHigh << " " << eAuErrHigh << endl;
    }

    auto leg_eAu = new TLegend(0.32,0.95,0.38,0.96);
    leg_eAu->SetBorderSize(0);
    leg_eAu->SetTextSize(0.045);
    leg_eAu->SetTextFont(42);
    leg_eAu->SetMargin(0.0);
    leg_eAu->AddEntry((TObject*)0,"BeAGLE eAu @ 9x100 GeV, Q^{2} > 1","");
    
    auto leg_ep = new TLegend(0.35,0.95,0.41,0.96);
    leg_ep->SetBorderSize(0);
    leg_ep->SetTextSize(0.045);
    leg_ep->SetTextFont(42);
    leg_ep->SetMargin(0.0);
    leg_ep->AddEntry((TObject*)0,"PYTHIA ep @ 9x275 GeV, Q^{2} > 1","");

    auto leg_epic = new TLegend(0.14,0.85,0.15,0.88);
    leg_epic->SetBorderSize(0);
    leg_epic->SetTextSize(0.05);
    leg_epic->SetTextFont(42);
    leg_epic->SetMargin(0.5);
    leg_epic->AddEntry((TObject*)0,"#bf{ePIC} #it{Simulation}","");
    
    double yval = 1.0;
    double cap = 0.01;
   
    TLine *leftCapDs[nPtbins]; 
    TLine *rightCapDs[nPtbins]; 
    TLine *leftCapD0[nPtbins]; 
    TLine *rightCapD0[nPtbins]; 
    for (int ig = 0; ig < nPtbins; ig++) {

        double xx = xS[ig];

        // vertical caps at both ends
        leftCapDs[ig] = new TLine(xx-xlow[ig], yval-0.01,
                                   xx-xlow[ig], yval+0.01);
        rightCapDs[ig] = new TLine(xx+xhigh[ig], yval-0.01,
                                    xx+xhigh[ig], yval+0.01);
        leftCapD0[ig] = new TLine(xx-xlow[ig], yval-0.003,
                                   xx-xlow[ig], yval+0.003);
        rightCapD0[ig] = new TLine(xx+xhigh[ig], yval-0.003,
                                    xx+xhigh[ig], yval+0.003);

    }

    //ep D0 uncertainty 
    TCanvas* c0 = new TCanvas("c0", "New Canvas", 50, 50, 600, 600);
    TPad* pad0 = new TPad("pad0","",0.,0.,1.0,1.0);

    pad0->Draw();
    pad0->cd();
    pad0->SetTopMargin(0.06);
    pad0->SetLeftMargin(0.11);
    gPad->SetLogx(1.0);

    TH1F* hframe0 = pad0->DrawFrame(0.5,0.85,22,1.15);
    hframe0->GetXaxis()->SetTitle(Form("p_{T} (GeV)"));
    hframe0->GetYaxis()->SetTitle(Form("D^{0} counts"));
    hframe0->GetXaxis()->CenterTitle(true);
    hframe0->GetYaxis()->CenterTitle(true);
    hframe0->GetXaxis()->SetTitleOffset(1.1);
    hframe0->GetYaxis()->SetTitleOffset(1.0);
    hframe0->GetXaxis()->SetTitleSize(0.052);
    hframe0->GetXaxis()->SetLabelSize(0.045);
    hframe0->GetYaxis()->SetTitleSize(0.05);
    hframe0->GetYaxis()->SetLabelSize(0.045);
    hframe0->GetXaxis()->SetNdivisions(509,kTRUE);
    hframe0->GetYaxis()->SetNdivisions(506,kTRUE);
    hframe0->GetXaxis()->SetTickLength(0.025);
    hframe0->GetYaxis()->SetMaxDigits(3);
    hframe0->Draw("AXIS");
    gStyle->SetErrorX(0);
    
    gStatUncepD0->Draw("psame");
    
    leg_ep->Draw();
    leg_epic->Draw();
    for (int ig = 0; ig < nPtbins; ig++) {
        leftCapD0[ig]->Draw();
        rightCapD0[ig]->Draw();
    }

    c0->SaveAs(Form("../plots/Dmeson/Statistical_projection_ep_D0.pdf"));
    c0->Close();
    delete c0;

    //ep Ds uncertainty 
    TCanvas* c1 = new TCanvas("c1", "New Canvas", 50, 50, 600, 600);
    TPad* pad1 = new TPad("pad1","",0.,0.,1.0,1.0);

    pad1->Draw();
    pad1->cd();
    pad1->SetTopMargin(0.06);
    pad1->SetLeftMargin(0.11);
    gPad->SetLogx(1.0);

    TH1F* hframe1 = pad1->DrawFrame(0.5,0.5,22,1.5);
    hframe1->GetXaxis()->SetTitle(Form("p_{T} (GeV)"));
    hframe1->GetYaxis()->SetTitle(Form("D_{s} counts"));
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
    
    gStatUncepDs->Draw("psame");
    
    leg_ep->Draw();
    leg_epic->Draw();
    for (int ig = 0; ig < nPtbins; ig++) {
        leftCapDs[ig]->Draw();
        rightCapDs[ig]->Draw();
    }

    c1->SaveAs(Form("../plots/Dmeson/Statistical_projection_ep_Ds.pdf"));
    c1->Close();
    delete c1;

    //eAu D0 uncertainty 
    TCanvas* c2 = new TCanvas("c2", "New Canvas", 50, 50, 600, 600);
    TPad* pad2 = new TPad("pad2","",0.,0.,1.0,1.0);

    pad2->Draw();
    pad2->cd();
    pad2->SetTopMargin(0.06);
    pad2->SetLeftMargin(0.11);
    gPad->SetLogx(1.0);

    TH1F* hframe2 = pad2->DrawFrame(0.5,0.85,22,1.15);
    hframe2->GetXaxis()->SetTitle(Form("p_{T} (GeV)"));
    hframe2->GetYaxis()->SetTitle(Form("D^{0} counts"));
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
    
    gStatUnceAuD0->Draw("psame");
    
    leg_eAu->Draw();
    leg_epic->Draw();
    for (int ig = 0; ig < nPtbins; ig++) {
        leftCapD0[ig]->Draw();
        rightCapD0[ig]->Draw();
    }

    c2->SaveAs(Form("../plots/Dmeson/Statistical_projection_eAu_D0.pdf"));
    c2->Close();
    delete c2;
    
    //eAu Ds uncertainty 
    TCanvas* c3 = new TCanvas("c3", "New Canvas", 50, 50, 600, 600);
    TPad* pad3 = new TPad("pad3","",0.,0.,1.0,1.0);

    pad3->Draw();
    pad3->cd();
    pad3->SetTopMargin(0.06);
    pad3->SetLeftMargin(0.11);
    gPad->SetLogx(1.0);

    TH1F* hframe3 = pad3->DrawFrame(0.5,0.5,22,1.5);
    hframe3->GetXaxis()->SetTitle(Form("p_{T} (GeV)"));
    hframe3->GetYaxis()->SetTitle(Form("D_{s} counts"));
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
    
    gStatUnceAuDs->Draw("psame");
    
    leg_eAu->Draw();
    leg_epic->Draw();
    for (int ig = 0; ig < nPtbins; ig++) {
        leftCapDs[ig]->Draw();
        rightCapDs[ig]->Draw();
    }

    c3->SaveAs(Form("../plots/Dmeson/Statistical_projection_eAu_Ds.pdf"));
    c3->Close();
    delete c3;
    
    //ep Ds/D0 uncertainty
    TCanvas* c4 = new TCanvas("c4", "New Canvas", 50, 50, 600, 600);
    TPad* pad4 = new TPad("pad4","",0.,0.,1.0,1.0);

    pad4->Draw();
    pad4->cd();
    pad4->SetTopMargin(0.06);
    pad4->SetLeftMargin(0.11);
    gPad->SetLogx(1.0);

    TH1F* hframe4 = pad4->DrawFrame(0.5,0.5,22,1.5);
    hframe4->GetXaxis()->SetTitle(Form("p_{T} (GeV)"));
    hframe4->GetYaxis()->SetTitle(Form("D_{s}/D^{0}"));
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
    
    gStatUncepDsD0Ratio->Draw("psame");

    leg_ep->Draw();
    leg_epic->Draw();
    
    for (int ig = 0; ig < nPtbins; ig++) {
        leftCapDs[ig]->Draw();
        rightCapDs[ig]->Draw();
    }

    c4->SaveAs(Form("../plots/Dmeson/Statistical_projection_ep_DsD0Ratio.pdf"));
    c4->Close();
    delete c4;
    
    //eAu Ds/D0 uncertainty
    TCanvas* c5 = new TCanvas("c5", "New Canvas", 50, 50, 600, 600);
    TPad* pad5 = new TPad("pad5","",0.,0.,1.0,1.0);

    pad5->Draw();
    pad5->cd();
    pad5->SetTopMargin(0.06);
    pad5->SetLeftMargin(0.11);
    gPad->SetLogx(1.0);

    TH1F* hframe5 = pad5->DrawFrame(0.5,0.5,22,1.5);
    hframe5->GetXaxis()->SetTitle(Form("p_{T} (GeV)"));
    hframe5->GetYaxis()->SetTitle(Form("D_{s}/D^{0}"));
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
    
    gStatUnceAuDsD0Ratio->Draw("psame");
    
    leg_eAu->Draw();
    leg_epic->Draw();

    for (int ig = 0; ig < nPtbins; ig++) {
        leftCapDs[ig]->Draw();
        rightCapDs[ig]->Draw();
    }

    c5->SaveAs(Form("../plots/Dmeson/Statistical_projection_eAu_DsD0Ratio.pdf"));
    c5->Close();
    delete c5;
    

}
