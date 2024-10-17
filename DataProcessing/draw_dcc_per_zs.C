void draw_dcc_per_zs()
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    // Float_t text_size = 28.;
    // gStyle->SetTextSize(text_size);
    // gStyle->SetLegendTextSize(text_size);
    // gStyle->SetLabelSize(text_size, "XYZ");
    // gStyle->SetTitleSize(text_size, "XYZ");

    // TARGET and SAFE 
    TFile *fin_LTH4_HTH8 = new TFile("jobs/DoSR_LTH4p0_HTH8p0_NoZS_MIEB8p0_MIEE8p0/srvalid_HI2023_R375790.root");
    TH1F *h_eb_lth4_hth8_zs8p0 = (TH1F *) fin_LTH4_HTH8->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_eb_lth4_hth8_zs8p0");

    TFile *fin_LTH8_HTH10 = new TFile("jobs/DoSR_LTH8p0_HTH10p0_NoZS_MIEB8p0_MIEE8p0/srvalid_HI2023_R375790.root");
    TH1F *h_eb_lth8_hth10_zs8p0 = (TH1F *) fin_LTH8_HTH10->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_eb_lth8_hth10_zs8p0");
    TH1F *h_ee_lth8_hth10_zs8p0 = (TH1F *) fin_LTH8_HTH10->Get("ecalSelectiveReadoutValidation/hDccVol5")->Clone("h_ee_lth8_hth10_zs8p0");

    TFile *fin_LTH10_HTH14 = new TFile("jobs/DoSR_LTH10p0_HTH14p0_NoZS_MIEB8p0_MIEE8p0/srvalid_HI2023_R375790.root");
    TH1F *h_ee_lth10_hth14_zs8p0 = (TH1F *) fin_LTH10_HTH14->Get("ecalSelectiveReadoutValidation/hDccVol5")->Clone("h_ee_lth10_hth14_zs8p0");

    // ZS threshold study
    TFile *fin_LTH4_HTH8_ZS8p5 = new TFile("jobs/DoSR_LTH4p0_HTH8p0_NoZS_MIEB8p5_MIEE8p5/srvalid_HI2023_R375790.root");
    TH1F *h_eb_lth4_hth8_zs8p5 = (TH1F *) fin_LTH4_HTH8_ZS8p5->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_eb_lth4_hth8_zs8p5");

    TFile *fin_LTH8_HTH10_ZS8p5 = new TFile("jobs/DoSR_LTH8p0_HTH10p0_NoZS_MIEB8p5_MIEE8p5/srvalid_HI2023_R375790.root");
    TH1F *h_eb_lth8_hth10_zs8p5 = (TH1F *) fin_LTH8_HTH10_ZS8p5->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_eb_lth8_hth10_zs8p5");
    TH1F *h_ee_lth8_hth10_zs8p5 = (TH1F *) fin_LTH8_HTH10_ZS8p5->Get("ecalSelectiveReadoutValidation/hDccVol5")->Clone("h_ee_lth8_hth10_zs8p5");

    TFile *fin_LTH10_HTH14_ZS8p5 = new TFile("jobs/DoSR_LTH10p0_HTH14p0_NoZS_MIEB8p5_MIEE8p5/srvalid_HI2023_R375790.root");
    TH1F *h_ee_lth10_hth14_zs8p5 = (TH1F *) fin_LTH10_HTH14_ZS8p5->Get("ecalSelectiveReadoutValidation/hDccVol5")->Clone("h_ee_lth10_hth14_zs8p5");

    TFile *fin_LTH4_HTH8_ZS9p0 = new TFile("jobs/DoSR_LTH4p0_HTH8p0_NoZS_MIEB9p0_MIEE9p0/srvalid_HI2023_R375790.root");
    TH1F *h_eb_lth4_hth8_zs9p0 = (TH1F *) fin_LTH4_HTH8_ZS9p0->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_eb_lth4_hth8_zs9p0");

    TFile *fin_LTH8_HTH10_ZS9p0 = new TFile("jobs/DoSR_LTH8p0_HTH10p0_NoZS_MIEB9p0_MIEE9p0/srvalid_HI2023_R375790.root");
    TH1F *h_eb_lth8_hth10_zs9p0 = (TH1F *) fin_LTH8_HTH10_ZS9p0->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_eb_lth8_hth10_zs9p0");
    TH1F *h_ee_lth8_hth10_zs9p0 = (TH1F *) fin_LTH8_HTH10_ZS9p0->Get("ecalSelectiveReadoutValidation/hDccVol5")->Clone("h_ee_lth8_hth10_zs9p0");

    TFile *fin_LTH10_HTH14_ZS9p0 = new TFile("jobs/DoSR_LTH10p0_HTH14p0_NoZS_MIEB9p0_MIEE9p0/srvalid_HI2023_R375790.root");
    TH1F *h_ee_lth10_hth14_zs9p0 = (TH1F *) fin_LTH10_HTH14_ZS9p0->Get("ecalSelectiveReadoutValidation/hDccVol5")->Clone("h_ee_lth10_hth14_zs9p0");

    // Normalize
    for (auto h : {
        h_eb_lth4_hth8_zs8p0, h_eb_lth8_hth10_zs8p0,
        h_ee_lth8_hth10_zs8p0, h_ee_lth10_hth14_zs8p0,
        h_eb_lth4_hth8_zs8p5, h_eb_lth8_hth10_zs8p5,
        h_ee_lth8_hth10_zs8p5, h_ee_lth10_hth14_zs8p5,
        h_eb_lth4_hth8_zs9p0, h_eb_lth8_hth10_zs9p0,
        // h_ee_lth8_hth10_zs9p0, h_ee_lth10_hth14_zs9p0,
        }) {
        h->Scale(1/h->Integral());
    }

    // Format

    // target eb : blue circles
    h_eb_lth4_hth8_zs8p0->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_eb_lth4_hth8_zs8p0->GetXaxis()->SetTitle("Event size (kB)");
    h_eb_lth4_hth8_zs8p0->SetMarkerColor(kBlue);
    h_eb_lth4_hth8_zs8p0->SetLineColor(kBlue);
    h_eb_lth4_hth8_zs8p0->SetMarkerStyle(kFullCircle);
    h_eb_lth4_hth8_zs8p0->GetYaxis()->SetRangeUser(0.00001, 1.1);

    h_eb_lth4_hth8_zs8p5->SetMarkerColor(kBlue);
    h_eb_lth4_hth8_zs8p5->SetLineColor(kBlue);
    h_eb_lth4_hth8_zs8p5->SetMarkerStyle(kOpenCircle);

    h_eb_lth4_hth8_zs9p0->SetMarkerColor(kBlue);
    h_eb_lth4_hth8_zs9p0->SetLineColor(kBlue);
    h_eb_lth4_hth8_zs9p0->SetMarkerStyle(kOpenDiamond);

    // safe eb : red squares
    h_eb_lth8_hth10_zs8p0->SetMarkerColor(kRed);
    h_eb_lth8_hth10_zs8p0->SetLineColor(kRed);
    h_eb_lth8_hth10_zs8p0->SetMarkerStyle(kFullSquare);

    h_eb_lth8_hth10_zs8p5->SetMarkerColor(kRed);
    h_eb_lth8_hth10_zs8p5->SetLineColor(kRed);
    h_eb_lth8_hth10_zs8p5->SetMarkerStyle(kOpenSquare);

    h_eb_lth8_hth10_zs9p0->SetMarkerColor(kRed);
    h_eb_lth8_hth10_zs9p0->SetLineColor(kRed);
    h_eb_lth8_hth10_zs9p0->SetMarkerStyle(kOpenStar);

    // h_ee_lth8_hth10_zs8->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    // h_ee_lth8_hth10_zs8->GetXaxis()->SetTitle("Event size (kB)");
    // h_ee_lth8_hth10_zs8->SetMarkerColor(kBlue);
    // h_ee_lth8_hth10_zs8->SetLineColor(kBlue);
    // h_ee_lth8_hth10_zs8->SetMarkerStyle(kOpenCircle);

    // h_ee_lth10_hth14_zs8->SetMarkerColor(kRed);
    // h_ee_lth10_hth14_zs8->SetLineColor(kRed);
    // h_ee_lth10_hth14_zs8->SetMarkerStyle(kOpenSquare);

    // h_ee_lth8_hth10_zs8_8p5->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    // h_ee_lth8_hth10_zs8_8p5->GetXaxis()->SetTitle("Event size (kB)");
    // h_ee_lth8_hth10_zs8_8p5->SetMarkerColor(kBlue);
    // h_ee_lth8_hth10_zs8_8p5->SetLineColor(kBlue);
    // h_ee_lth8_hth10_zs8_8p5->SetMarkerStyle(kOpenCircle);

    // h_ee_lth10_hth14_zs8p5->SetMarkerColor(kRed);
    // h_ee_lth10_hth14_zs8p5->SetLineColor(kRed);
    // h_ee_lth10_hth14_zs8p5->SetMarkerStyle(kOpenSquare);

    // Draw 
    TLegend *leg_barrel = new TLegend(0.4, 0.58, 0.85, 0.95);
    leg_barrel->SetBorderSize(0);
    // leg_barrel->SetHeader("Barrel");
    leg_barrel->AddEntry(h_eb_lth4_hth8_zs8p0, Form("LTH=4, HTH=8, ZS=8; #mu=%.2f", h_eb_lth4_hth8_zs8p0->GetMean()));
    leg_barrel->AddEntry(h_eb_lth4_hth8_zs8p5, Form("LTH=4, HTH=8, ZS=8.5; #mu=%.2f", h_eb_lth4_hth8_zs8p5->GetMean()));
    leg_barrel->AddEntry(h_eb_lth4_hth8_zs9p0, Form("LTH=4, HTH=8, ZS=9; #mu=%.2f", h_eb_lth4_hth8_zs9p0->GetMean()));
    leg_barrel->AddEntry(h_eb_lth8_hth10_zs8p0, Form("LTH=8, HTH=10, ZS=8; #mu=%.2f", h_eb_lth8_hth10_zs8p0->GetMean()));
    leg_barrel->AddEntry(h_eb_lth8_hth10_zs8p5, Form("LTH=8, HTH=10, ZS=8.5; #mu=%.2f", h_eb_lth8_hth10_zs8p5->GetMean()));
    leg_barrel->AddEntry(h_eb_lth8_hth10_zs9p0, Form("LTH=8, HTH=10, ZS=9; #mu=%.2f", h_eb_lth8_hth10_zs9p0->GetMean()));

    TCanvas *c_barrel = new TCanvas("c_barrel", "", 800, 600);
    c_barrel->SetLeftMargin(0.15);
    c_barrel->SetTopMargin(0.01);
    c_barrel->SetLogy();
    h_eb_lth4_hth8_zs8p0->Draw();
    h_eb_lth4_hth8_zs8p5->Draw("same");
    h_eb_lth4_hth8_zs9p0->Draw("same");
    h_eb_lth8_hth10_zs8p0->Draw("same");
    h_eb_lth8_hth10_zs8p5->Draw("same");
    h_eb_lth8_hth10_zs9p0->Draw("same");
    leg_barrel->Draw();

    // TLegend *leg_endcap = new TLegend(0.5, 0.6, 0.8, 0.85);
    // leg_endcap->SetBorderSize(0);
    // leg_endcap->SetHeader("Endcap");
    // leg_endcap->AddEntry(h_ee_lth8_hth10_zs8, Form("Target; #mu=%.2f", h_ee_lth8_hth10_zs8->GetMean()));
    // leg_endcap->AddEntry(h_ee_lth10_hth14_zs8, Form("Safe; #mu=%.2f", h_ee_lth10_hth14_zs8->GetMean()));

    // TCanvas *c_endcap = new TCanvas("c_endcap", "", 800, 600);
    // c_endcap->SetLeftMargin(0.15);
    // c_endcap->SetTopMargin(0.01);
    // c_endcap->SetLogy();
    // h_ee_lth8_hth10_zs8->Draw();
    // h_ee_lth10_hth14_zs8->Draw("same");
    // leg_endcap->Draw();
}