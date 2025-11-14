void draw_dcc_per_year()
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    // Float_t text_size = 28.;
    // gStyle->SetTextSize(text_size);
    // gStyle->SetLegendTextSize(text_size);
    // gStyle->SetLabelSize(text_size, "XYZ");
    // gStyle->SetTitleSize(text_size, "XYZ");

    // --- BARREL --- 

    TFile *fin_2023_barrel = new TFile("/afs/cern.ch/user/l/lkalipol/private/ecal-hin/CMSSW_13_2_4/src/DataProcessing/jobs/DoSR_LTH4p0_HTH8p0_NoZS_MIEB8p0_MIEE8p0/srvalid_HI2023_R375513.root");
    TH1F *h_2023_barrel = (TH1F *) fin_2023_barrel->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_2023_barrel");

    TFile *fin_2024_barrel = new TFile("jobs/DoSR_LTH4p0_HTH8p0_NoZS_MIEB8p0_MIEE8p0/srvalid_hi2024_r387855.root");
    TH1F *h_2024_barrel = (TH1F *) fin_2024_barrel->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_2024_barrel");

    TFile *fin_2025_barrel = new TFile("jobs/DoSR_LTH4p0_HTH8p0_NoZS_MIEB8p0_MIEE8p0/srvalid_hi2025setup_r399002.root");
    TH1F *h_2025_barrel = (TH1F *) fin_2025_barrel->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_2025_barrel");

    // TFile *fin_2024_zs9 = new TFile("jobs/DoSR_LTH4p0_HTH8p0_NoZS_MIEB9p0_MIEE8p0/srvalid_hi2024_r387855.root");
    // TH1F *h_2024_barrel_zs9 = (TH1F *) fin_2024_zs9->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_2024_barrel_zs9");

    // TFile *fin_2024_zs10 = new TFile("jobs/DoSR_LTH4p0_HTH8p0_NoZS_MIEB10p0_MIEE8p0/srvalid_hi2024_r387855.root");
    // TH1F *h_2024_barrel_zs10 = (TH1F *) fin_2024_zs10->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_2024_barrel_zs10");

    // Normalize
    for (auto h : {
        h_2023_barrel, h_2024_barrel, h_2025_barrel,
        // h_2024_barrel_zs9, h_2024_barrel_zs10
        }) {
            // h->GetXaxis()->SetRange(1, h->GetNbinsX()+1);
            // std::cout << h->GetBinContent(101) << std::endl;
            h->Scale(1/h->Integral());
    }

    // Format
    h_2023_barrel->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_2023_barrel->GetXaxis()->SetTitle("Event size (kB)");
    h_2023_barrel->SetMarkerColor(kBlack);
    h_2023_barrel->SetLineColor(kBlack);
    h_2023_barrel->SetMarkerStyle(kFullCircle);

    h_2024_barrel->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_2024_barrel->GetXaxis()->SetTitle("Event size (kB)");
    h_2024_barrel->SetMarkerColor(kBlue);
    h_2024_barrel->SetLineColor(kBlue);
    h_2024_barrel->SetMarkerStyle(kFullCircle);

    h_2025_barrel->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_2025_barrel->GetXaxis()->SetTitle("Event size (kB)");
    h_2025_barrel->SetMarkerColor(kRed);
    h_2025_barrel->SetLineColor(kRed);
    h_2025_barrel->SetMarkerStyle(kFullCircle);

    // h_2024_barrel_zs9->SetMarkerColor(kRed);
    // h_2024_barrel_zs9->SetLineColor(kRed);
    // h_2024_barrel_zs9->SetMarkerStyle(kOpenCircle);

    // h_2024_barrel_zs10->SetMarkerColor(kGreen);
    // h_2024_barrel_zs10->SetLineColor(kGreen);
    // h_2024_barrel_zs10->SetMarkerStyle(kOpenSquare);

    // Draw 
    TLegend *leg_barrel = new TLegend(0.5, 0.6, 0.8, 0.85);
    leg_barrel->SetBorderSize(0);
    leg_barrel->SetHeader("Barrel, LTH=4, HTH=8");
    leg_barrel->AddEntry(h_2023_barrel, Form("2023 (375513) ZS=8; #mu=%.2f", h_2023_barrel->GetMean()));
    leg_barrel->AddEntry(h_2024_barrel, Form("2024 (387855) ZS=8; #mu=%.2f", h_2024_barrel->GetMean()));
    leg_barrel->AddEntry(h_2025_barrel, Form("2025 (399002) ZS=8; #mu=%.2f", h_2025_barrel->GetMean()));
    // leg_barrel->AddEntry(h_2024_barrel_zs9, Form("387855 ZS=9; #mu=%.2f", h_2024_barrel_zs9->GetMean()));
    // leg_barrel->AddEntry(h_2024_barrel_zs10, Form("387855 ZS=10; #mu=%.2f", h_2024_barrel_zs10->GetMean()));

    TCanvas *c_barrel = new TCanvas("c_barrel", "", 800, 600);
    c_barrel->SetLeftMargin(0.15);
    c_barrel->SetTopMargin(0.01);
    c_barrel->SetLogy();
    h_2023_barrel->Draw();
    h_2024_barrel->Draw("same");
    h_2025_barrel->Draw("same");
    // h_2024_barrel_zs9->Draw("same");
    // h_2024_barrel_zs10->Draw("same");
    leg_barrel->Draw();

    // --- ENDCAP --- 

    TFile *fin_2023_endcap = new TFile("/afs/cern.ch/user/l/lkalipol/private/ecal-hin/CMSSW_13_2_4/src/DataProcessing/jobs/DoSR_LTH8p0_HTH10p0_NoZS_MIEB8p0_MIEE8p0/srvalid_HI2023_R375513.root");
    TH1F *h_2023_endcap = (TH1F *) fin_2023_endcap->Get("ecalSelectiveReadoutValidation/hDccVol5")->Clone("h_2023_endcap");

    TFile *fin_2024_endcap = new TFile("jobs/DoSR_LTH8p0_HTH10p0_NoZS_MIEB8p0_MIEE8p0/srvalid_hi2024_r387855.root");
    TH1F *h_2024_endcap = (TH1F *) fin_2024_endcap->Get("ecalSelectiveReadoutValidation/hDccVol5")->Clone("h_2024_endcap");

    TFile *fin_2025_endcap = new TFile("jobs/DoSR_LTH8p0_HTH10p0_NoZS_MIEB8p0_MIEE8p0/srvalid_hi2025setup_r399002.root");
    TH1F *h_2025_endcap = (TH1F *) fin_2025_endcap->Get("ecalSelectiveReadoutValidation/hDccVol5")->Clone("h_2025_endcap");

    // TFile *fin_2024_zs9 = new TFile("jobs/DoSR_LTH8p0_HTH10p0_NoZS_MIEB9p0_MIEE8p0/srvalid_hi2024_r387855.root");
    // TH1F *h_2024_endcap_zs9 = (TH1F *) fin_2024_zs9->Get("ecalSelectiveReadoutValidation/hDccVol5")->Clone("h_2024_endcap_zs9");

    // TFile *fin_2024_zs10 = new TFile("jobs/DoSR_LTH8p0_HTH10p0_NoZS_MIEB10p0_MIEE8p0/srvalid_hi2024_r387855.root");
    // TH1F *h_2024_endcap_zs10 = (TH1F *) fin_2024_zs10->Get("ecalSelectiveReadoutValidation/hDccVol5")->Clone("h_2024_endcap_zs10");

    // Normalize
    for (auto h : {
        h_2023_endcap, h_2024_endcap, h_2025_endcap,
        // h_2024_endcap_zs9, h_2024_endcap_zs10
        }) {
            // h->GetXaxis()->SetRange(1, h->GetNbinsX()+1);
            // std::cout << h->GetBinContent(101) << std::endl;
            h->Scale(1/h->Integral());
    }

    // Format
    h_2023_endcap->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_2023_endcap->GetXaxis()->SetTitle("Event size (kB)");
    h_2023_endcap->SetMarkerColor(kBlack);
    h_2023_endcap->SetLineColor(kBlack);
    h_2023_endcap->SetMarkerStyle(kFullCircle);

    h_2024_endcap->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_2024_endcap->GetXaxis()->SetTitle("Event size (kB)");
    h_2024_endcap->SetMarkerColor(kBlue);
    h_2024_endcap->SetLineColor(kBlue);
    h_2024_endcap->SetMarkerStyle(kFullCircle);

    h_2025_endcap->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_2025_endcap->GetXaxis()->SetTitle("Event size (kB)");
    h_2025_endcap->SetMarkerColor(kRed);
    h_2025_endcap->SetLineColor(kRed);
    h_2025_endcap->SetMarkerStyle(kFullCircle);

    // h_2024_endcap_zs9->SetMarkerColor(kRed);
    // h_2024_endcap_zs9->SetLineColor(kRed);
    // h_2024_endcap_zs9->SetMarkerStyle(kOpenCircle);

    // h_2024_endcap_zs10->SetMarkerColor(kGreen);
    // h_2024_endcap_zs10->SetLineColor(kGreen);
    // h_2024_endcap_zs10->SetMarkerStyle(kOpenSquare);

    // Draw 
    TLegend *leg_endcap = new TLegend(0.5, 0.6, 0.8, 0.85);
    leg_endcap->SetBorderSize(0);
    leg_endcap->SetHeader("endcap, LTH=8, HTH=10");
    leg_endcap->AddEntry(h_2023_endcap, Form("2023 (375513) ZS=8; #mu=%.2f", h_2023_endcap->GetMean()));
    leg_endcap->AddEntry(h_2024_endcap, Form("2024 (387855) ZS=8; #mu=%.2f", h_2024_endcap->GetMean()));
    leg_endcap->AddEntry(h_2025_endcap, Form("2025 (399002) ZS=8; #mu=%.2f", h_2025_endcap->GetMean()));
    // leg_endcap->AddEntry(h_2024_endcap_zs9, Form("387855 ZS=9; #mu=%.2f", h_2024_endcap_zs9->GetMean()));
    // leg_endcap->AddEntry(h_2024_endcap_zs10, Form("387855 ZS=10; #mu=%.2f", h_2024_endcap_zs10->GetMean()));

    TCanvas *c_endcap = new TCanvas("c_endcap", "", 800, 600);
    c_endcap->SetLeftMargin(0.15);
    c_endcap->SetTopMargin(0.01);
    c_endcap->SetLogy();
    h_2023_endcap->Draw();
    h_2024_endcap->Draw("same");
    h_2025_endcap->Draw("same");
    // h_2024_endcap_zs9->Draw("same");
    // h_2024_endcap_zs10->Draw("same");
    leg_endcap->Draw();
}