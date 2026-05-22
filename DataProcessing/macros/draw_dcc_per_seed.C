void draw_dcc_per_seed()
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    // Float_t text_size = 28.;
    // gStyle->SetTextSize(text_size);
    // gStyle->SetLegendTextSize(text_size);
    // gStyle->SetLabelSize(text_size, "XYZ");
    // gStyle->SetTitleSize(text_size, "XYZ");

    // ---------- EB files --------------
    TFile *fin_R375513_HFandZDC = new TFile("/afs/cern.ch/user/l/lkalipol/private/ecal-hin/CMSSW_13_2_4/src/DataProcessing/jobs/DoSR_LTH4p0_HTH8p0_NoZS_MIEB8p0_MIEE8p0/srvalid_HI2023_R375513.root");
    TFile *fin_R387855_HF = new TFile("../jobs/DoSR_LTH4p0_HTH8p0_NoZS_MIEB8p0_MIEE8p0/srvalid_hi2024_r387855_HF.root");
    TFile *fin_R387855_HFandZDC = new TFile("../jobs/DoSR_LTH4p0_HTH8p0_NoZS_MIEB8p0_MIEE8p0/srvalid_hi2024_r387855_HFandZDC.root");
    TFile *fin_R387973_HF = new TFile("../jobs/DoSR_LTH4p0_HTH8p0_NoZS_MIEB8p0_MIEE8p0/srvalid_hi2024_r387973_HF.root");
    TFile *fin_R387973_HFandZDC = new TFile("../jobs/DoSR_LTH4p0_HTH8p0_NoZS_MIEB8p0_MIEE8p0/srvalid_hi2024_r387973_HFandZDC.root");
    TFile *fin_R388006_HFandZDC = new TFile("../jobs/DoSR_LTH4p0_HTH8p0_NoZS_MIEB8p0_MIEE8p0/srvalid_hi2024_r388006.root");
    TFile *fin_R388122_HFandZDC = new TFile("../jobs/DoSR_LTH4p0_HTH8p0_NoZS_MIEB8p0_MIEE8p0/srvalid_hi2024_r388122.root");

    // ---------- EE files --------------
    TFile *fin_R375513_HFandZDC_ee = new TFile("/afs/cern.ch/user/l/lkalipol/private/ecal-hin/CMSSW_13_2_4/src/DataProcessing/jobs/DoSR_LTH8p0_HTH10p0_NoZS_MIEB8p0_MIEE8p0/srvalid_HI2023_R375513.root");
    TFile *fin_R387855_HF_ee = new TFile("../jobs/DoSR_LTH8p0_HTH10p0_NoZS_MIEB8p0_MIEE8p0/srvalid_hi2024_r387855_HF.root");
    TFile *fin_R387855_HFandZDC_ee = new TFile("../jobs/DoSR_LTH8p0_HTH10p0_NoZS_MIEB8p0_MIEE8p0/srvalid_hi2024_r387855_HFandZDC.root");
    TFile *fin_R387973_HF_ee = new TFile("../jobs/DoSR_LTH8p0_HTH10p0_NoZS_MIEB8p0_MIEE8p0/srvalid_hi2024_r387973_HF.root");
    TFile *fin_R387973_HFandZDC_ee = new TFile("../jobs/DoSR_LTH8p0_HTH10p0_NoZS_MIEB8p0_MIEE8p0/srvalid_hi2024_r387973_HFandZDC.root");
    TFile *fin_R388006_HFandZDC_ee = new TFile("../jobs/DoSR_LTH8p0_HTH10p0_NoZS_MIEB8p0_MIEE8p0/srvalid_hi2024_r388006.root");
    TFile *fin_R388122_HFandZDC_ee = new TFile("../jobs/DoSR_LTH8p0_HTH10p0_NoZS_MIEB8p0_MIEE8p0/srvalid_hi2024_r388122.root");

    // ----------- EB histograms ----------
    TH1F *h_960b_hfANDzdc_eb = (TH1F *) fin_R375513_HFandZDC->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_960b_hfANDzdc_eb");
    TH1F *h_58b_hf_eb = (TH1F *) fin_R387855_HF->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_58b_hf_eb");
    TH1F *h_58b_hfANDzdc_eb = (TH1F *) fin_R387855_HFandZDC->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_58b_hfANDzdc_eb");
    TH1F *h_640b_hf_eb = (TH1F *) fin_R387973_HF->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_640b_hf_eb");
    TH1F *h_640b_hfANDzdc_eb = (TH1F *) fin_R387973_HFandZDC->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_640b_hfANDzdc_eb");
    TH1F *h_1032b_hfANDzdc_eb = (TH1F *) fin_R388006_HFandZDC->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_1032b_hfANDzdc_eb");

    // ------------ EE histograms ----------
    TH1F *h_960b_hfANDzdc_ee = (TH1F *) fin_R375513_HFandZDC_ee->Get("ecalSelectiveReadoutValidation/hDccVol5")->Clone("h_960b_hfANDzdc_ee");
    TH1F *h_58b_hf_ee = (TH1F *) fin_R387855_HF_ee->Get("ecalSelectiveReadoutValidation/hDccVol5")->Clone("h_58b_hf_ee");
    TH1F *h_58b_hfANDzdc_ee = (TH1F *) fin_R387855_HFandZDC_ee->Get("ecalSelectiveReadoutValidation/hDccVol5")->Clone("h_58b_hfANDzdc_ee");
    TH1F *h_640b_hf_ee = (TH1F *) fin_R387973_HF_ee->Get("ecalSelectiveReadoutValidation/hDccVol5")->Clone("h_640b_hf_ee");
    TH1F *h_640b_hfANDzdc_ee = (TH1F *) fin_R387973_HFandZDC_ee->Get("ecalSelectiveReadoutValidation/hDccVol5")->Clone("h_640b_hfANDzdc_ee");
    TH1F *h_1032b_hfANDzdc_ee = (TH1F *) fin_R388006_HFandZDC_ee->Get("ecalSelectiveReadoutValidation/hDccVol5")->Clone("h_1032b_hfANDzdc_ee");

    // Normalize
    for (auto h : {
        h_960b_hfANDzdc_eb,
        h_58b_hf_eb, h_58b_hfANDzdc_eb,
        h_640b_hf_eb, h_640b_hfANDzdc_eb,
        h_1032b_hfANDzdc_eb,
        h_960b_hfANDzdc_ee,
        h_58b_hf_ee, h_58b_hfANDzdc_ee,
        h_640b_hf_ee, h_640b_hfANDzdc_ee,
        h_1032b_hfANDzdc_ee
        }) {
            // h->GetXaxis()->SetRange(1, h->GetNbinsX()+1);
            // std::cout << h->GetBinContent(101) << std::endl;
            h->Scale(1/h->Integral());
            h->GetYaxis()->SetRangeUser(1e-4, 1.1);
    }

    // --------- Format EB ----------
    h_960b_hfANDzdc_eb->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_960b_hfANDzdc_eb->GetXaxis()->SetTitle("Event size (kB)");
    h_960b_hfANDzdc_eb->SetMarkerColor(kViolet);
    h_960b_hfANDzdc_eb->SetLineColor(kViolet);
    h_960b_hfANDzdc_eb->SetMarkerStyle(0);
    h_960b_hfANDzdc_eb->SetLineWidth(4);

    h_58b_hf_eb->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_58b_hf_eb->GetXaxis()->SetTitle("Event size (kB)");
    h_58b_hf_eb->SetMarkerColor(kBlue);
    h_58b_hf_eb->SetLineColor(kBlue);
    h_58b_hf_eb->SetMarkerStyle(kOpenSquare);

    h_58b_hfANDzdc_eb->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_58b_hfANDzdc_eb->GetXaxis()->SetTitle("Event size (kB)");
    h_58b_hfANDzdc_eb->SetMarkerColor(kRed);
    h_58b_hfANDzdc_eb->SetLineColor(kRed);
    h_58b_hfANDzdc_eb->SetMarkerStyle(kFullSquare);

    h_640b_hf_eb->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_640b_hf_eb->GetXaxis()->SetTitle("Event size (kB)");
    h_640b_hf_eb->SetMarkerColor(kGreen);
    h_640b_hf_eb->SetLineColor(kGreen);
    h_640b_hf_eb->SetMarkerStyle(kOpenCircle);

    h_640b_hfANDzdc_eb->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_640b_hfANDzdc_eb->GetXaxis()->SetTitle("Event size (kB)");
    h_640b_hfANDzdc_eb->SetMarkerColor(kBlack);
    h_640b_hfANDzdc_eb->SetLineColor(kBlack);
    h_640b_hfANDzdc_eb->SetMarkerStyle(kFullCircle);

    h_1032b_hfANDzdc_eb->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_1032b_hfANDzdc_eb->GetXaxis()->SetTitle("Event size (kB)");
    h_1032b_hfANDzdc_eb->SetMarkerColor(kOrange);
    h_1032b_hfANDzdc_eb->SetLineColor(kOrange);
    h_1032b_hfANDzdc_eb->SetMarkerStyle(kFullStar);

    // ----------- Draw EB ----------------
    TLegend *leg_barrel = new TLegend(0.45, 0.55, 0.85, 0.95);
    leg_barrel->SetBorderSize(0);
    leg_barrel->SetHeader("Barrel");
    leg_barrel->AddEntry(h_960b_hfANDzdc_eb, Form("2023 960b, HF and ZDC; #mu=%.2f", h_960b_hfANDzdc_eb->GetMean()));
    leg_barrel->AddEntry(h_58b_hf_eb, Form("58b, HF only; #mu=%.2f", h_58b_hf_eb->GetMean()));
    leg_barrel->AddEntry(h_58b_hfANDzdc_eb, Form("58b, HF and ZDC; #mu=%.2f", h_58b_hfANDzdc_eb->GetMean()));
    leg_barrel->AddEntry(h_640b_hf_eb, Form("640b, HF only; #mu=%.2f", h_640b_hf_eb->GetMean()));
    leg_barrel->AddEntry(h_640b_hfANDzdc_eb, Form("640b, HF and ZDC; #mu=%.2f", h_640b_hfANDzdc_eb->GetMean()));
    leg_barrel->AddEntry(h_1032b_hfANDzdc_eb, Form("1032b, HF and ZDC; #mu=%.2f", h_1032b_hfANDzdc_eb->GetMean()));

    TCanvas *c_barrel = new TCanvas("c_barrel", "", 800, 600);
    c_barrel->SetLeftMargin(0.15);
    c_barrel->SetTopMargin(0.01);
    c_barrel->SetLogy();
    h_960b_hfANDzdc_eb->Draw("hist");
    h_58b_hf_eb->Draw("pe same");
    h_58b_hfANDzdc_eb->Draw("pe same");
    h_640b_hf_eb->Draw("pe same");
    h_640b_hfANDzdc_eb->Draw("pe same");
    h_1032b_hfANDzdc_eb->Draw("pe same");
    leg_barrel->Draw();

    // --------- Format EE ----------
    h_960b_hfANDzdc_ee->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_960b_hfANDzdc_ee->GetXaxis()->SetTitle("Event size (kB)");
    h_960b_hfANDzdc_ee->SetMarkerColor(kViolet);
    h_960b_hfANDzdc_ee->SetLineColor(kViolet);
    h_960b_hfANDzdc_ee->SetMarkerStyle(0);
    h_960b_hfANDzdc_ee->SetLineWidth(4);

    h_58b_hf_ee->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_58b_hf_ee->GetXaxis()->SetTitle("Event size (kB)");
    h_58b_hf_ee->SetMarkerColor(kBlue);
    h_58b_hf_ee->SetLineColor(kBlue);
    h_58b_hf_ee->SetMarkerStyle(kOpenSquare);

    h_58b_hfANDzdc_ee->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_58b_hfANDzdc_ee->GetXaxis()->SetTitle("Event size (kB)");
    h_58b_hfANDzdc_ee->SetMarkerColor(kRed);
    h_58b_hfANDzdc_ee->SetLineColor(kRed);
    h_58b_hfANDzdc_ee->SetMarkerStyle(kFullSquare);

    h_640b_hf_ee->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_640b_hf_ee->GetXaxis()->SetTitle("Event size (kB)");
    h_640b_hf_ee->SetMarkerColor(kGreen);
    h_640b_hf_ee->SetLineColor(kGreen);
    h_640b_hf_ee->SetMarkerStyle(kOpenCircle);

    h_640b_hfANDzdc_ee->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_640b_hfANDzdc_ee->GetXaxis()->SetTitle("Event size (kB)");
    h_640b_hfANDzdc_ee->SetMarkerColor(kBlack);
    h_640b_hfANDzdc_ee->SetLineColor(kBlack);
    h_640b_hfANDzdc_ee->SetMarkerStyle(kFullCircle);

    h_1032b_hfANDzdc_ee->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_1032b_hfANDzdc_ee->GetXaxis()->SetTitle("Event size (kB)");
    h_1032b_hfANDzdc_ee->SetMarkerColor(kOrange);
    h_1032b_hfANDzdc_ee->SetLineColor(kOrange);
    h_1032b_hfANDzdc_ee->SetMarkerStyle(kFullStar);

    // ----------- Draw EE ----------------
    TLegend *leg_endcap = new TLegend(0.45, 0.55, 0.85, 0.95);
    leg_endcap->SetBorderSize(0);
    leg_endcap->SetHeader("Endcap");
    leg_endcap->AddEntry(h_960b_hfANDzdc_ee, Form("2023 960b, HF and ZDC; #mu=%.2f", h_960b_hfANDzdc_ee->GetMean()));
    leg_endcap->AddEntry(h_58b_hf_ee, Form("58b, HF only; #mu=%.2f", h_58b_hf_ee->GetMean()));
    leg_endcap->AddEntry(h_58b_hfANDzdc_ee, Form("58b, HF and ZDC; #mu=%.2f", h_58b_hfANDzdc_ee->GetMean()));
    leg_endcap->AddEntry(h_640b_hf_ee, Form("640b, HF only; #mu=%.2f", h_640b_hf_ee->GetMean()));
    leg_endcap->AddEntry(h_640b_hfANDzdc_ee, Form("640b, HF and ZDC; #mu=%.2f", h_640b_hfANDzdc_ee->GetMean()));
    leg_endcap->AddEntry(h_1032b_hfANDzdc_ee, Form("1032b, HF and ZDC; #mu=%.2f", h_1032b_hfANDzdc_ee->GetMean()));

    TCanvas *c_endcap = new TCanvas("c_endcap", "", 800, 600);
    c_endcap->SetLeftMargin(0.15);
    c_endcap->SetTopMargin(0.01);
    c_endcap->SetLogy();
    h_960b_hfANDzdc_ee->Draw("hist same");
    h_58b_hf_ee->Draw("pe same");
    h_58b_hfANDzdc_ee->Draw("pe same");
    h_640b_hf_ee->Draw("pe same");
    h_640b_hfANDzdc_ee->Draw("pe same");
    h_1032b_hfANDzdc_ee->Draw("pe same");
    leg_endcap->Draw();
}