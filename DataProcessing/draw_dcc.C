void draw_dcc()
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
    TH1F *h_eb_target = (TH1F *) fin_LTH4_HTH8->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_eb_target");

    TFile *fin_LTH8_HTH10 = new TFile("jobs/DoSR_LTH8p0_HTH10p0_NoZS_MIEB8p0_MIEE8p0/srvalid_HI2023_R375790.root");
    TH1F *h_eb_safe = (TH1F *) fin_LTH8_HTH10->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_eb_safe");
    TH1F *h_ee_target = (TH1F *) fin_LTH8_HTH10->Get("ecalSelectiveReadoutValidation/hDccVol5")->Clone("h_ee_target");

    TFile *fin_LTH10_HTH14 = new TFile("jobs/DoSR_LTH10p0_HTH14p0_NoZS_MIEB8p0_MIEE8p0/srvalid_HI2023_R375790.root");
    TH1F *h_ee_safe = (TH1F *) fin_LTH10_HTH14->Get("ecalSelectiveReadoutValidation/hDccVol5")->Clone("h_ee_safe");

    // Normalize
    for (auto h : {
        h_eb_target, h_eb_safe,
        h_ee_target, h_ee_safe
        }) {
        h->Scale(1/h->Integral());
    }

    // Format
    h_eb_target->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_eb_target->GetXaxis()->SetTitle("Event size (kB)");
    h_eb_target->SetMarkerColor(kBlue);
    h_eb_target->SetLineColor(kBlue);
    h_eb_target->SetMarkerStyle(kOpenCircle);

    h_eb_safe->SetMarkerColor(kRed);
    h_eb_safe->SetLineColor(kRed);
    h_eb_safe->SetMarkerStyle(kOpenSquare);

    h_ee_target->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_ee_target->GetXaxis()->SetTitle("Event size (kB)");
    h_ee_target->SetMarkerColor(kBlue);
    h_ee_target->SetLineColor(kBlue);
    h_ee_target->SetMarkerStyle(kOpenCircle);

    h_ee_safe->SetMarkerColor(kRed);
    h_ee_safe->SetLineColor(kRed);
    h_ee_safe->SetMarkerStyle(kOpenSquare);

    // Draw 
    TLegend *leg_barrel = new TLegend(0.5, 0.6, 0.8, 0.85);
    leg_barrel->SetBorderSize(0);
    leg_barrel->SetHeader("Barrel");
    leg_barrel->AddEntry(h_eb_target, Form("Target; #mu=%.2f", h_eb_target->GetMean()));
    leg_barrel->AddEntry(h_eb_safe, Form("Safe; #mu=%.2f", h_eb_safe->GetMean()));

    TCanvas *c_barrel = new TCanvas("c_barrel", "", 800, 600);
    c_barrel->SetLeftMargin(0.15);
    c_barrel->SetTopMargin(0.01);
    c_barrel->SetLogy();
    h_eb_target->Draw();
    h_eb_safe->Draw("same");
    leg_barrel->Draw();

    TLegend *leg_endcap = new TLegend(0.5, 0.6, 0.8, 0.85);
    leg_endcap->SetBorderSize(0);
    leg_endcap->SetHeader("Endcap");
    leg_endcap->AddEntry(h_ee_target, Form("Target; #mu=%.2f", h_ee_target->GetMean()));
    leg_endcap->AddEntry(h_ee_safe, Form("Safe; #mu=%.2f", h_ee_safe->GetMean()));

    TCanvas *c_endcap = new TCanvas("c_endcap", "", 800, 600);
    c_endcap->SetLeftMargin(0.15);
    c_endcap->SetTopMargin(0.01);
    c_endcap->SetLogy();
    h_ee_target->Draw();
    h_ee_safe->Draw("same");
    leg_endcap->Draw();
}