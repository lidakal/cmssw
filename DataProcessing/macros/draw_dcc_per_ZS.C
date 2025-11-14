void draw_dcc_per_ZS()
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    // Float_t text_size = 28.;
    // gStyle->SetTextSize(text_size);
    // gStyle->SetLegendTextSize(text_size);
    // gStyle->SetLabelSize(text_size, "XYZ");
    // gStyle->SetTitleSize(text_size, "XYZ");

    TFile *fin_zs8 = new TFile("jobs/DoSR_LTH4p0_HTH8p0_NoZS_MIEB8p0_MIEE8p0/srvalid_hi2024setup_r387377.root");
    TH1F *h_zs8 = (TH1F *) fin_zs8->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_zs8");

    TFile *fin_zs9 = new TFile("jobs/DoSR_LTH4p0_HTH8p0_NoZS_MIEB9p0_MIEE8p0/srvalid_hi2024setup_r387377.root");
    TH1F *h_zs9 = (TH1F *) fin_zs9->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_zs9");

    TFile *fin_zs10 = new TFile("jobs/DoSR_LTH4p0_HTH8p0_NoZS_MIEB10p0_MIEE8p0/srvalid_hi2024setup_r387377.root");
    TH1F *h_zs10 = (TH1F *) fin_zs10->Get("ecalSelectiveReadoutValidation/hDccVol20")->Clone("h_zs10");

    // Normalize
    for (auto h : {
        h_zs8, h_zs9, h_zs10
        }) {
            h->GetXaxis()->SetRange(1, h->GetNbinsX()+1);
            std::cout << h->GetBinContent(101) << std::endl;
            h->Scale(1/h->Integral());
    }

    // Format
    h_zs8->GetYaxis()->SetTitle("1/N_{tot} * N_{evt}");
    h_zs8->GetXaxis()->SetTitle("Event size (kB)");
    h_zs8->SetMarkerColor(kBlack);
    h_zs8->SetLineColor(kBlack);
    h_zs8->SetMarkerStyle(kFullCircle);

    h_zs9->SetMarkerColor(kRed);
    h_zs9->SetLineColor(kRed);
    h_zs9->SetMarkerStyle(kOpenCircle);

    h_zs10->SetMarkerColor(kBlue);
    h_zs10->SetLineColor(kBlue);
    h_zs10->SetMarkerStyle(kOpenSquare);

    // Draw 
    TLegend *leg_barrel = new TLegend(0.5, 0.6, 0.8, 0.85);
    leg_barrel->SetBorderSize(0);
    leg_barrel->SetHeader("Barrel, LTH=4, HTH=8");
    leg_barrel->AddEntry(h_zs8, Form("ZS=8; #mu=%.2f", h_zs8->GetMean()));
    leg_barrel->AddEntry(h_zs9, Form("ZS=9; #mu=%.2f", h_zs9->GetMean()));
    leg_barrel->AddEntry(h_zs10, Form("ZS=10; #mu=%.2f", h_zs10->GetMean()));

    TCanvas *c_barrel = new TCanvas("c_barrel", "", 800, 600);
    c_barrel->SetLeftMargin(0.15);
    c_barrel->SetTopMargin(0.01);
    // c_barrel->SetLogy();
    h_zs8->Draw();
    h_zs9->Draw("same");
    h_zs10->Draw("same");
    leg_barrel->Draw();
}