//sseeds 10.16.23 - script to aggregate various emin varied tdiff plots and obtain SNR vs emin

#include "TH1D.h"
#include "TCanvas.h"

const int nev[8] = {10012,10013,10014,10015,10016,10017,10018,10019};
const double emin[8] = {0.050,0.020,0.010,0.250,0.005,0.100,0.150,0.200};
const double tdiff_lcut = -9.;
const double tdiff_hcut = 9.;

void DrawAndIntegrate(TH1D* histogram, double x1, double x2, TCanvas* canvas, double &intA, double &intB, double &nev) {
  // Check if the histogram is valid
  if (!histogram) {
    std::cerr << "Invalid histogram!" << std::endl;
    return;
  }

  //Get number of events
  nev = histogram->GetEntries();

  // Create vertical lines
  TLine *line1 = new TLine(x1, 0, x1, histogram->GetBinContent(histogram->FindBin(x1)));
  line1->SetLineColor(kGreen);
  line1->SetLineWidth(2);
  TLine *line2 = new TLine(x2, 0, x2, histogram->GetBinContent(histogram->FindBin(x2)));
  line2->SetLineColor(kGreen);
  line2->SetLineWidth(2);

  // Calculate integrals within and outside the lines
  double integralWithin = histogram->Integral(histogram->FindBin(x1), histogram->FindBin(x2));
  intA = integralWithin;
  double integralOutside = histogram->Integral(1, histogram->GetNbinsX()) - integralWithin;
  intB = integralOutside;

  // Create a legend
  auto legend = new TLegend(0.7, 0.7, 0.89, 0.89);
  legend->AddEntry(line1, Form("Integral Within (%.2f to %.2f): %.2f", x1, x2, integralWithin), "l");
  legend->AddEntry(line2, Form("Integral Outside: %.2f", integralOutside), "l");

  // Draw the histogram and lines on the canvas
  canvas->cd();
  histogram->Draw();
  line1->Draw("same");
  line2->Draw("same");

  // Draw the legend
  legend->Draw();

  canvas->Update();
}

void test_tmax_agg(){
  
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);

  TCanvas *c1 = new TCanvas("c1",Form("tdiff vs emin=%0.3f, nev=%d",emin[0],nev[0]),2200,1200);
  TCanvas *c2 = new TCanvas("c2",Form("tdiff vs emin=%0.3f, nev=%d",emin[1],nev[1]),2200,1200);
  TCanvas *c3 = new TCanvas("c3",Form("tdiff vs emin=%0.3f, nev=%d",emin[2],nev[2]),2200,1200);
  TCanvas *c4 = new TCanvas("c4",Form("tdiff vs emin=%0.3f, nev=%d",emin[3],nev[3]),2200,1200);
  TCanvas *c5 = new TCanvas("c5",Form("tdiff vs emin=%0.3f, nev=%d",emin[4],nev[4]),2200,1200);
  TCanvas *c6 = new TCanvas("c6",Form("tdiff vs emin=%0.3f, nev=%d",emin[5],nev[5]),2200,1200);
  TCanvas *c7 = new TCanvas("c7",Form("tdiff vs emin=%0.3f, nev=%d",emin[6],nev[6]),2200,1200);
  TCanvas *c8 = new TCanvas("c8",Form("tdiff vs emin=%0.3f, nev=%d",emin[7],nev[7]),2200,1200);

  vector<string> file_paths;

  for( int i=0; i<8; ++i )
    file_paths.push_back(Form("root_out/tmaxcompOUT_control_run12916_Nev%d.root",nev[i]));

  // open files
  TFile *f0 = TFile::Open(file_paths[0].c_str());
  TFile *f1 = TFile::Open(file_paths[1].c_str());
  TFile *f2 = TFile::Open(file_paths[2].c_str());
  TFile *f3 = TFile::Open(file_paths[3].c_str());
  TFile *f4 = TFile::Open(file_paths[4].c_str());
  TFile *f5 = TFile::Open(file_paths[5].c_str());
  TFile *f6 = TFile::Open(file_paths[6].c_str());
  TFile *f7 = TFile::Open(file_paths[7].c_str());

  //get tdiff histos
  TH1D *hf0_tdiff = (TH1D*)f0->Get("ht_all_diff");
  hf0_tdiff->SetTitle(Form("adct pblk-sblk, emin=%0.3f",emin[0]));

  TH1D *hf1_tdiff = (TH1D*)f1->Get("ht_all_diff");
  hf1_tdiff->SetTitle(Form("adct pblk-sblk, emin=%0.3f",emin[1]));

  TH1D *hf2_tdiff = (TH1D*)f2->Get("ht_all_diff");
  hf2_tdiff->SetTitle(Form("adct pblk-sblk, emin=%0.3f",emin[2]));

  TH1D *hf3_tdiff = (TH1D*)f3->Get("ht_all_diff");
  hf3_tdiff->SetTitle(Form("adct pblk-sblk, emin=%0.3f",emin[3]));

  TH1D *hf4_tdiff = (TH1D*)f4->Get("ht_all_diff");
  hf4_tdiff->SetTitle(Form("adct pblk-sblk, emin=%0.3f",emin[4]));

  TH1D *hf5_tdiff = (TH1D*)f5->Get("ht_all_diff");
  hf5_tdiff->SetTitle(Form("adct pblk-sblk, emin=%0.3f",emin[5]));

  TH1D *hf6_tdiff = (TH1D*)f6->Get("ht_all_diff");
  hf6_tdiff->SetTitle(Form("adct pblk-sblk, emin=%0.3f",emin[6]));

  TH1D *hf7_tdiff = (TH1D*)f7->Get("ht_all_diff");
  hf7_tdiff->SetTitle(Form("adct pblk-sblk, emin=%0.3f",emin[7]));

  //set up integral doubles for TGraph
  double f0_in;
  double f0_out;
  double f1_in;
  double f1_out;
  double f2_in;
  double f2_out;
  double f3_in;
  double f3_out;
  double f4_in;
  double f4_out;
  double f5_in;
  double f5_out;
  double f6_in;
  double f6_out;
  double f7_in;
  double f7_out;

  //set up NEV doubles for TGraph
  double f0_nev;
  double f1_nev;
  double f2_nev;
  double f3_nev;
  double f4_nev;
  double f5_nev;
  double f6_nev;
  double f7_nev;

  c1->cd();
  DrawAndIntegrate(hf0_tdiff,tdiff_lcut,tdiff_hcut,c1,f0_in,f0_out,f0_nev);
  c2->cd();
  DrawAndIntegrate(hf1_tdiff,tdiff_lcut,tdiff_hcut,c2,f1_in,f1_out,f1_nev);
  c3->cd();
  DrawAndIntegrate(hf2_tdiff,tdiff_lcut,tdiff_hcut,c3,f2_in,f2_out,f2_nev);
  c4->cd();
  DrawAndIntegrate(hf3_tdiff,tdiff_lcut,tdiff_hcut,c4,f3_in,f3_out,f3_nev);
  c5->cd();
  DrawAndIntegrate(hf4_tdiff,tdiff_lcut,tdiff_hcut,c5,f4_in,f4_out,f4_nev);
  c6->cd();
  DrawAndIntegrate(hf5_tdiff,tdiff_lcut,tdiff_hcut,c6,f5_in,f5_out,f5_nev);
  c7->cd();
  DrawAndIntegrate(hf6_tdiff,tdiff_lcut,tdiff_hcut,c7,f6_in,f6_out,f6_nev);
  c8->cd();
  DrawAndIntegrate(hf7_tdiff,tdiff_lcut,tdiff_hcut,c8,f7_in,f7_out,f7_nev);

  //calculate SNR
  double f0_snr = f0_in/f0_out;
  double f1_snr = f1_in/f1_out;
  double f2_snr = f2_in/f2_out;
  double f3_snr = f3_in/f3_out;
  double f4_snr = f4_in/f4_out;
  double f5_snr = f5_in/f5_out;
  double f6_snr = f6_in/f6_out;
  double f7_snr = f7_in/f7_out;

  double snr_max=0;
  if(snr_max<f0_snr)
    snr_max=f0_snr;
  if(snr_max<f1_snr)
    snr_max=f1_snr;
  if(snr_max<f2_snr)
    snr_max=f2_snr;
  if(snr_max<f3_snr)
    snr_max=f3_snr;
  if(snr_max<f4_snr)
    snr_max=f4_snr;
  if(snr_max<f5_snr)
    snr_max=f5_snr;
  if(snr_max<f6_snr)
    snr_max=f6_snr;
  if(snr_max<f7_snr)
    snr_max=f7_snr;

  //draw the canvas with the tgraph
  TCanvas *c9 = new TCanvas("c9","tgraph",2200,1200);
  c9->cd();

  TGraph* graph = new TGraph();
  graph->SetPoint(0, emin[0], f0_snr);
  graph->SetPoint(1, emin[1], f1_snr);
  graph->SetPoint(2, emin[2], f2_snr);
  graph->SetPoint(3, emin[3], f3_snr);
  graph->SetPoint(4, emin[4], f4_snr);
  graph->SetPoint(5, emin[5], f5_snr);
  graph->SetPoint(6, emin[6], f6_snr);
  graph->SetPoint(7, emin[7], f7_snr);
  
  // Customize the point size and style
  graph->SetMarkerStyle(20); // Set the point style (20 is a circle)
  graph->SetMarkerSize(1.5); // Set the point size
  graph->GetXaxis()->SetTitle("emin setting (both seed and blk) in GeV");
  graph->GetYaxis()->SetTitle("Signal To Noise");
  graph->GetYaxis()->SetRangeUser(0,1.1*snr_max);

  graph->Draw("AP");

  c9->Update();

  //scale hint1 to the pad coordinates

  double nev_max=0;
  if(nev_max<f0_nev)
    nev_max=f0_nev;
  if(nev_max<f1_nev)
    nev_max=f1_nev;
  if(nev_max<f2_nev)
    nev_max=f2_nev;
  if(nev_max<f3_nev)
    nev_max=f3_nev;
  if(nev_max<f4_nev)
    nev_max=f4_nev;
  if(nev_max<f5_nev)
    nev_max=f5_nev;
  if(nev_max<f6_nev)
    nev_max=f6_nev;
  if(nev_max<f7_nev)
    nev_max=f7_nev;

  double rightmax = 1.1*nev_max;
  double scale = gPad->GetUymax()/rightmax;

  f0_nev*=scale;
  f1_nev*=scale;
  f2_nev*=scale;
  f3_nev*=scale;
  f4_nev*=scale;
  f5_nev*=scale;
  f6_nev*=scale;
  f7_nev*=scale;

  TGraph* graph_overlay = new TGraph();
  graph_overlay->SetPoint(0, emin[0], f0_nev);
  graph_overlay->SetPoint(1, emin[1], f1_nev);
  graph_overlay->SetPoint(2, emin[2], f2_nev);
  graph_overlay->SetPoint(3, emin[3], f3_nev);
  graph_overlay->SetPoint(4, emin[4], f4_nev);
  graph_overlay->SetPoint(5, emin[5], f5_nev);
  graph_overlay->SetPoint(6, emin[6], f6_nev);
  graph_overlay->SetPoint(7, emin[7], f7_nev);

  graph_overlay->SetMarkerStyle(21); // Square
  graph_overlay->SetMarkerSize(1.5);
  graph_overlay->SetMarkerColor(kRed);

  graph_overlay->Draw("P same");

  //draw an axis on the right side
  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
			    gPad->GetUxmax(),gPad->GetUymax(),0,rightmax,510,"+L");
  axis->SetLineColor(kRed);
  axis->SetLabelColor(kRed);
  //axis->SetTextSize(0.002);
  axis->SetTextColor(kRed);
  axis->SetTitle("N Events");
  axis->SetTitleOffset(1.1);

  axis->Draw();

  c9->Update();

}
