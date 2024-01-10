#include <TCanvas.h>
#include <TLatex.h>
#include <vector>
#include <string>
#include <sstream>
#include <TText.h>
#include "../include/sbs.h"

bool verb = false;

// Define the struct
// struct suppset {
//   std::string exper;
//   int kine;
//   int pass;
//   std::string date;
//   int runb;
//   int rune;
//   int runexb;
//   int runexe;
//   std::string targ;
// };

void TimingReport(const std::vector<caldiag>& data, const suppset& supplement) {
    // Create a new ROOT file to save the canvases
    TFile *file = new TFile(".root", "RECREATE");

    for (size_t i = 0; i < data.size(); ++i) {
        // Create a new canvas for each MyData struct
        TCanvas *canvas = new TCanvas(Form("c%zu", i+1), Form("Canvas %zu", i+1), 1000, 900);
        canvas->cd();

        double yPos = 0.9; // Starting Y position for text, top of the canvas
        double yStep = 0.025; // Step size for each line of text
	double textSize = 0.02; //size of text
	double textSize_h = 0.025; //size of header text

        std::stringstream ss;
        ss << "General Set ADCt Alignment Info";
        TLatex *genLatex_l1 = new TLatex(0.1, yPos, ss.str().c_str());
        genLatex_l1->SetTextColor(kRed-5);
        genLatex_l1->SetTextSize(textSize_h);
        genLatex_l1->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "Experiment: " << supplement.exper 
	   << ", Configuration: " << supplement.kine 
	   << ", Pass: " << supplement.pass;
        TLatex *headerLatex_l1 = new TLatex(0.1, yPos, ss.str().c_str());
        headerLatex_l1->SetTextColor(kBlack);
        headerLatex_l1->SetTextSize(textSize);
        headerLatex_l1->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "Creation Date: " << supplement.date;
        TLatex *headerLatex_l2 = new TLatex(0.1, yPos, ss.str().c_str());
        headerLatex_l2->SetTextColor(kBlack);
        headerLatex_l2->SetTextSize(textSize);
        headerLatex_l2->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "Run Range: " << supplement.runb << " - " << supplement.rune;
        TLatex *headerLatex_l3 = new TLatex(0.1, yPos, ss.str().c_str());
        headerLatex_l3->SetTextColor(kBlack);
        headerLatex_l3->SetTextSize(textSize);
        headerLatex_l3->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "Exclusion Range: " << supplement.runexb << " - " << supplement.runexe;
        TLatex *headerLatex_l4 = new TLatex(0.1, yPos, ss.str().c_str());
        headerLatex_l4->SetTextColor(kBlack);
        headerLatex_l4->SetTextSize(textSize);
        headerLatex_l4->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "Target(s) Used: " << supplement.targ;
        TLatex *headerLatex_l5 = new TLatex(0.1, yPos, ss.str().c_str());
        headerLatex_l5->SetTextColor(kBlack);
        headerLatex_l5->SetTextSize(textSize);
        headerLatex_l5->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "";
        TLatex *genLatex_l2 = new TLatex(0.1, yPos, ss.str().c_str());
        genLatex_l2->SetTextColor(kRed-5);
        genLatex_l2->SetTextSize(textSize);
        genLatex_l2->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "Electron Arm Elastic Cuts";
        TLatex *genLatex_l3 = new TLatex(0.1, yPos, ss.str().c_str());
        genLatex_l3->SetTextColor(kRed-5);
        genLatex_l3->SetTextSize(textSize_h);
        genLatex_l3->Draw();
        yPos -= yStep;

        // Now display the data from the caldiag struct
        ss.str("");
        ss << "Target: " << data[i].target;
        TLatex *caldiagLatex_la = new TLatex(0.1, yPos, ss.str().c_str());
        caldiagLatex_la->SetTextColor(kBlue);
        caldiagLatex_la->SetTextSize(textSize);
        caldiagLatex_la->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "SBS Field: " << data[i].field;
        TLatex *caldiagLatex_lb = new TLatex(0.1, yPos, ss.str().c_str());
        caldiagLatex_lb->SetTextColor(kBlue);
        caldiagLatex_lb->SetTextSize(textSize);
        caldiagLatex_lb->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "Global Elastic Cuts: " << data[i].gcut;
        TLatex *caldiagLatex_l1 = new TLatex(0.1, yPos, ss.str().c_str());
        caldiagLatex_l1->SetTextColor(kBlue);
        caldiagLatex_l1->SetTextSize(textSize);
        caldiagLatex_l1->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "W2 Min: " << data[i].W2_min;
        TLatex *caldiagLatex_l2 = new TLatex(0.1, yPos, ss.str().c_str());
        caldiagLatex_l2->SetTextColor(kBlue);
        caldiagLatex_l2->SetTextSize(textSize);
        caldiagLatex_l2->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "W2 Max: " << data[i].W2_max;
        TLatex *caldiagLatex_l3 = new TLatex(0.1, yPos, ss.str().c_str());
        caldiagLatex_l3->SetTextColor(kBlue);
        caldiagLatex_l3->SetTextSize(textSize);
        caldiagLatex_l3->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "";
        TLatex *genLatex_l4 = new TLatex(0.1, yPos, ss.str().c_str());
        genLatex_l4->SetTextColor(kRed-5);
        genLatex_l4->SetTextSize(textSize);
        genLatex_l4->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "Hadron Arm Elastic Cuts";
        TLatex *genLatex_l5 = new TLatex(0.1, yPos, ss.str().c_str());
        genLatex_l5->SetTextColor(kRed-5);
        genLatex_l5->SetTextSize(textSize_h);
        genLatex_l5->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "dy mean: " << data[i].dy0;
        TLatex *caldiagLatex_l4 = new TLatex(0.1, yPos, ss.str().c_str());
        caldiagLatex_l4->SetTextColor(kBlue);
        caldiagLatex_l4->SetTextSize(textSize);
        caldiagLatex_l4->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "dy sigma: " << data[i].dy_sig;
        TLatex *caldiagLatex_l5 = new TLatex(0.1, yPos, ss.str().c_str());
        caldiagLatex_l5->SetTextColor(kBlue);
        caldiagLatex_l5->SetTextSize(textSize);
        caldiagLatex_l5->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "dx proton mean: " << data[i].dx0_p;
        TLatex *caldiagLatex_l6 = new TLatex(0.1, yPos, ss.str().c_str());
        caldiagLatex_l6->SetTextColor(kBlue);
        caldiagLatex_l6->SetTextSize(textSize);
        caldiagLatex_l6->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "dx proton sigma: " << data[i].dx_sig_p;
        TLatex *caldiagLatex_l7 = new TLatex(0.1, yPos, ss.str().c_str());
        caldiagLatex_l7->SetTextColor(kBlue);
        caldiagLatex_l7->SetTextSize(textSize);
        caldiagLatex_l7->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "dx neutron mean: " << data[i].dx0_n;
        TLatex *caldiagLatex_l8 = new TLatex(0.1, yPos, ss.str().c_str());
        caldiagLatex_l8->SetTextColor(kBlue);
        caldiagLatex_l8->SetTextSize(textSize);
        caldiagLatex_l8->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "dx neutron sigma: " << data[i].dx_sig_n;
        TLatex *caldiagLatex_l9 = new TLatex(0.1, yPos, ss.str().c_str());
        caldiagLatex_l9->SetTextColor(kBlue);
        caldiagLatex_l9->SetTextSize(textSize);
        caldiagLatex_l9->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "";
        TLatex *genLatex_l6 = new TLatex(0.1, yPos, ss.str().c_str());
        genLatex_l6->SetTextColor(kRed-5);
        genLatex_l6->SetTextSize(textSize);
        genLatex_l6->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "Other Cuts";
        TLatex *genLatex_l7 = new TLatex(0.1, yPos, ss.str().c_str());
        genLatex_l7->SetTextColor(kRed-5);
        genLatex_l7->SetTextSize(textSize);
        genLatex_l7->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "Minimum Events Per Cell: " << data[i].min_ev;
        TLatex *caldiagLatex_l10 = new TLatex(0.1, yPos, ss.str().c_str());
        caldiagLatex_l10->SetTextColor(kBlue);
        caldiagLatex_l10->SetTextSize(textSize);
        caldiagLatex_l10->Draw();
        yPos -= yStep;

        // Write the canvas to the file
        canvas->Write();
        canvas->Update();
    }

    // Close the file
    file->Close();
}

void saveManyCanvas(const char *experiment = "gmn", Int_t kine=4, Int_t pass=1) {

  // Get the date
  std::string date = util::getDate();

  std::vector<caldiag> data;

  std::string struct_dir = Form("../config/%s_p%d/",experiment,pass);

  if(kine!=4){
    cout << "ERROR: Test only works with kine==4 on current version" << endl;
    return;
  }
  
  std::string exper = experiment;
  int runb = 0;
  int rune = 0;
  int runexb = 0;
  int runexe = 0;
  std::string targ = "All Available";
  suppset supplement;

  supplement.exper = exper;
  supplement.kine = kine;
  supplement.pass = pass;
  supplement.date = date;
  supplement.runb = runb;
  supplement.rune = rune;
  supplement.runexb = runexb;
  supplement.runexe = runexe;
  supplement.targ = targ;
  

  //Hardcode SBS4 field sets and targets
  vector<caldiag> cut_f1;
  Int_t mag_f1 = 0;
  std::string targ_f1 = "lh2";
  util::ReadDiagnosticCutList(struct_dir,experiment,kine,targ_f1,mag_f1,verb,cut_f1); 
  data.push_back(cut_f1[0]);

  vector<caldiag> cut_f2;
  Int_t mag_f2 = 30;
  std::string targ_f2 = "lh2";
  util::ReadDiagnosticCutList(struct_dir,experiment,kine,targ_f2,mag_f2,verb,cut_f2); 
  data.push_back(cut_f2[0]);

  vector<caldiag> cut_f3;
  Int_t mag_f3 = 50;
  std::string targ_f3 = "lh2";
  util::ReadDiagnosticCutList(struct_dir,experiment,kine,targ_f3,mag_f3,verb,cut_f3); 
  data.push_back(cut_f3[0]);

  vector<caldiag> cut_f4;
  Int_t mag_f4 = 0;
  std::string targ_f4 = "ld2";
  util::ReadDiagnosticCutList(struct_dir,experiment,kine,targ_f4,mag_f4,verb,cut_f4); 
  data.push_back(cut_f4[0]);

  vector<caldiag> cut_f5;
  Int_t mag_f5 = 30;
  std::string targ_f5 = "ld2";
  util::ReadDiagnosticCutList(struct_dir,experiment,kine,targ_f5,mag_f5,verb,cut_f5); 
  data.push_back(cut_f5[0]);

  vector<caldiag> cut_f6;
  Int_t mag_f6 = 50;
  std::string targ_f6 = "ld2";
  util::ReadDiagnosticCutList(struct_dir,experiment,kine,targ_f6,mag_f6,verb,cut_f6);  
  data.push_back(cut_f6[0]);

  TimingReport(data,supplement);
}
