#include "../include/sbs.h"
#include <vector>
#include <iostream>
#include <algorithm>

void dbgettest(){
  
  std::string DBpath = gSystem->Getenv("DB_DIR");
  std::string tdcOffsetPath = DBpath + "/db_sbs.hcal.dat";
  std::string tdcCalibPath = DBpath + "/db_sbs.hcal.dat";
  std::string adctOffsetPath = DBpath + "/db_sbs.hcal.dat";
  std::string newtwPath = "../timing/parameters/tdctw_class_gmn_conf4_pass0.txt";
  std::string testPath = "../timing/parameters/test.txt";
  std::string new_adcgain_path = "parameters/adcgaincoeff_w2update_gmn_conf4_pass0.txt";
  std::string again_path = "../energy/parameters/adcgaincoeff_gmn_conf4_pass0.txt";
  std::string new_adctoffset_path = "../timing/parameters/adctsetoffsets_newclus_gmn_conf4_pass0_0_to_0_exclude_0_to_0.txt";

  cout << "Looking for constants in file: " << again_path << endl;

  Double_t calib_const = -1;
  Double_t tdcoffsets[hcal::maxHCalChan] = {0.};
  vector<Double_t> tdctw;
  Double_t again[hcal::maxHCalChan] = {1.};
  Double_t adctoffsets[hcal::maxHCalChan] = {0.};

  std::string none_timestamp = "none";
  std::string current_timestamp = "--------[ 2021-10-21 00:00:00 ]";
  //std::string current_timestamp = "-------[ 2021-10-24 04:30:00 ]";
  std::string tw_timestamp = "-------[ 2021-10-24 04:30:00 ]";
  //std::string current_timestamp = "--------[ 2024-07-26 00:00:00 ]";
  std::string again_timestamp = "-------[ 2021-10-24 04:30:00 ]";
  std::string tdc_calib_type = "sbs.hcal.tdc.calib";
  std::string tdc_offset_type = "sbs.hcal.tdc.offset";
  std::string timewalkA_type = "sbs.hcal.tdc.chan_tw_a";
  std::string timewalkB_type = "sbs.hcal.tdc.chan_tw_b";
  std::string timewalk_type = "sbs.hcal.tdc.tw";
  std::string again_type = "sbs.hcal.adc.gain";
  std::string nev_type = "#Number of events available for calibration";
  std::string adct_type = "(inverted) sbs.hcal.adc.timeoffset";

  //util::readDB(tdcCalibPath,current_timestamp,tdc_calib_type,calib_const);
  //util::readDB(tdcOffsetPath,current_timestamp,tdc_offset_type,tdcoffsets);
  //util::readDB(newtwPath,current_timestamp,timewalkB_type,tdcoffsets);
  //util::readDB(newtwPath,tw_timestamp,timewalk_type,tdctw);
  //util::readDB(again_path,again_timestamp,nev_type,again);
  //util::readDB(new_adctoffset_path,none_timestamp,adct_type,adctoffsets);
  util::readDB(adctOffsetPath,none_timestamp,tdc_offset_type,tdcoffsets);

  //cout << calib_const << endl;

  // for( Int_t i=0; i<tdctw.size(); i++ )
  //   cout << tdctw[i] << " ";
  // cout << endl << endl;

  for( Int_t r=0; r<hcal::maxHCalRows; r++){
    for( Int_t c=0; c<hcal::maxHCalCols; c++){
      Int_t i = r*hcal::maxHCalCols+c;
      cout << adctoffsets[i] << " ";
    }
    cout << endl;
  }

}
