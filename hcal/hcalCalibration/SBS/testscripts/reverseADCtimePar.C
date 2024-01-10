//sseeds 10.23.23 Script to add the adct offsets instead of subtract them

#include "../include/sbs.h"
#include <vector>
#include <iostream>
#include <algorithm>

void coutParams( Double_t params[hcal::maxHCalChan] ){
  Int_t col = 0;
  for( Int_t i=0; i<hcal::maxHCalChan; i++ ){
    std::cout << params[i] << " ";
    if(col==11){
      col=0;
      std::cout << std::endl;
      continue;
    }
    col++;
  }
  std::cout << std::endl << std::endl;
}

void subtractOffsets( Double_t array1[hcal::maxHCalChan], Double_t array2[hcal::maxHCalChan], Double_t array3[hcal::maxHCalChan] ){
  Int_t col = 0;
  for( Int_t i=0; i<hcal::maxHCalChan; i++ ){
    array3[i] = array1[i] - array2[i];
    std::cout << array3[i] << " ";
    if(col==11){
      col=0;
      std::cout << std::endl;
      continue;
    }
    col++;
  }
  std::cout << std::endl << std::endl;
}

void addOffsets( Double_t array1[hcal::maxHCalChan], Double_t array2[hcal::maxHCalChan], Double_t array3[hcal::maxHCalChan] ){
  Int_t col = 0;
  for( Int_t i=0; i<hcal::maxHCalChan; i++ ){
    array3[i] = array1[i] + array2[i];
    std::cout << array3[i] << " ";
    if(col==11){
      col=0;
      std::cout << std::endl;
      continue;
    }
    col++;
  }
  std::cout << std::endl << std::endl;
}

void writeFixedOffsets( Double_t array1[hcal::maxHCalChan], const char *path, const char *db_var ){
  
  ofstream fo;

  fo.open( path );
  fo << "#fixed adct offsets from testscripts/reverseADCtimePar.C" << std::endl;
  fo << db_var << " =" << std::endl;

  Int_t col = 0;
  for( Int_t i=0; i<hcal::maxHCalChan; i++ ){
    fo << array1[i] << " ";
    if(col==11){
      col=0;
      fo << std::endl;
      continue;
    }
    col++;
  }
  fo << std::endl << std::endl;
}


void reverseADCtimePar(Int_t config = 11, Int_t pass = 1, bool writeopt = true){
  
  //if desired, write the output to the timing/parameters directory
  std::string fixed_offset_path = Form("../timing/parameters/adctoffsets_pass2_fixed_gmn_conf%d_pass%d_12450_to_12860_exclude_0_to_0.txt",config,pass);

  //Get the old offsets and write them to screen
  Double_t new_adct_offsets[hcal::maxHCalChan] = {0.};
  Double_t old_adct_offsets[hcal::maxHCalChan] = {0.};
  Double_t addback_offsets[hcal::maxHCalChan] = {0.};
  Double_t fixed_offsets[hcal::maxHCalChan] = {0.};
  std::string db_adctoffset_variable = "sbs.hcal.adc.timeoffset";

  //std::string db_path = gSystem->Getenv("DB_DIR");
  std::string db_path = "/work/halla/sbs/seeds/alt_sbsreplay/SBS-replay/DB";

  std::string old_db_path = db_path + "/db_sbs.hcal.dat";
  std::string old_adctoffset_timestamp = "none"; //should make sure this is right

  util::readDB( old_db_path, old_adctoffset_timestamp, db_adctoffset_variable, old_adct_offsets );  

  cout << "OLD ADCt offsets from path: " << old_db_path << endl;
  coutParams(old_adct_offsets);

  //Get the new offsets and write them to screen
  std::string new_db_path = "../timing/parameters/adctsetoffsets_gmn_conf11_pass1_12450_to_12860_exclude_0_to_0.txt";
  std::string new_adctoffset_timestamp = "none"; //should make sure this is right

  util::readDB( new_db_path, new_adctoffset_timestamp, db_adctoffset_variable, new_adct_offsets ); 

  cout << "NEW ADCt offsets from path: " << new_db_path << endl;
  coutParams(new_adct_offsets);

  //Perform subtraction
  cout << "Difference between these sets per block: " << endl;
  subtractOffsets(old_adct_offsets,new_adct_offsets,addback_offsets);

  //Add the difference to reverse the calibration about the original point
  cout << "New offset parameters, assuming reversal is needed: " << endl;
  addOffsets(old_adct_offsets,addback_offsets,fixed_offsets);

  if( !writeopt )
    fixed_offset_path = "/dev/null"; //Safety to prevent overwriting constants if quasireplay not checked
  else
    writeFixedOffsets(fixed_offsets,fixed_offset_path.c_str(),db_adctoffset_variable.c_str());

  if( writeopt )
    std::cout << "writeopt true, writing fixed offsets to " << fixed_offset_path << std::endl;
  
}
