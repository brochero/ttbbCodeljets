// https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagCalibration
#include "BTagSFUtil.h"
#include "BTagCalibrationStandalone.h"

BTagSFUtil::BTagSFUtil(string BTagAlgorithm, TString OperatingPoint, int SystematicIndex, TString FastSimDataset, int Seed) {

  std::cout << FastSimDataset <<std::endl;

  //rand_ = new TRandom3(Seed);
  rand_ = new TRandom3(123);

  string CSVFileName = "/afs/cern.ch/user/b/brochero/brochero_WorkArea/CATTuples_Dic/TopCodeljets/TopTools/SF_btag/" + BTagAlgorithm + ".csv";
  BTagCalibration calib(BTagAlgorithm, CSVFileName);

  string SystematicFlagBC = "central";
  string SystematicFlagL  = "central";

  if (abs(SystematicIndex)<10) {
    if (SystematicIndex==-1 || SystematicIndex==-2) SystematicFlagBC = "down";
    if (SystematicIndex==+1 || SystematicIndex==+2) SystematicFlagBC = "up";
    if (SystematicIndex==-3) SystematicFlagL = "down";
    if (SystematicIndex==+3) SystematicFlagL = "up";
  }

  TaggerCut = -1;
  TaggerName = BTagAlgorithm;

  if (OperatingPoint=="Loose")  {
    if (TaggerName=="CSV") TaggerCut = 0.605;
    reader_bc = new BTagCalibrationReader(&calib, BTagEntry::OP_LOOSE, "comb", SystematicFlagBC);
    reader_l  = new BTagCalibrationReader(&calib, BTagEntry::OP_LOOSE, "comb", SystematicFlagL);
  } 
  else if (OperatingPoint=="Medium")  {
    if (TaggerName=="CSV") TaggerCut = 0.890;
    reader_bc = new BTagCalibrationReader(&calib, BTagEntry::OP_MEDIUM, "mujets", SystematicFlagBC);
    reader_l  = new BTagCalibrationReader(&calib, BTagEntry::OP_MEDIUM, "comb", SystematicFlagL);
  } 
  else if (OperatingPoint=="Tight")  {
    if (TaggerName=="CSV")  TaggerCut = 0.970;
    reader_bc = new BTagCalibrationReader(&calib, BTagEntry::OP_TIGHT, "mujets", SystematicFlagBC);
    reader_l  = new BTagCalibrationReader(&calib, BTagEntry::OP_TIGHT, "comb", SystematicFlagL);
  } 

  if (TaggerCut<0) 
    cout << " " << TaggerName << " not supported for " << OperatingPoint << " WP" << endl;


  if (TaggerCut<0.) 
    cout << "BTagSFUtil: " << BTagAlgorithm << " not a supported b-tagging algorithm" << endl;

}

BTagSFUtil::~BTagSFUtil() {

  delete rand_;

}

float BTagSFUtil::JetTagEfficiency(int JetFlavor, float JetPt, float JetEta) {

  if      (abs(JetFlavor)==5) return TagEfficiencyB(JetPt, JetEta);
  else if (abs(JetFlavor)==4) return TagEfficiencyC(JetPt, JetEta);
  else                        return TagEfficiencyLight(JetPt, JetEta);

}

float BTagSFUtil::GetJetSF(int JetFlavor, float JetPt, float JetEta) {

  float Btag_SF;
  float ThisJetPt = JetPt;
  if (abs(JetFlavor)==4 || abs(JetFlavor)==5) {
    if (JetPt>669.99) ThisJetPt = 669.99;
  }  else if (JetPt>999.99) ThisJetPt = 999.99;

  if (abs(JetFlavor)==5) 
    Btag_SF = reader_bc->eval(BTagEntry::FLAV_B, JetEta, ThisJetPt);
  else if (abs(JetFlavor)==4) 
    Btag_SF = reader_bc->eval(BTagEntry::FLAV_C, JetEta, ThisJetPt);
  else 
    Btag_SF = reader_l->eval(BTagEntry::FLAV_UDSG, JetEta, ThisJetPt);
    
  return Btag_SF;

}

bool BTagSFUtil::IsTagged(float JetDiscriminant, int JetFlavor, float JetPt, float JetEta) {
  bool isBTagged = JetDiscriminant > TaggerCut;
  if (JetFlavor==-999999) return isBTagged;// Data: no correction needed

  bool newBTag = isBTagged;
  float Btag_SF = GetJetSF(JetFlavor, JetPt, JetEta);
  if (Btag_SF == 1) return newBTag; //no correction needed 

  //throw dice
  float coin = rand_->Uniform(1.);    
 
  if(Btag_SF > 1){  // use this if SF>1
  

    if( !isBTagged ) {

      float Btag_eff = JetTagEfficiency(JetFlavor, JetPt, fabs(JetEta));

      //fraction of jets that need to be upgraded
      float mistagPercent = (1.0 - Btag_SF) / (1.0 - (1.0/Btag_eff) );
      
      //upgrade to tagged
      if( coin < mistagPercent ) {
	newBTag = true;
      }
    }
  }
  else{  // use this if SF<1
      
    //downgrade tagged to untagged
    if( isBTagged && coin > Btag_SF ) {
      newBTag = false;
    }

  }

  return newBTag;

}

float BTagSFUtil::TagEfficiencyB(float JetPt, float JetEta) {
  // CSVMedium b efficiency
  if (TaggerCut > 0.8 && TaggerCut < 0.9){
    if (JetPt > 20 && JetPt <= 40){
      if      (fabs(JetEta) > 0   && fabs(JetEta) <= 0.6) return 0.615984;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.60655;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.561642;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.468765;
    }
    else if (JetPt > 40 && JetPt <= 60){
      if      (fabs(JetEta) > 0   && fabs(JetEta) <= 0.6) return 0.66877;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.662889;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.614205;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.535718;
    }
    else if (JetPt > 60 && JetPt <= 80){
      if      (fabs(JetEta) > 0   && fabs(JetEta) <= 0.6) return 0.684861;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.687402;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.639179;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.566327;
    }
    else if (JetPt > 80 && JetPt <= 100){
      if      (fabs(JetEta) > 0   && fabs(JetEta) <= 0.6) return 0.698976;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.698125;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.653923;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.567678;
    }
    else if (JetPt > 100 && JetPt <= 120){
      if      (fabs(JetEta) > 0   && fabs(JetEta) <= 0.6) return 0.69307;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.694449;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.636424;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.545673;
    }
    else if (JetPt > 120 && JetPt <= 1000){
      if      (fabs(JetEta) > 0   && fabs(JetEta) <= 0.6) return 0.701633;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.700276;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.630122;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.555133;
    }
    else return 0.0;
  }

  // CSVTight b efficiency
  else if (TaggerCut > 0.9 && TaggerCut < 1.0){
    if (JetPt > 20 && JetPt <= 40){
      if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.379539;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.345202;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.344406;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.251914;
    }
    else if (JetPt > 40 && JetPt <= 60){
      if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.45233;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.439536;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.411741;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.335119;
    }
    else if (JetPt > 60 && JetPt <= 80){
      if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.470329;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.469562;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.432281;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.361218;
    }
    else if (JetPt > 80 && JetPt <= 100){
      if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.47489;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.475878;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.435468;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.357269;
    }
    else if (JetPt > 100 && JetPt <= 120){
      if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.463199;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.468148;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.403141;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.346495;
    }
    else if (JetPt > 120 && JetPt <= 140){
      if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.471426;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.467574;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.402546;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.339674;
    }
    else return 0.0;
  }
}

float BTagSFUtil::TagEfficiencyC(float JetPt, float JetEta) {
  // CSVMedium c efficiency
  if (TaggerCut > 0.8 && TaggerCut < 0.9){
    if (JetPt > 20 && JetPt <= 40){
      if      (fabs(JetEta) > 0   && fabs(JetEta) <= 0.6) return 0.140064;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.128012;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.13083;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.0988929;
    }
    else if (JetPt > 40 && JetPt <= 60){
      if      (fabs(JetEta) > 0   && fabs(JetEta) <= 0.6) return 0.131538;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.126563;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.129976;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.106902;
    }
    else if (JetPt > 60 && JetPt <= 80){
      if      (fabs(JetEta) > 0   && fabs(JetEta) <= 0.6) return 0.139282;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.129809;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.128481;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.0969443;
    }
    else if (JetPt > 80 && JetPt <= 100){
      if      (fabs(JetEta) > 0   && fabs(JetEta) <= 0.6) return 0.142206;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.13893;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.140348;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.110479;
    }
    else if (JetPt > 100 && JetPt <= 120){
      if      (fabs(JetEta) > 0   && fabs(JetEta) <= 0.6) return 0.134041;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.141719;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.140077;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.119207;
    }
    else if (JetPt > 120 && JetPt <= 1000){
      if      (fabs(JetEta) > 0   && fabs(JetEta) <= 0.6) return 0.14852;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.145738;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.140673;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.0899382;
    }
    else return 0.0;
  }

  // CSVTight c efficiency
  if (TaggerCut > 0.9 && TaggerCut < 1.0){
    if (JetPt > 20 && JetPt <= 40){
      if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.0186641;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.0150869;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.0205993;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.0112705;
    }
    else if (JetPt > 40 && JetPt <= 60){
      if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.0235886;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.0218171;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.0307017;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.0204706;
    }
    else if (JetPt > 60 && JetPt <= 80){
      if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.0254829;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.0227086;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.030127;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.0225034;
    }
    else if (JetPt > 80 && JetPt <= 100){
      if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.0260879;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.0204292;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.0334066;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.023649;
    }
    else if (JetPt > 100 && JetPt <= 120){
      if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.0216412;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.0292853;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.0319888;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.0176128;
    }
    else if (JetPt > 120 && JetPt <= 140){
      if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.0307386;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.0321199;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.0291277;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.013191;
    }  
    else return 0.0;
  }
}

float BTagSFUtil::TagEfficiencyLight(float JetPt, float JetEta) {
  // CSVMedium l efficiency
  if (TaggerCut > 0.8 && TaggerCut < 0.9){
    if (JetPt > 20 && JetPt <= 40){
      if      (fabs(JetEta) > 0   && fabs(JetEta) <= 0.6) return 0.00842663;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.00912602;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.0102513;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.0100411;
    }
    else if (JetPt > 40 && JetPt <= 60){
      if      (fabs(JetEta) > 0   && fabs(JetEta) <= 0.6) return 0.00855914;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.00842964;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.00969896;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.00958788;
    }
    else if (JetPt > 60 && JetPt <= 80){
      if      (fabs(JetEta) > 0   && fabs(JetEta) <= 0.6) return 0.00841776;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.00798238;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.0101725;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.0083111;
    }
    else if (JetPt > 80 && JetPt <= 100){
      if      (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.0110561;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.0121583;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.0115468;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.0121956;
    }
    else if (JetPt > 100 && JetPt <= 120){
      if      (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.0132316;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.0132753;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.0131566;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.0113801;
    }
    else if (JetPt > 120 && JetPt <= 1000){
      if      (fabs(JetEta) > 0   && fabs(JetEta) <= 0.6) return 0.016496;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.0139797;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.00908977;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.00984587;
    }
    else return 0.0;
  }
  
  // CSVTight l efficiency
  if (TaggerCut > 0.9 && TaggerCut < 1.0){
    if (JetPt > 20 && JetPt <= 40){
      if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.00141487;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.0011681;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.00116283;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.000603457;
    }
    else if (JetPt > 40 && JetPt <= 60){
      if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.00231471;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.00223418;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.00147519;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.00187843;
    }
    else if (JetPt > 60 && JetPt <= 80){
      if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.00392123;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.00326615;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.00283517;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.00165348;
    }
    else if (JetPt > 80 && JetPt <= 100){
      if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.00422838;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.00536069;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.00402618;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.00446061;
    }
    else if (JetPt > 100 && JetPt <= 120){
      if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.00518832;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.00491131;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.00569144;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.00384262;
    }
    else if (JetPt > 120 && JetPt <= 140){
      if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.00776751;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.00515171;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.00419529;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.000508386;
    }
    else return 0.0;
  }
}
