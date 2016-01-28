/***************************************************************************************************************************************
   ttbar Categorization from:
   https://twiki.cern.ch/twiki/pub/CMSPublic/GenHFHadronMatcher/GenTtbarCategorizer.cc


     The classification scheme returns an ID per event, and works as follows:
     
     All jets in the following need to be in the acceptance as given by the config parameters |eta|, pt.
     
     First, b jets from top are identified, i.e. jets containing a b hadron from t->b decay
     They are encoded in the ID as numberOfBjetsFromTop*100, i.e.
     0xx: no b jets from top in acceptance
     1xx: 1 b jet from top in acceptance
     2xx: both b jets from top in acceptance
     
     Then, b jets from W are identified, i.e. jets containing a b hadron from W->b decay
     They are encoded in the ID as numberOfBjetsFromW*1000, i.e.
     0xxx: no b jets from W in acceptance
     1xxx: 1 b jet from W in acceptance
     2xxx: 2 b jets from W in acceptance
     
     Then, c jets from W are identified, i.e. jets containing a c hadron from W->c decay, but no b hadrons
     They are encoded in the ID as numberOfCjetsFromW*10000, i.e.
     0xxxx: no c jets from W in acceptance
     1xxxx: 1 c jet from W in acceptance
     2xxxx: 2 c jets from W in acceptance
     
     From the remaining jets, the ID is formed based on the additional b jets (IDs x5x) and c jets (IDs x4x) in the following order:
     x55: at least 2 additional b jets with two of them having >= 2 b hadrons
     x54: at least 2 additional b jets with one of them having >= 2 b hadrons, the other having =1 b hadron
     x53: at least 2 additional b jets with all having =1 b hadron
     x52: exactly 1 additional b jet having >=2 b hadrons
     x51: exactly 1 additional b jet having =1 b hadron
     x45: at least 2 additional c jets with two of them having >= 2 c hadrons
     x44: at least 2 additional c jets with one of them having >= 2 c hadrons, the other having =1 c hadron
     x43: at least 2 additional c jets with all having =1 c hadron
     x42: exactly 1 additional c jet having >=2 c hadrons
     x41: exactly 1 additional c jet having =1 c hadron
     x00: No additional b or c jet, i.e. only light flavour jets

***************************************************************************************************************************************/

#include "TString.h"

bool ttbar_category(TString ttbarCatName, int ttbarGenCatId, int NGenJets){

  bool tt_Evt   = false;
  bool ttb_Evt  = false;
  bool ttc_Evt  = false;
  bool ttjj_Evt = false;
  bool ttbb_Evt = false;
  bool ttcc_Evt = false;

  int h, tu, t, u;

  h  = ttbarGenCatId / 100;
  tu = ttbarGenCatId % 100;
  t  = tu / 10;
  u  = tu % 10;
    
  if (tu == 0){
    // // TEMPORAL (Inside Acceptance (pT>25,|eta|<2.4))
    // if ((h == 2 || h==102) && NGenJets==6)      ttjj_Evt = true;
    // else if ((h == 1 || h==101) && NGenJets==5) ttjj_Evt = true;
    // else if ((h == 0 || h==100) && NGenJets==4) ttjj_Evt = true;
    // else tt_Evt   = true;
    tt_Evt   = true;
  }

  else if (tu == 55 || tu == 54 || tu == 53) ttbb_Evt = true;
  else if (tu == 52 || tu == 51)             ttb_Evt  = true;
  else if (tu == 45 || tu == 44 || tu == 43) ttcc_Evt = true;
  else if (tu == 42 || tu == 41)             ttc_Evt  = true;
  
  if (ttbarCatName == "tt")   return tt_Evt;
  if (ttbarCatName == "ttc")  return ttc_Evt;
  if (ttbarCatName == "ttb")  return ttb_Evt;
  if (ttbarCatName == "ttjj") return ttjj_Evt;
  if (ttbarCatName == "ttcc") return ttcc_Evt;
  if (ttbarCatName == "ttbb") return ttbb_Evt;    
  
}

//bool ttbar_conecategory(TString ttbarCatName, int ttbarGenCatId){
  // [0]: Decay mode. "semiLeptonic(0)" Just as a Cross Check!
  // [1]: Number of Jets "NJets20()"
  // [2]: Number of b-Jets "NbJets20()"
  // [3]: Number of c-Jets "NcJets20()"
  // [4]: Number of b-Jets Not comming from the top "NbJets20NoTop()"
  // [5]: Number of add Jets "NaddJets20()"
  // [6]: Number of add b-Jets "NaddbJets20()"



//}
