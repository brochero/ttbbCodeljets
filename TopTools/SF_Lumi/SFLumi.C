#include<vector>
#include<iostream>
#include "TString.h"


/***************************
    SF: Normalization
***************************/
      
float SFLumi(TString filename,
	     float Lumi,
	     float NGenEvents){
  
  /************************
    Normalization Weights
  *************************/

  float NormWeight = 0.0;
  // NormWeight = Lumi*(1.0/N_Gen_events)*(Xsec)*(Br)

  if(filename.Contains("QCD_MuEnr_20to30"))     NormWeight = Lumi * (1.0/NGenEvents) * (1273190000.0) * (0.003); // [pb] (cross section) * (Filter Eff)
  if(filename.Contains("QCD_MuEnr_20to30"))     NormWeight = Lumi * (1.0/NGenEvents) * (558528000.0) * (0.0053); // [pb] (cross section) * (Filter Eff)
  if(filename.Contains("QCD_MuEnr_30to50"))     NormWeight = Lumi * (1.0/NGenEvents) * (139803000.0) * (0.01182); // [pb] (cross section) * (Filter Eff)
  if(filename.Contains("QCD_MuEnr_50to80"))     NormWeight = Lumi * (1.0/NGenEvents) * (19222500.0) * (0.02276); // [pb] (cross section) * (Filter Eff)
  if(filename.Contains("QCD_MuEnr_80to120"))    NormWeight = Lumi * (1.0/NGenEvents) * (2758420.0) * (0.03844);  // [pb] (cross section) * (Filter Eff)
  if(filename.Contains("QCD_MuEnr_120to170"))   NormWeight = Lumi * (1.0/NGenEvents) * (469797.0) * (0.05362);   // [pb] (cross section) * (Filter Eff)
  if(filename.Contains("QCD_MuEnr_170to300"))   NormWeight = Lumi * (1.0/NGenEvents) * (117989.0) * (0.07335);   // [pb] (cross section) * (Filter Eff)
  if(filename.Contains("QCD_MuEnr_300to470"))   NormWeight = Lumi * (1.0/NGenEvents) * (7820.25) * (0.10196);    // [pb] (cross section) * (Filter Eff)
  if(filename.Contains("QCD_MuEnr_470to600"))   NormWeight = Lumi * (1.0/NGenEvents) * (645.528) * (0.12242);    // [pb] (cross section) * (Filter Eff)
  if(filename.Contains("QCD_MuEnr_800to1000"))  NormWeight = Lumi * (1.0/NGenEvents) * (32.3486) * (0.14552);    // [pb] (cross section) * (Filter Eff)
  if(filename.Contains("QCD_MuEnr_1000toInf"))  NormWeight = Lumi * (1.0/NGenEvents) * (10.4305) * (0.15544);    // [pb] (cross section) * (Filter Eff)

  if(filename.Contains("QCD_EGEnr_15to20"))     NormWeight = Lumi * (1.0/NGenEvents) * (1279000000.0) * (0.0018); // [pb] (cross section) * (Filter Eff)
  if(filename.Contains("QCD_EGEnr_20to30"))     NormWeight = Lumi * (1.0/NGenEvents) * (557600000.0) * (0.0096); // [pb] (cross section) * (Filter Eff)
  if(filename.Contains("QCD_EGEnr_30to50"))     NormWeight = Lumi * (1.0/NGenEvents) * (136000000.0) * (0.073); // [pb] (cross section) * (Filter Eff)
  if(filename.Contains("QCD_EGEnr_50to80"))     NormWeight = Lumi * (1.0/NGenEvents) * (19800000.0) * (0.146); // [pb] (cross section) * (Filter Eff)
  if(filename.Contains("QCD_EGEnr_80to120"))    NormWeight = Lumi * (1.0/NGenEvents) * (2800000.0) * (0.125); // [pb] (cross section) * (Filter Eff)
  if(filename.Contains("QCD_EGEnr_120to170"))   NormWeight = Lumi * (1.0/NGenEvents) * (477000.0) * (0.132); // [pb] (cross section) * (Filter Eff)
  if(filename.Contains("QCD_EGEnr_170to300"))   NormWeight = Lumi * (1.0/NGenEvents) * (114000.0) * (0.165); // [pb] (cross section) * (Filter Eff)
  if(filename.Contains("QCD_EGEnr_300toInf"))   NormWeight = Lumi * (1.0/NGenEvents) * (9000.0) * (0.15); // [pb] (cross section) * (Filter Eff)


  if(filename.Contains("ZJets_M50"))      NormWeight = Lumi * (1.0/NGenEvents) * (6025.2);  // [pb]
  if(filename.Contains("ZJets_M10to50"))  NormWeight = Lumi * (1.0/NGenEvents) * (18610.0);// [pb]
  if(filename.Contains("WJets"))          NormWeight = Lumi * (1.0/NGenEvents) * (61526.7); // [pb]
  if(filename.Contains("tW"))             NormWeight = Lumi * (1.0/NGenEvents) * (35.6);    // [pb]
  if(filename.Contains("tbarW"))          NormWeight = Lumi * (1.0/NGenEvents) * (35.6);    // [pb]
  if(filename.Contains("t_tchannel"))     NormWeight = Lumi * (1.0/NGenEvents) * (44.33);  // [pb]
  if(filename.Contains("tbar_tchannel"))  NormWeight = Lumi * (1.0/NGenEvents) * (26.38);   // [pb]
  if(filename.Contains("WW"))             NormWeight = Lumi * (1.0/NGenEvents) * (110.8);   // [pb]
  if(filename.Contains("WZ"))             NormWeight = Lumi * (1.0/NGenEvents) * (47.13);    // [pb]
  if(filename.Contains("ZZ"))             NormWeight = Lumi * (1.0/NGenEvents) * (16.52);    // [pb]                                                                                            
  
  if(filename.Contains("ttbar_PowhegPythia")) NormWeight = Lumi * (1.0/NGenEvents) * (831.76); // [pb] Br = (leptonic) * Hadronic = (0.1086*3) * (0.67)
  if(filename.Contains("ttbar_MCatNLO"))      NormWeight = Lumi * (1.0/NGenEvents) * (831.76); // * (0.1086*3.0*3.0); // Br correction
  if(filename.Contains("ttbar_Madgraph"))     NormWeight = Lumi * (1.0/NGenEvents) * (831.76); 
  
  if(filename.Contains("Data"))               NormWeight = 1.0;
  
  return NormWeight;
}

