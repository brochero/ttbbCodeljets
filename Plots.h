#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TObject.h"
#include "TStopwatch.h"
#include "TLegend.h"
#include <vector>
#include <string>
#include <fstream>
#include <vector>
#include <sys/stat.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <set>
#include "tdrstyle.C"
 

typedef struct histos{
  TH1F     *hist[3]; 
  THStack  *mc[3];
} histos;


/***********************
      Functions
************************/
std::vector<histos> loadhistograms (TString plots, TString namefile);
std::vector<histos> addhistograms  (std::vector<histos> histos_0, std::vector<histos> histos_1);
void setuphistograms               (std::vector<histos> histos, int color);
std::vector<histos> addstack       (std::vector<histos> stack_0, std::vector<histos> histos_0);

void overwritehistograms (std::vector<histos> newhistos, TString plots, TString namefile);

/***********************
 Files and Directories
************************/
TString dirnameIn= "TopResults/";
// TString fl  = "hSF-PrApr_v0_Tree_LepJets_v7-6-2_v0_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-CSVv2T-v0_Tree_LepJets_v7-6-2_v0_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-CSVv2M-v0_Tree_LepJets_v7-6-2_v0_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-TestJetMass_Tree_LepJets_v7-6-2_v0_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-FullPhSpCat_Tree_LepJets_v7-6-3_AddcJets_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-FullPhSpCat_Tree_LepJets_v7-6-3_AddcJets_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-minDeltaInvMassjj_Tree_LepJets_v7-6-2_v0_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-ANUpdate_Tree_LepJets_v7-6-3_ANsu_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-ANUpdateClassic_Tree_LepJets_v7-6-3_ANsu_Spring15-bx25_2260pb-1";  
// TString fl  = "hAcc-PureWjj_Tree_LepJets_v7-6-3_WGENJets_Spring15-bx25_2260pb-1";  
// TString fl  = "hAcc-Purettjj_With_Wjj_Tree_LepJets_v7-6-3_WGENJets_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-CSVPlots-NoW_Tree_LepJets_v7-6-3_WMatched_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-PrintTest_Tree_LepJets_v7-6-3_WMatched_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-WjjEnd_Tree_LepJets_v7-6-3_WMatched_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-NJets6_MVAWithTop_Tree_LepJets_v7-6-4_DRAddJets_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-NJets6_MVAAllJets_Tree_LepJets_v7-6-4_DRAddJets_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-NoMVAOrder_Tree_LepJets_v7-6-4_ctag_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-NoMVAOrder_Tree_LepJets_v7-6-4_ctag_Spring15-bx25_2260pb-1";  

// TString fl  = "hSF-MVASce1_Tree_LepJets_v7-6-4_ctag_Spring15-bx25_2260pb-1";  
TString fl  = "hSF-NoMVASce_Tree_LepJets_v7-6-4_ctag_Spring15-bx25_2260pb-1";  
