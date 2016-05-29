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

// Plot Style
#include "tdrstyle.C"

struct HistoFit{ 
  TH1F     *hist1D[2][3]; 
  TH2F     *hist2D[3]; 
  float    events[3]; 
}; 


TString dirnameIn= "TopResults/";
// TString fl  = "hSF-PrApr_v0_Tree_LepJets_v7-6-2_v0_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-CSVv2T-v0_Tree_LepJets_v7-6-2_v0_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-CSVv2M-v0_Tree_LepJets_v7-6-2_v0_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-FullPhSpCat_Tree_LepJets_v7-6-3_AddcJets_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-ANUpdateClassic_Tree_LepJets_v7-6-3_ANsu_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-ANUpdate_Tree_LepJets_v7-6-3_ANsu_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-NewSF_Tree_LepJets_v7-6-3_ANsu_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-NewCSVBin_Tree_LepJets_v7-6-3_ANsu_Spring15-bx25_2260pb-1";  

// TString fl  = "hSF-NJets6_MVAAllJets_Tree_LepJets_v7-6-4_DRAddJets_Spring15-bx25_2260pb-1";  
// TString fl  = "hSF-Nj6-GENOrder_Tree_LepJets_v7-6-4_ctag_Spring15-bx25_2260pb-1";  
TString fl  = "hSF-NoMVASce_Tree_LepJets_v7-6-4_ctag_Spring15-bx25_2260pb-1";  

enum FileSample{data, 
		ttbb, ttb, ttcc, ttLF, ttccLF, Bkgtt, ttjj, 
		WJets, ZJets, SingleTop, VV, QCD, BkgOther, BkgFull};

HistoFit LoadSample(TString FileName);
