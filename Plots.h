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
//TString fl  = "hSF-Lumi2260-v2_Tree_LepJets_v7-6-1_NewSF_btag_Spring15-bx25_2170pb-1";  
TString fl  = "hSF-paJet25GeV-v0_Tree_LepJets_v7-6-1_systunc_Spring15-bx25_2170pb-1";  

