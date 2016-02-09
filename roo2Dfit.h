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


typedef struct HistoFit{
  TH1F     *hist1D[2][3];
  TH2F     *hist2D[3];
  float     events[3];
} HistoFit;


TString dirnameIn= "TopResults/";
// TString fl  = "hSF-Lumi2260-v2_Tree_LepJets_v7-6-1_NewSF_btag_Spring15-bx25_2170pb-1";  
TString fl  = "hSF-paJet25GeV-v0_Tree_LepJets_v7-6-1_systunc_Spring15-bx25_2170pb-1";  

enum FileSample{data, 
		ttbb, ttb, ttcc, ttLF, ttccLF, Bkgtt, ttjj, 
		WJets, ZJets, SingleTop, VV, QCD, BkgOther, BkgFull};

HistoFit LoadSample(TString FileName);
