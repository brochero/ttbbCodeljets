#ifndef __CINT__

#include<string>
#include<iostream>
#include<sstream>
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<set>
#include<vector>

#endif

// Root
#include "TString.h"
#include "TRegexp.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TTree.h"
#include <sys/stat.h>
#include "TError.h"

// TopCode
  
#ifndef __CINT__

void display_usage()
{
  std::cout << "\033[1;37musage:\033[1;m skimfile cutindex [options]" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "    -i   inputfile  Input file without .root" << std::endl;
  std::cout << "    -o   name in the output file \"h_\"" << std::endl;
  std::cout << "    -cat ttbar categorization" << std::endl;
  std::cout << "    -d   Input file directory. Default directory: InputTrees" << std::endl;
  std::cout << "    -s create a file with the systematic uncertainty yields" << std::endl;
  std::cout << "    -h   displays this help message and exits " << std::endl;
  std::cout << "" << std::endl;
}


// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const TString currentDateTime() {

  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[80];
  tstruct = *localtime(&now);
  // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
  // for more information about date/time format
  strftime(buf, sizeof(buf), "%Y-%m-%d at %X", &tstruct);

  return buf;
}

void print_progress(int TreeEntries, Long64_t ievt)
{
  int step = TreeEntries/50;
  if (ievt%(step) == 0){ 
    float progress=(ievt)/(TreeEntries*1.0);
    int barWidth = 50;
    
    std::cout << "[";
    int pos = barWidth * progress;
    
    for (int i = 0; i < barWidth; ++i) {
      if (i < pos) std::cout << "=";
      else if (i == pos) std::cout << ">";
      else std::cout << " ";
    }
    
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
  }  
}

int main(int argc, const char* argv[]){

  gErrorIgnoreLevel = kError;

  gSystem->Load("libTree");

  bool   _ttbar_cat = false;
  bool   _syst      = false;
  const char * _output   = 0;
  const char * _input    = 0;
  // TopTrees directory
  const char * _dir      = "../Files_v7-6-3/";
  const char * _syst_var = 0;
  const char * _ttbar_id = 0;

  // Arguments used
  //std::set<int> usedargs;
  //Parsing input options
  if(argc == 1){
    display_usage();
    return -1;
  }

  else{
      //Argumet 1 must be a valid input fileName
      for (int i = 1; i < argc; i++){
	if( strcmp(argv[i],"-i") == 0 ){
	  _input = argv[i+1];
	  i++;
	}
	if( strcmp(argv[i],"-d") == 0 ){
	  _dir = argv[i+1];
	  i++;
	}
	if( strcmp(argv[i],"-o") == 0 ){
	  _output= argv[i+1];
	  i++;
	}
	if( strcmp(argv[i],"-cat") == 0 ){
	  _ttbar_cat = true;
	  _ttbar_id  = argv[i+1];
	  i++;
	}
	if( strcmp(argv[i],"-s") == 0 ){
	  _syst= true;
	  _syst_var = argv[i+1];
	}
	if( strcmp(argv[i],"-h") == 0 ||
	    strcmp(argv[i],"--help") == 0 ){
	  display_usage();
	  return 0;
	}
      }
  }//else
  if( _input ==0 ){
    std::cerr << "\033[1;31mskimfile ERROR:\033[1;m The '-i' option is mandatory!"
	      << std::endl;
    display_usage();
    return -1;
  }
  
  // reassigning
  TString fname(_input);
  TString hname(_output);
  TString fdir(_dir);
  TString ttbar_id(_ttbar_id);
  TString syst_varname(_syst_var);
  
  TChain theTree("ttbbLepJets/gentree"); 
  
  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  std::cout << "Signal: ";
  std::cout << fname + ".root" << std::endl;

  theTree.Add(fdir + fname + ".root");

  int Channel;
  float GENWeight; 
  std::vector<float> *ScaleWeight=0;
  float Lep_pT, Lep_eta;
  std::vector<float> *Jet_pT=0;
  std::vector<int>   *Jet_partonFlavour=0;

  // Categorization
  int  GenCat_ID;
  std::vector<int> *GenConeCat=0; 
 /*********************************
           Tree Branches
  **********************************/
  
  theTree.SetBranchAddress( "genweight", &GENWeight );
  theTree.SetBranchAddress( "scaleweight",  &ScaleWeight);
  theTree.SetBranchAddress( "genchannel",   &Channel);
  theTree.SetBranchAddress( "gencatid",     &GenCat_ID);
  theTree.SetBranchAddress( "genconecatid", &GenConeCat);


  theTree.SetBranchAddress( "genlepton_pT",  &Lep_pT);
  theTree.SetBranchAddress( "genlepton_eta", &Lep_eta);
  theTree.SetBranchAddress( "genjet_pT",     &Jet_pT);

    
  /*********************************
             Histograms
  **********************************/

  //Correct Statistical Uncertainty Treatment
  TH1::SetDefaultSumw2(kTRUE);  
  
  TH1F *hNJets[2], *hNBtagJets[2];

  TString namech[2];
  namech[0]="mujets";
  namech[1]="ejets";
  
  TString titlenamech[2];
  titlenamech[0]="#mu+Jets";
  titlenamech[1]="e+Jets";
  
  for(int i=0; i<2; i++){ // Channel
    hNJets[i]      = new TH1F("hNJets_" + namech[i], "Jet multiplicity " + titlenamech[i],9,-0.5,8.5);
    hNJets[i]->GetXaxis()->SetTitle("Number of jets");      
    hNJets[i]->GetXaxis()->SetBinLabel(1,"0");
    hNJets[i]->GetXaxis()->SetBinLabel(2,"1");
    hNJets[i]->GetXaxis()->SetBinLabel(3,"2");
    hNJets[i]->GetXaxis()->SetBinLabel(4,"3");
    hNJets[i]->GetXaxis()->SetBinLabel(5,"4");
    hNJets[i]->GetXaxis()->SetBinLabel(6,"5");
    hNJets[i]->GetXaxis()->SetBinLabel(7,"6");
    hNJets[i]->GetXaxis()->SetBinLabel(8,"7");
    hNJets[i]->GetXaxis()->SetBinLabel(9,"#geq 8");
    
    hNBtagJets[i]  = new TH1F("hNBtagJets_"+namech[i], "b-tag jet multiplicity " + titlenamech[i],6,-0.5,5.5);
    hNBtagJets[i]->GetXaxis()->SetTitle("Number of b-jets");
    hNBtagJets[i]->GetXaxis()->SetBinLabel(1,"0");
    hNBtagJets[i]->GetXaxis()->SetBinLabel(2,"1");
    hNBtagJets[i]->GetXaxis()->SetBinLabel(3,"2");
    hNBtagJets[i]->GetXaxis()->SetBinLabel(4,"3");
    hNBtagJets[i]->GetXaxis()->SetBinLabel(5,"4");
    hNBtagJets[i]->GetXaxis()->SetBinLabel(6,"#geq 5");
  }//for(i)
  
  
  TStopwatch sw;
  sw.Start(kTRUE);

  ///////////////////////////////////////
  // Please, IGNORE. Temporal solution //
  ///////////////////////////////////////
  TCanvas *mydummycanvas=new TCanvas();// 
  ///////////////////////////////////////
  // Number de events for acceptance
  //          [Channel]
  float fAccEvent_full_ttbb[2]={0.0,0.0};
  float fAccEvent_full_ttjj[2]={0.0,0.0};
  
  float fAccEvent_vis_ttbb[2]={0.0,0.0};
  float fAccEvent_vis_ttjj[2]={0.0,0.0};

  int AccEvent_full_ttbb[2]={0,0};
  int AccEvent_full_ttjj[2]={0,0};
  
  int AccEvent_vis_ttbb[2]={0,0};
  int AccEvent_vis_ttjj[2]={0,0};

  // Uncertainties file name
  if(_syst) fname += "_SYS_" + syst_varname;

  /********************************
             Event Loop
  ********************************/
  std::cout << "--- Processing: " << theTree.GetEntries() << " events" << std::endl;
  
  for (Long64_t ievt=0; ievt<theTree.GetEntries();ievt++) {
    
    theTree.GetEntry(ievt);  
    print_progress(theTree.GetEntries(), ievt);
    
    int NJets, NBtagJets;

    NBtagJets = 0;    
    // Jets 
    NJets = 0;
    for(int ijet=0; ijet < Jet_pT->size(); ijet++){
      
      float jet_pT = (*Jet_pT)[ijet];
      
      if(jet_pT>20) NJets++; // Number of GEN-Jets
      
    }// for(jets)
    
    
    /***************************
        Categorization 1
    ***************************/
    // int h, tu, t, u;
    
    // bool tt_Evt   = false;
    // bool ttb_Evt  = false;
    // bool ttc_Evt  = false;
    // bool ttjj_Evt = false;
    // bool ttbb_Evt = false;
    // bool ttcc_Evt = false;
    
    // h  = GenCat_ID / 100;
    // tu = GenCat_ID % 100;
    // t  = tu / 10;
    // u  = tu % 10;

    // if (tu == 0) tt_Evt   = true;
    // else if (tu == 55 || tu == 54 || tu == 53) ttbb_Evt = true;
    // else if (tu == 52 || tu == 51)             ttb_Evt  = true;
    // else if (tu == 45 || tu == 44 || tu == 43) ttcc_Evt = true;
    // else if (tu == 42 || tu == 41)             ttc_Evt  = true;
    
    // NBtagJets  = 0;
    
    // if (h == 2 || h==102)      NBtagJets += 2;
    // else if (h == 1 || h==101) NBtagJets += 1;
    // else if (h == 0 || h==100) NBtagJets  = 0;
    
    
    // if (!ttbb_Evt) continue;
    // if (ttbb_Evt) NBtagJets += 2;
    
    /***************************
        Categorization 2
    ***************************/

    int cone_channel = (*GenConeCat)[0];

    // Visible Phase Space:
    // pT(jet) > 20GeV && |eta(Jet)| < 2.5
    int cone_NJets  = (*GenConeCat)[1];
    int cone_NbJets = (*GenConeCat)[2];
    int cone_NcJets = (*GenConeCat)[3];
    
    int cone_NbJetsNoTop = (*GenConeCat)[4];
    
    // Full Phase Space:
    // pT(jet) > 20GeV && |eta(Jet)| < 2.5
    int cone_NaddJets  = (*GenConeCat)[5];
    int cone_NaddbJets = (*GenConeCat)[6];


    /******************
      Scale Weights
    ******************/
    int scaleSysPar;
    if     (_syst && syst_varname.Contains("ScaleRnF_Up"))   scaleSysPar = 0; // muR=Nom,  muF=Up
    else if(_syst && syst_varname.Contains("ScaleRnF_Down")) scaleSysPar = 1; // muR=Nom,  muF=Down
    else if(_syst && syst_varname.Contains("ScaleRuF_Nom"))  scaleSysPar = 2; // muR=Up,   muF=Nom
    else if(_syst && syst_varname.Contains("ScaleRuF_Up"))   scaleSysPar = 3; // muR=Up,   muF=Up
    else if(_syst && syst_varname.Contains("ScaleRdF_Nom"))  scaleSysPar = 4; // muR=Down, muF=Nom
    else if(_syst && syst_varname.Contains("ScaleRdF_Down")) scaleSysPar = 5; // muR=Down, muF=Down

    float EvtStep = GENWeight;
    
    if (_syst && syst_varname.Contains("ScaleR"))
      EvtStep = EvtStep*(*ScaleWeight)[scaleSysPar];

     /******************
         Acceptace
    ******************/
    if(cone_NaddJets  > 1) fAccEvent_full_ttjj[Channel]+=EvtStep;
    if(cone_NaddbJets > 1) fAccEvent_full_ttbb[Channel]+=EvtStep;

    if(Lep_pT > 30 && abs(Lep_eta) < 2.4){
      if(cone_NbJets > 1 && cone_NJets > 5) fAccEvent_vis_ttjj[Channel]+=EvtStep;
      if(cone_NbJets > 3 && cone_NJets > 5) fAccEvent_vis_ttbb[Channel]+=EvtStep;    
    }
    

    /******************
        Histograms
    ******************/
    hNJets[Channel]->Fill(NJets); 
  
  }//for(events)

  AccEvent_full_ttjj[0] = fAccEvent_full_ttjj[0];
  AccEvent_full_ttbb[0] = fAccEvent_full_ttbb[0];

  AccEvent_vis_ttjj[0] = fAccEvent_vis_ttjj[0];
  AccEvent_vis_ttbb[0] = fAccEvent_vis_ttbb[0];

  AccEvent_full_ttjj[1] = fAccEvent_full_ttjj[1];
  AccEvent_full_ttbb[1] = fAccEvent_full_ttbb[1];

  AccEvent_vis_ttjj[1] = fAccEvent_vis_ttjj[1];
  AccEvent_vis_ttbb[1] = fAccEvent_vis_ttbb[1];
  
  // Get elapsed time
  sw.Stop();
  std::cout << "==================================================] 100% " << std::endl;
  std::cout << "--- End of event loop: "; sw.Print();
  
  
  //Acceptance-Efficiency
  std::cout << "--------  Acceptace Full Ph-Sp  --------" << std::endl;
  std::cout << "Number of RAW-mu+Jets events:" << std::endl;
  std::cout << "ttjj Acceptance Full Ph-Sp: " << AccEvent_full_ttjj[0] << std::endl;
  std::cout << "ttbb Acceptance Full Ph-Sp: " << AccEvent_full_ttbb[0] << std::endl;
  std::cout << "ttbb/ttjj Full Ph-Sp: " << 1.0*AccEvent_full_ttbb[0]/AccEvent_full_ttjj[0] << std::endl;
  std::cout << "-----------------------------" << std::endl;
  std::cout << "Number of RAW-e+Jets events:" << std::endl;
  std::cout << "ttjj Acceptance Full Ph-Sp: " << AccEvent_full_ttjj[1] << std::endl;
  std::cout << "ttbb Acceptance: Full Ph-Sp: " << AccEvent_full_ttbb[1] << std::endl;
  std::cout << "ttbb/ttjj Full Ph-Sp: " << 1.0*AccEvent_full_ttbb[1]/AccEvent_full_ttjj[1] << std::endl;
  std::cout << "-----------------------------" << std::endl;
  std::cout << "Number of RAW-l+Jets events:" << std::endl;
  std::cout << "ttjj Acceptance Full Ph-Sp: " << AccEvent_full_ttjj[0] + AccEvent_full_ttjj[1] << std::endl;
  std::cout << "ttbb Acceptance Full Ph-Sp: " << AccEvent_full_ttbb[0] + AccEvent_full_ttbb[1] << std::endl;
  std::cout << "ttbb/ttjj Full Ph-Sp: " << 1.0*(AccEvent_full_ttbb[0] + AccEvent_full_ttbb[1])/(AccEvent_full_ttjj[0] + AccEvent_full_ttjj[1]) << std::endl;
  std::cout << "-----------------------------" << std::endl;
  std::cout << "-----------------------------" << std::endl;
  std::cout << "-----------------------------" << std::endl;
  std::cout << "--------  Acceptace Visible Ph-Sp  --------" << std::endl;
  std::cout << "Number of RAW-mu+Jets events:" << std::endl;
  std::cout << "ttjj Acceptance Visible Ph-Sp: " << AccEvent_vis_ttjj[0] << std::endl;
  std::cout << "ttbb Acceptance Visible Ph-Sp: " << AccEvent_vis_ttbb[0] << std::endl;
  std::cout << "ttbb/ttjj Visible Ph-Sp: " << 1.0*AccEvent_vis_ttbb[0]/AccEvent_vis_ttjj[0] << std::endl;
  std::cout << "-----------------------------" << std::endl;
  std::cout << "Number of RAW-e+Jets events:" << std::endl;
  std::cout << "ttjj Acceptance Visible Ph-Sp: " << AccEvent_vis_ttjj[1] << std::endl;
  std::cout << "ttbb Acceptance: Visible Ph-Sp: " << AccEvent_vis_ttbb[1] << std::endl;
  std::cout << "ttbb/ttjj Visible Ph-Sp: " << 1.0*AccEvent_vis_ttbb[1]/AccEvent_vis_ttjj[1] << std::endl;
  std::cout << "-----------------------------" << std::endl;
  std::cout << "Number of RAW-l+Jets events:" << std::endl;
  std::cout << "ttjj Acceptance Visible Ph-Sp: " << AccEvent_vis_ttjj[0] + AccEvent_vis_ttjj[1] << std::endl;
  std::cout << "ttbb Acceptance Visible Ph-Sp: " << AccEvent_vis_ttbb[0] + AccEvent_vis_ttbb[1] << std::endl;
  std::cout << "ttbb/ttjj Visible Ph-Sp: " << 1.0*(AccEvent_vis_ttbb[0] + AccEvent_vis_ttbb[1])/(AccEvent_vis_ttjj[0] + AccEvent_vis_ttjj[1]) << std::endl;
  std::cout << "-----------------------------" << std::endl;
  std::cout << "-----------------------------" << std::endl;
  std::cout << "-----------------------------" << std::endl;
  std::cout << "--------  ttbb Acceptace  --------" << std::endl;
  std::cout << "ttbb Acceptance (mu+Jets): " << 1.0*AccEvent_vis_ttbb[0]/AccEvent_full_ttbb[0] << std::endl;
  std::cout << "ttbb Acceptance (e+Jets): "  << 1.0*AccEvent_vis_ttbb[1]/AccEvent_full_ttbb[1] << std::endl;
  std::cout << "ttbb Acceptance (l+Jets): "  << 1.0*(AccEvent_vis_ttbb[0] + AccEvent_vis_ttbb[1])/(AccEvent_full_ttbb[0] + AccEvent_full_ttbb[1]) << std::endl;
  std::cout << "--------  ttjj Acceptace  --------" << std::endl;
  std::cout << "ttjj Acceptance (mu+Jets): " << 1.0*AccEvent_vis_ttjj[0]/AccEvent_full_ttjj[0] << std::endl;
  std::cout << "ttjj Acceptance (e+Jets): "  << 1.0*AccEvent_vis_ttjj[1]/AccEvent_full_ttjj[1] << std::endl;
  std::cout << "ttjj Acceptance (l+Jets): "  << 1.0*(AccEvent_vis_ttjj[0] + AccEvent_vis_ttjj[1])/(AccEvent_full_ttjj[0] + AccEvent_full_ttjj[1]) << std::endl;

  //Output Dir
  TString dirname="TopResults";   
  // make a dir if it does not exist!!
  struct stat st;
  if(stat(dirname,&st) != 0) system("mkdir " + dirname);
  

  // Sample name identification
  TString samplename="";
  bool matchsamplename=false;
  
  for(int i=0; i<fname.Sizeof(); i++){
    if (i>2){
      if (fname[i-3]=='-' && 
	  fname[i-2]=='1' && 
	  fname[i-1]=='_') matchsamplename=true;
    }
    if (matchsamplename) samplename.Append(fname[i]);
  }
    
  // --- Write histograms
  
  TString outfname=dirname + "/hAcc-" + hname + "_" + fname + ".root";
  TFile *target  = new TFile(outfname,"RECREATE" );  
  
  for(int i=0; i<2; i++){    
    hNJets[i]->Write();
    hNBtagJets[i]->Write();    
  }//for(i)
  
  std::cout << "File saved as " << outfname << std::endl;
  
}


#endif

