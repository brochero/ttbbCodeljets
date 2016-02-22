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
//#include "TKey.h"
//#include "TPrint.h"
//#include <exception>
#include <sys/stat.h>
#include "TError.h"

// TopCode
#include <SFIDISOTrigger.h> // SF_ID-ISO-Trigger
#include <ttbar_category.h> // Event Categorization
#include <BTagSFUtil.h>     // btag SF
#include <SFLumi.h>         // Normalization SF

#ifndef __CINT__

void display_usage()
{
  std::cout << "\033[1;37musage:\033[1;m skimfile cutindex [options]" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "    -i inputfile  Input file without .root" << std::endl;
  std::cout << "    -o name in the output file \"h_\"" << std::endl;
  std::cout << "    -s create a file with the systematic uncertainty yields" << std::endl;
  std::cout << "    -tr SF Trigger Uncertainty" << std::endl;
  std::cout << "    -idiso SF ID/ISO Uncertainty" << std::endl;
  std::cout << "    -cat ttbar categorization" << std::endl;
  std::cout << "    -d Input file directory. Default directory: InputTrees" << std::endl;
  std::cout << "    -h                 displays this help message and exits " << std::endl;
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


//enum class btagUnc:unsigned int{CENTRAL,
struct btagUnc{
  enum:unsigned int{CENTRAL,
      JES_UP,     JES_DN,
      LF_UP,       LF_DN,
      HF_UP,       HF_DN,
      HFSTAT1_UP,  HFSTAT1_DN, 
      HFSTAT2_UP,  HFSTAT2_DN,
      LFSTAT1_UP,  LFSTAT1_DN, 
      LFSTAT2_UP,  LFSTAT2_DN,
      CFERR1_UP,   CFERR1_DN, 
      CFERR2_UP,   CFERR2_DN};
};


int main(int argc, const char* argv[]){

  gErrorIgnoreLevel = kError;

  gSystem->Load("libTree");
  gROOT->ProcessLine("#include <vector>");

  bool   _ttbar_cat = false;
  bool   _syst      = false;
  bool	 _tr_unc    = false;
  bool	 _idiso_unc = false;
  const char * _output   = 0;
  const char * _input    = 0;
  // TopTrees directory
  const char * _dir      = "../Files_v7-6-2/";
  const char * _syst_var = 0;
  const char * _tr       = 0;
  const char * _idiso    = 0;
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
	if( strcmp(argv[i],"-s") == 0 ){
	  _syst= true;
	  _syst_var = argv[i+1];
	}
	if( strcmp(argv[i],"-tr") == 0 ){
	  _tr_unc= true;
	  _tr= argv[i+1];
	  i++;
	}
	if( strcmp(argv[i],"-idiso") == 0 ){
	  _idiso_unc= true;
	  _idiso= argv[i+1];
	  i++;
	}
	if( strcmp(argv[i],"-cat") == 0 ){
	  _ttbar_cat = true;
	  _ttbar_id  = argv[i+1];
	  i++;
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
  TString syst_varname(_syst_var);
  TString fdir(_dir);
  TString TrUnc(_tr);
  TString IDISOUnc(_idiso);
  TString ttbar_id(_ttbar_id);
  
  TChain theTree("ttbbLepJets/tree"); 
  
  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  std::cout << "Signal: ";
  std::cout << fname + ".root" << std::endl;

  theTree.Add(fdir + fname + ".root");

  if(!theTree.GetFile()){
    std::cerr << "Input file not found!!!"  << std::endl;
    std::exit(0);
  }


  int Event,Run,Channel, GoodPV;
  float PUWeight, GENWeight; 
  std::vector<float> *PUWeight_sys=0;


  float MET,MET_Phi;

  float Lep_px, Lep_py, Lep_pz, Lep_E;
  std::vector<float> *Lep_SF=0;

  std::vector<float> *Jet_px=0, *Jet_py=0, *Jet_pz=0, *Jet_E=0;
  std::vector<int>   *Jet_partonFlavour=0;
  std::vector<float> *Jet_CSV=0;
  std::vector<float> *Jet_SF_CSV=0;
  std::vector<float> *Jet_JER_Up=0, *Jet_JER_Nom=0, *Jet_JER_Down=0;
  std::vector<float> *Jet_JES_Up=0, *Jet_JES_Down=0;

  // GEN Info
  std::vector<int> *GenConeCat=0;
  float GenLep_pT;
  std::vector<float> *GenJet_pT=0;

  // Scale Syst. Unc.
  std::vector<float> *ScaleWeight=0;
   
 /*********************************
           Tree Branches
  **********************************/
  
  theTree.SetBranchAddress( "event",    &Event );
  theTree.SetBranchAddress( "run",      &Run );

  theTree.SetBranchAddress( "PUWeight",  &PUWeight_sys );
  theTree.SetBranchAddress( "channel",   &Channel );

  theTree.SetBranchAddress( "GoodPV",  &GoodPV );

  theTree.SetBranchAddress( "MET",     &MET );
  theTree.SetBranchAddress( "MET_phi", &MET_Phi );

  theTree.SetBranchAddress( "lepton_px", &Lep_px );
  theTree.SetBranchAddress( "lepton_py", &Lep_py );
  theTree.SetBranchAddress( "lepton_pz", &Lep_pz );
  theTree.SetBranchAddress( "lepton_E",  &Lep_E );

  theTree.SetBranchAddress( "lepton_SF",  &Lep_SF );

  theTree.SetBranchAddress( "jet_px", &Jet_px );
  theTree.SetBranchAddress( "jet_py", &Jet_py );
  theTree.SetBranchAddress( "jet_pz", &Jet_pz );
  theTree.SetBranchAddress( "jet_E",  &Jet_E );

  theTree.SetBranchAddress( "jet_CSV",  &Jet_CSV );
  theTree.SetBranchAddress( "jet_SF_CSV",  &Jet_SF_CSV );
  theTree.SetBranchAddress( "jet_partonFlavour",  &Jet_partonFlavour );


  if(!fname.Contains("Data")){
    theTree.SetBranchAddress( "jet_JES_Up",  &Jet_JES_Up );
    theTree.SetBranchAddress( "jet_JES_Down",  &Jet_JES_Down );
    
    theTree.SetBranchAddress( "jet_JER_Up",  &Jet_JER_Up );
    theTree.SetBranchAddress( "jet_JER_Nom",  &Jet_JER_Nom );
    theTree.SetBranchAddress( "jet_JER_Down",  &Jet_JER_Down );
  }
  
  // Number of Events and Weights (MC@NLO)
  TFile *fileEntries = NULL;
  fileEntries = TFile::Open(fdir + fname + ".root");
  TH1F *h_NumberEvt;
  h_NumberEvt = (TH1F*)fileEntries->Get("ttbbLepJets/EventInfo");

  float NTotal_Event, NTotal_Weight, nNorm_Event, NTotal_ScalemuRF_Weight[6];
  NTotal_Event  = h_NumberEvt->GetBinContent(1);
  NTotal_Weight = h_NumberEvt->GetBinContent(2);
  for (unsigned int ibin = 0; ibin< 6; ibin++) NTotal_ScalemuRF_Weight[ibin]= h_NumberEvt->GetBinContent(ibin + 9);


  /************************
      MCatNLO Weights
  *************************/
  if(fname.Contains("MCatNLO")){
    theTree.SetBranchAddress( "genweight", &GENWeight );
    nNorm_Event = NTotal_Weight;    
  }
  else{
    GENWeight = 1.0;
    nNorm_Event = NTotal_Event;
  }

  // ttbar event categorization
  if(fname.Contains("ttbar") && !fname.Contains("Bkg")){
    theTree.SetBranchAddress("scaleweight",  &ScaleWeight );
    theTree.SetBranchAddress("genconecatid", &GenConeCat);
    theTree.SetBranchAddress("genlepton_pT", &GenLep_pT);
    theTree.SetBranchAddress("genjet_pT",    &GenJet_pT);
  }

  /*********************************
             Histograms
  **********************************/

  //Correct Statistical Uncertainty Treatment
  TH1::SetDefaultSumw2(kTRUE);  
  
  TH1F *hPV[4][2];

  TH1F *hMET[4][2], *hMET_Phi[4][2];

  TH1F *hLepPt[4][2], *hLepEta[4][2], *hLepPhi[4][2];

  TH1F *hmT[4][2];

  TH1F *hNJets[4][2], *hHT[4][2], *hNBtagJets[4][2];
  TH1F *hCSV[4][4][2], *hJetPt[4][4][2], *hJetpTUncVar[4][4][2];
  TH2F *h2DCSV_23Jet[4][2];

  TH1F *hSFpT[4][2], *hSFpTError[4][2];
  TH1F *hSFIDISO[4][2], *hSFIDISOError[4][2];
  TH1F *hSFTrigger[4][2], *hSFTriggerError[4][2];
  TH2F *h2DSFbtag_Global[4][2];
  TH1F *hSFbtag_Global[4][2], *hSFbtag_Global_var[4][2];
  TH2F *h2DSFbtag_b[4][2], *h2DSFbtag_c[4][2], *h2DSFbtag_l[4][2], *h2DSFbtag_btag_b[4][2], *h2DSFbtag_btag_c[4][2], *h2DSFbtag_btag_l[4][2]; 

  TH1F *hEvtCatego[4][2];

  TString namech[2];
  namech[0]="mujets";
  namech[1]="ejets";
  
  TString namecut[4];
  namecut[0]="lepton";
  namecut[1]="6Jets";
  namecut[2]="2btag";
  namecut[3]="4btag";
  
  TString titlenamech[2];
  titlenamech[0]="#mu+Jets";
  titlenamech[1]="e+Jets";
  
  for(int j=0; j<4; j++){   // Cut
    for(int i=0; i<2; i++){ // Channel
      hPV[j][i]         = new TH1F("hPV_"+namech[i]+"_"+namecut[j],"PV Distribution  " + titlenamech[i],15,0,30);
      hPV[j][i]->GetXaxis()->SetTitle("PV");      
      hMET[j][i]        = new TH1F("hMET_"+ namech[i]+"_"+namecut[j],"#slash{E}_{T} " + titlenamech[i],10,0,200);
      hMET[j][i]->GetXaxis()->SetTitle("#slash{E}_{T}[GeV]");      
      hMET_Phi[j][i]    = new TH1F("hMET_Phi_"+ namech[i]+"_"+namecut[j],"#Phi_{#slash{E}_{T}} " + titlenamech[i],160,-4,4);
      hMET_Phi[j][i]->GetXaxis()->SetTitle("#Phi_{#slash{E}_{T}}[rad]");      

      hmT[j][i]    = new TH1F("hmT_"+ namech[i]+"_"+namecut[j],"transverse Mass Lepton/MET " + titlenamech[i],40,0,160);
      hmT[j][i]->GetXaxis()->SetTitle("m_{T}[GeV]");      
      
      hLepPt [j][i]    = new TH1F("hLepPt_"  + namech[i] + "_" + namecut[j], "Lepton p_{T} " + titlenamech[i],20,0.0,200.0);
      hLepPt[j][i]->GetXaxis()->SetTitle("Lepton p_{T}[GeV]");      
      hLepEta[j][i]    = new TH1F("hLepEta_" + namech[i] + "_" + namecut[j], "#eta_{Lep} " + titlenamech[i],12,0.0,2.2);
      hLepEta[j][i]->GetXaxis()->SetTitle("Lepton #eta");      
      hLepPhi[j][i]    = new TH1F("hLepPhi_" + namech[i] + "_" + namecut[j], "#phi_{Lep} " + titlenamech[i],16,0.0,3.2);
      hLepPhi[j][i]->GetXaxis()->SetTitle("lepton #Phi[rad]");      
      
      hNJets[j][i]      = new TH1F("hNJets_" + namech[i] + "_" + namecut[j], "Jet multiplicity " + titlenamech[i],9,-0.5,8.5);
      hNJets[j][i]->GetXaxis()->SetTitle("Number of jets");      
      hNJets[j][i]->GetXaxis()->SetBinLabel(1,"0");
      hNJets[j][i]->GetXaxis()->SetBinLabel(2,"1");
      hNJets[j][i]->GetXaxis()->SetBinLabel(3,"2");
      hNJets[j][i]->GetXaxis()->SetBinLabel(4,"3");
      hNJets[j][i]->GetXaxis()->SetBinLabel(5,"4");
      hNJets[j][i]->GetXaxis()->SetBinLabel(6,"5");
      hNJets[j][i]->GetXaxis()->SetBinLabel(7,"6");
      hNJets[j][i]->GetXaxis()->SetBinLabel(8,"7");
      hNJets[j][i]->GetXaxis()->SetBinLabel(9,"#geq 8");
      if(j>0) hNJets[j][i]->GetXaxis()->SetRange(7,9);

      hNBtagJets[j][i]  = new TH1F("hNBtagJets_"+namech[i]+"_"+namecut[j],"b-tag jet multiplicity " + titlenamech[i],6,-0.5,5.5);
      hNBtagJets[j][i]->GetXaxis()->SetTitle("Number of b-jets");      
      hNBtagJets[j][i]->GetXaxis()->SetBinLabel(1,"0");
      hNBtagJets[j][i]->GetXaxis()->SetBinLabel(2,"1");
      hNBtagJets[j][i]->GetXaxis()->SetBinLabel(3,"2");
      hNBtagJets[j][i]->GetXaxis()->SetBinLabel(4,"3");
      hNBtagJets[j][i]->GetXaxis()->SetBinLabel(5,"4");
      hNBtagJets[j][i]->GetXaxis()->SetBinLabel(6,"#geq 5");
      if(j == 2) hNBtagJets[j][i]->GetXaxis()->SetRange(3,6);
      if(j == 3) hNBtagJets[j][i]->GetXaxis()->SetRange(5,6);
      
      h2DCSV_23Jet[j][i] = new TH2F("h2DCSV_23Jet_" + namech[i] + "_" + namecut[j], "CSVv2 Discriminant for 3rd and 4th Jets " + titlenamech[i], 10,0,1,10,0,1);

      hHT[j][i]         = new TH1F("hHT_"+namech[i]+"_"+namecut[j],"H_{T} " + titlenamech[i],100,0,600);
      hHT[j][i]->GetXaxis()->SetTitle("HT[GeV]");      

      /***************************
          SF(ID,ISO & Trigger)
      ***************************/
      hSFIDISO[j][i]           = new TH1F("hSFIDISO_"+namech[i]+"_"+namecut[j],"SF_{ID,ISO} " + titlenamech[i],400,0.8,1.2);    
      hSFIDISOError[j][i]      = new TH1F("hSFIDISOError_"+namech[i]+"_"+namecut[j],"#Delta SF_{ID,ISO} " + titlenamech[i],400,0,0.05); 
      hSFTrigger[j][i]         = new TH1F("hSFTrigger_"+namech[i]+"_"+namecut[j],"SF^{Trigger} " + titlenamech[i],400,0.8,1.2);    
      hSFTriggerError[j][i]    = new TH1F("hSFTriggerError_"+namech[i]+"_"+namecut[j],"#Delta SF^{Trigger} " + titlenamech[i],400,0,0.05);


      /***************************
              SF(b-tag)
      ***************************/
      h2DSFbtag_Global[j][i]     = new TH2F("h2DSFbtag_Global_"+namech[i]+"_"+namecut[j], "Global SF_{b-tag} Vs  #Delta SF_{b-tag} " + titlenamech[i], 40, 0.0, 4.0, 50, 0.0, 0.5);
      hSFbtag_Global[j][i]     = new TH1F("hSFbtag_Global_"+namech[i]+"_"+namecut[j], "Global SF_{b-tag} " + titlenamech[i],40.0, 0.0, 4.0);
      hSFbtag_Global_var[j][i] = new TH1F("hSFbtag_Global_var_"+namech[i]+"_"+namecut[j], "Global #Delta SF_{b-tag} " + titlenamech[i],20.0, 0.0, 1.0);

      h2DSFbtag_b[j][i]    = new TH2F("hSFbtag_b_"+namech[i]+"_"+namecut[j], "N^{b}(p_{T} vs #eta) " + titlenamech[i],7,0.0,140.0,4,0.0,2.4);
      h2DSFbtag_b[j][i]->GetXaxis()->SetTitle("p_{T}[GeV]");
      h2DSFbtag_b[j][i]->GetYaxis()->SetTitle("#eta");
      h2DSFbtag_c[j][i]    = new TH2F("hSFbtag_c_"+namech[i]+"_"+namecut[j], "N^{c}(p_{T} vs #eta) " + titlenamech[i],7,0.0,140.0,4,0.0,2.4);
      h2DSFbtag_c[j][i]->GetXaxis()->SetTitle("p_{T}[GeV]");
      h2DSFbtag_c[j][i]->GetYaxis()->SetTitle("#eta");
      h2DSFbtag_l[j][i]    = new TH2F("hSFbtag_l_"+namech[i]+"_"+namecut[j], "N^{l}(p_{T} vs #eta) " + titlenamech[i],7,0.0,140.0,4,0.0,2.4);
      h2DSFbtag_l[j][i]->GetXaxis()->SetTitle("p_{T}[GeV]");
      h2DSFbtag_l[j][i]->GetYaxis()->SetTitle("#eta");

      h2DSFbtag_btag_b[j][i]    = new TH2F("hSFbtag_btag_b_"+namech[i]+"_"+namecut[j], "N_{btag}^{b}(p_{T} vs #eta) " + titlenamech[i],7,0.0,140.0,4,0.0,2.4);
      h2DSFbtag_btag_b[j][i]->GetXaxis()->SetTitle("p_{T}[GeV]");
      h2DSFbtag_btag_b[j][i]->GetYaxis()->SetTitle("#eta");
      h2DSFbtag_btag_c[j][i]    = new TH2F("hSFbtag_btag_c_"+namech[i]+"_"+namecut[j], "N_{btag}^{c}(p_{T} vs #eta) " + titlenamech[i],7,0.0,140.0,4,0.0,2.4);
      h2DSFbtag_btag_c[j][i]->GetXaxis()->SetTitle("p_{T}[GeV]");
      h2DSFbtag_btag_c[j][i]->GetYaxis()->SetTitle("#eta");
      h2DSFbtag_btag_l[j][i]    = new TH2F("hSFbtag_btag_l_"+namech[i]+"_"+namecut[j], "N_{btag}^{l}(p_{T} vs #eta) " + titlenamech[i],7,0.0,140.0,4,0.0,2.4);
      h2DSFbtag_btag_l[j][i]->GetXaxis()->SetTitle("p_{T}[GeV]");
      h2DSFbtag_btag_l[j][i]->GetYaxis()->SetTitle("#eta");

            
      /***************************
            SF(pT Reweight)
      ***************************/
      hSFpT[j][i]           = new TH1F("hSFpT_"+namech[i]+"_"+namecut[j],"SF_{pT} " + titlenamech[i],500,0.4,1.4);    
      hSFpTError[j][i]      = new TH1F("hSFpTError_"+namech[i]+"_"+namecut[j],"#Delta SF_{pT} " + titlenamech[i],400,0,0.05); 
      
      
      TString jetn[4];
      jetn[0]= "Jet-0"; 
      jetn[1]= "Jet-1"; 
      jetn[2]= "Jet-2"; 
      jetn[3]= "Jet-3"; 
      
      for(int ij=0; ij<4; ij++){
	hCSV[ij][j][i] = new TH1F("hCSV_" + jetn[ij] + "_" + namech[i] + "_" + namecut[j],"CSV " + jetn[ij] + " " + titlenamech[i],10,0,1);
	hCSV[ij][j][i]->GetXaxis()->SetTitle("CSVv2");      
	hJetPt[ij][j][i] = new TH1F("hJetPt_" + jetn[ij] + "_" + namech[i] + "_" + namecut[j],"p_{T}^{Jet} " + jetn[ij] + " " + titlenamech[i],10,0,200);
	hJetPt[ij][j][i]->GetXaxis()->SetTitle("p_{T}[GeV]");      

	hJetpTUncVar[ij][j][i] = new TH1F("hJetpTUncVar_" + jetn[ij] + "_" + namech[i] + "_" + namecut[j], "#Delta pT^{Jet} " + jetn[ij] + " " + titlenamech[i], 20.0, 0.0, 2.0);
      }

      hEvtCatego[j][i]  = new TH1F("hEvtCatego_"+namech[i]+"_"+namecut[j],"ttbar Event Categorization " + titlenamech[i],4,-0.5,3.5);
      hEvtCatego[j][i]->GetXaxis()->SetTitle("ttbar");      

    }//for(i)
  }//for(j)
  

  TStopwatch sw;
  sw.Start(kTRUE);

  ///////////////////////////////////////
  // Please, IGNORE. Temporal solution //
  ///////////////////////////////////////
  TCanvas *mydummycanvas=new TCanvas();// 
  ///////////////////////////////////////
    
  /************************
     SF Parametrization
  *************************/

  TString fSFdir = "../TopTrees_CATuples/ScaleFactors/";
  
  TH2F *hmuIDISOSF, *hmuTriggerSF;
  TH2F *heIDISOSF,  *heTriggerSF;

  // Lepton SFs: ID and ISO with stat. + syst. Errors
  //TFile *MuSF = TFile::Open(fSFdir + "SF_muon_IDISO_13TeV_v2.root"); 
  TFile *MuSF = TFile::Open(fSFdir + "MuonSF_IDISO_Trigger_POG25ns.root"); 
  TFile *ElSF = TFile::Open(fSFdir + "ElectronSF_IDISO_Trigger_POG25ns.root"); 

  if(!MuSF || !ElSF){
    std::cerr << "ERROR [SF]: Could not open SF files!!!"  << std::endl;
    std::exit(0);
  }

  hmuIDISOSF = (TH2F*) MuSF->Get("GlobalSF")->Clone("muIDISOSF");
  hmuTriggerSF = (TH2F*) MuSF->Get("TriggerSF")->Clone("muTriggerSF"); 
  if(!hmuIDISOSF || !hmuTriggerSF){
    std::cerr << "ERROR [MuonSF]: Could not find histogram for SF reweighting" << std::endl;
  }

  heIDISOSF = (TH2F*) ElSF->Get("GlobalSF")->Clone("eIDISOSF");
  heTriggerSF = (TH2F*) ElSF->Get("TriggerSF")->Clone("eTriggerSF"); 
  if(!heIDISOSF || !heTriggerSF){
    std::cerr << "ERROR [ElectronSF]: Could not find histogram for SF reweighting" << std::endl;
  }

  // Trigger and ID-ISO uncertainties
  if(_idiso_unc){
    if(IDISOUnc=="Up")        fname += "_SYS_IDISO_Up";      
    else if(IDISOUnc=="Down") fname += "_SYS_IDISO_Down";
    else if(IDISOUnc=="Nom")  fname += "_SYS_IDISO_Nom";
  } // if(_idiso_unc)
  
  if(_tr_unc){
    if(TrUnc=="Up")        fname += "_SYS_Trigger_Up";
    else if(TrUnc=="Down") fname += "_SYS_Trigger_Down";
    else if(TrUnc=="Nom")  fname += "_SYS_Trigger_Nom";
  }// if(_tr_unc) 

  // Jet uncertainties (btag, JES and JER)
  if(_syst) fname += "_SYS_" + syst_varname;


  // Number de events for <pT Reweight>
  //          [Cut][Channel]
  float SF_pTweight[4][3]={0,0,0,0,
			   0,0,0,0,
			   0,0,0,0};
  // Number de events for acceptance
  //          [Cut][Channel]
  int AccEvent[4][3]={0,0,0,0,
		      0,0,0,0,
		      0,0,0,0};
  // Number de events for acceptance
  //          [Cut][Channel]
  float EffEvent[4][3]={0.0,0.0,0.0,0.0,
			0.0,0.0,0.0,0.0,
			0.0,0.0,0.0,0.0};
  
  /********************************
      pT Reweight; <SF_pT> 
   *******************************/
  // <SF_pT> per channel at single lepton level 
  float SF_tPt_mean[4][3];
  SF_tPt_mean[0][0]= 1.0;
  SF_tPt_mean[0][1]= 1.0;
  SF_tPt_mean[0][2]= 1.0;
  SF_tPt_mean[1][0]= 1.0;
  SF_tPt_mean[1][1]= 1.0;
  SF_tPt_mean[1][2]= 1.0;
  SF_tPt_mean[2][0]= 1.0;
  SF_tPt_mean[2][1]= 1.0;
  SF_tPt_mean[2][2]= 1.0;
  SF_tPt_mean[3][0]= 1.0;
  SF_tPt_mean[3][1]= 1.0;
  SF_tPt_mean[3][2]= 1.0;
  
  float SF_pT_TotalMean; // Over all cuts and channels
  float SF_n = 0.0;
  for(int SF_cut= 0; SF_cut<2; SF_cut++){
    for(int SF_ch= 0; SF_ch<3; SF_ch++){
      SF_pT_TotalMean +=  SF_tPt_mean[SF_cut][SF_ch];
      SF_n ++;
    }
  }
  if(SF_n != 0.0) SF_pT_TotalMean = SF_pT_TotalMean/SF_n; 
  

  /***************************
     ttbar Categorization
  ***************************/
  if(_ttbar_cat) fname += ttbar_id; // add in the sample name the ttbar category


  // New WP for 76X: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X
  float CSV_WP = 0.800; // Medium
  int btagSysPar = 0;
  // BTagSFUtil *fBTagSF;   //The BTag SF utility
  // BTagSFUtil *fBTagSFT;   //The BTag SF utility
  // fBTagSF = new BTagSFUtil("CSV", "Medium", btagSysPar); 
  // fBTagSF = new BTagSFUtil("CSV", "Tight", btagSysPar); 

  // Global SF uncertainty: 18 Components
  if     (_syst && syst_varname.Contains("btagjes_Up"))     btagSysPar = btagUnc::JES_UP;
  else if(_syst && syst_varname.Contains("btagjes_Down"))   btagSysPar = btagUnc::JES_DN;
  else if(_syst && syst_varname.Contains("btaglf_Up"))      btagSysPar = btagUnc::LF_UP;
  else if(_syst && syst_varname.Contains("btaglf_Down"))    btagSysPar = btagUnc::LF_DN;
  else if(_syst && syst_varname.Contains("btaghf_Up"))      btagSysPar = btagUnc::HF_UP;
  else if(_syst && syst_varname.Contains("btaghf_Down"))    btagSysPar = btagUnc::HF_DN;
  else if(_syst && syst_varname.Contains("btaghfsI_Up"))    btagSysPar = btagUnc::HFSTAT1_UP;
  else if(_syst && syst_varname.Contains("btaghfsI_Down"))  btagSysPar = btagUnc::HFSTAT1_DN;
  else if(_syst && syst_varname.Contains("btaghfsII_Up"))   btagSysPar = btagUnc::HFSTAT2_UP;
  else if(_syst && syst_varname.Contains("btaghfsII_Down")) btagSysPar = btagUnc::HFSTAT2_DN;
  else if(_syst && syst_varname.Contains("btaglfsI_Up"))    btagSysPar = btagUnc::LFSTAT1_UP;
  else if(_syst && syst_varname.Contains("btaglfsI_Down"))  btagSysPar = btagUnc::LFSTAT1_DN;
  else if(_syst && syst_varname.Contains("btaglfsII_Up"))   btagSysPar = btagUnc::LFSTAT2_UP;
  else if(_syst && syst_varname.Contains("btaglfsII_Down")) btagSysPar = btagUnc::LFSTAT2_DN;
  else if(_syst && syst_varname.Contains("btagcfI_Up"))     btagSysPar = btagUnc::CFERR1_UP;
  else if(_syst && syst_varname.Contains("btagcfI_Down"))   btagSysPar = btagUnc::CFERR1_DN;
  else if(_syst && syst_varname.Contains("btagcfII_Up"))    btagSysPar = btagUnc::CFERR2_UP;
  else if(_syst && syst_varname.Contains("btagcfII_Down"))  btagSysPar = btagUnc::CFERR2_DN;

  // Scale Uncertainty
  // From: https://indico.cern.ch/event/494682/contribution/2/attachments/1223578/1800218/mcaod-Feb15-2016.pdf
  int scaleSysPar;
  if     (_syst && syst_varname.Contains("ScaleRnF_Up"))   scaleSysPar = 0; // muR=Nom,  muF=Up
  else if(_syst && syst_varname.Contains("ScaleRnF_Down")) scaleSysPar = 1; // muR=Nom,  muF=Down
  else if(_syst && syst_varname.Contains("ScaleRuF_Nom"))  scaleSysPar = 2; // muR=Up,   muF=Nom
  else if(_syst && syst_varname.Contains("ScaleRuF_Up"))   scaleSysPar = 3; // muR=Up,   muF=Up
  else if(_syst && syst_varname.Contains("ScaleRdF_Nom"))  scaleSysPar = 4; // muR=Down, muF=Nom
  else if(_syst && syst_varname.Contains("ScaleRdF_Down")) scaleSysPar = 5; // muR=Down, muF=Down
  // Normalization for Scale Weights:
  if (_syst && syst_varname.Contains("Scale")){
    if(scaleSysPar < 6) nNorm_Event = NTotal_ScalemuRF_Weight[scaleSysPar];
    else{
      std::cerr << "No entry for Scale normalization! Check HISTO!"  << std::endl;
      std::exit(0);
    }
  }
  // PileUp Uncertainty  
  int pileupSysPar;
  if     (_syst && syst_varname.Contains("PileUp_Up"))   pileupSysPar = 1; // Up
  else if(_syst && syst_varname.Contains("PileUp_Down")) pileupSysPar = 2; // Down

  /************************
    Normalization Weights
  *************************/  
  float NormWeight = 0.0;
  // NormWeight = Lumi*(1.0/N_Gen_events)*(Xsec)
  NormWeight = SFLumi(fname, 2260, nNorm_Event);  

  std::cout << "-----------------------                                 -------------------------" << std::endl;
  std::cout << "Number of Events     = " << nNorm_Event << std::endl;
  std::cout << "Normalization Factor = " << NormWeight  << std::endl;
  std::cout << "---------------------------------------------------------------------------------" << std::endl;


  /********************************
             Event Loop
  ********************************/
  std::cout << "--- Processing: " << theTree.GetEntries() << " events" << std::endl;
  
  for (Long64_t ievt=0; ievt<theTree.GetEntries();ievt++) {
    
    theTree.GetEntry(ievt);  
    print_progress(theTree.GetEntries(), ievt);

    // PU reweight: Includes Syst. Unc.
    if (_syst && syst_varname.Contains("PileUp"))
      PUWeight = (*PUWeight_sys)[pileupSysPar]; // Up
    else PUWeight = (*PUWeight_sys)[0];
    
    // MCatNLO GEN Weights (For MC@NLO)
    PUWeight = PUWeight * GENWeight;
    // Normalization Weight
    PUWeight = PUWeight * NormWeight;
    
    // Scale reweight: Syst. Unc.
    if (_syst && syst_varname.Contains("ScaleRF"))
      PUWeight = PUWeight*(*ScaleWeight)[scaleSysPar];
    
    
    int NJets,NBtagJets, NBtagTJets;
    
    TLorentzVector Lep;
    std::vector<int>  JetIndex;
    std::vector<bool> bJet;
    
    Lep.SetPxPyPzE(Lep_px,Lep_py,Lep_pz,Lep_E);
    if(Lep.Pt() < 30)  continue; // Lep pT >30GeV
    
    // Transverse W Mass
    TLorentzVector METv;
    METv.SetPtEtaPhiE(MET,0.0,MET_Phi,MET);
    
    float mT = sqrt(2*MET*Lep.Pt()*(1.0-cos( Lep.DeltaPhi(METv) )));
    
    // Jets 
    NJets      = 0;
    NBtagJets  = 0;
    
    // Global SF_b-tag
    float btagUnc_val = 0.0;
    // From: https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration
    //btagUnc btagvar;
    if (!fname.Contains("Data")){
      if(_syst && syst_varname.Contains("btag")){
	if(syst_varname.Contains("Up"))
	  btagUnc_val = 1.0 * (*Jet_SF_CSV)[btagSysPar];
	if(syst_varname.Contains("Down")) 
	  btagUnc_val = -1.0 * (*Jet_SF_CSV)[btagSysPar];
	
	PUWeight = PUWeight * ((*Jet_SF_CSV)[btagUnc::CENTRAL] + btagUnc_val);
	
      } // if(_syst && btag)
      else {
	// SF estimated for jets with pT > 25GeV
	PUWeight = PUWeight * (*Jet_SF_CSV)[btagUnc::CENTRAL];
      }
    }// if(!data)
    
    for(int ijet=0; ijet < Jet_px->size(); ijet++){
      
      TLorentzVector jet;
      jet.SetPxPyPzE((*Jet_px)[ijet],(*Jet_py)[ijet],(*Jet_pz)[ijet],(*Jet_E)[ijet]);
      float jet_pT = jet.Pt(); 
      int JetFlav  = (*Jet_partonFlavour)[ijet];      
      
      float JetSystVar = 1.0;
      if(_syst){
	if(syst_varname.Contains("JES_Up")){
	  JetSystVar = (*Jet_JES_Up)[ijet];
	}
	else if(syst_varname.Contains("JES_Down")){
	  JetSystVar = (*Jet_JES_Down)[ijet];
	}
	else if(syst_varname.Contains("JER_Up")){
	  JetSystVar = (*Jet_JER_Up)[ijet];
	}
	else if(syst_varname.Contains("JER_Nom")){
	  JetSystVar = (*Jet_JER_Nom)[ijet];
	}
	else if(syst_varname.Contains("JER_Down")){
	  JetSystVar = (*Jet_JER_Down)[ijet];
	}
      }
      jet_pT = jet_pT*JetSystVar;

      if(jet_pT>25){ // Jet pT > 30GeV
	
	JetIndex.push_back(ijet);
	NJets++; // Number of jets

	/*******************************************
                       b-tagging
	*******************************************/    
	bool btagDisc = false;
		
	// OLD Method: 
	// b-tagging WP from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagging
	// Not recommended for analysis with observables that depends ontagging.
	// if (fname.Contains("Data")) btagDisc = fBTagSF->IsTagged((*Jet_CSV)[ijet], -999999, jet.Pt(), jet.Eta());
	// else btagDisc = fBTagSF->IsTagged((*Jet_CSV)[ijet], JetFlav, jet.Pt(), jet.Eta());
	
	// New Method (Event SF from tth group)
	// https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration
	btagDisc = (*Jet_CSV)[ijet] > CSV_WP;

	bJet.push_back(btagDisc);
	if(btagDisc) NBtagJets++; // Number of b-tagged jets
		
      } // if(Jet_pT)
    }// for(jets)

    // Number of GENJets
    int NGenJets = 0;
    if (_ttbar_cat){
      for(int igenjet=0; igenjet < GenJet_pT->size(); igenjet++){
	if( (*GenJet_pT)[igenjet] > 20 ) NGenJets++;
      }
    }

  /*******************************************
   Trigger,ID & ISO Scale Factors/bin(Pt,Eta)
  *******************************************/    
    std::vector<float> SF_ID_ISO_Tr;

    if (fname.Contains("Data")){
      SF_ID_ISO_Tr.push_back(1.0); // SF_ID_ISO_Tr    [0] 
      SF_ID_ISO_Tr.push_back(1.0); // SF_ID_ISO       [1] 
      SF_ID_ISO_Tr.push_back(1.0); // SF_ID_ISO_Error [2] 
      SF_ID_ISO_Tr.push_back(1.0); // SF_Tr           [3] 
      SF_ID_ISO_Tr.push_back(1.0); // SF_Tr_Error     [4] 
    }
    
    else {
      // Second Method: Taking SF from root file
      SFIDISOTrigger(SF_ID_ISO_Tr,
      		     Lep, Channel,
      		     hmuIDISOSF, hmuTriggerSF,
      		     heIDISOSF,  heTriggerSF);
      
      // Easiest adaptation: Create Lep_SF in the same way the Lep_SF vector
      // SF_ID_ISO_Tr = Lep_SF;

      if(_idiso_unc){
	if     (IDISOUnc == "Up")   PUWeight = PUWeight * (SF_ID_ISO_Tr[1] + SF_ID_ISO_Tr[2]);
	else if(IDISOUnc == "Down") PUWeight = PUWeight * (SF_ID_ISO_Tr[1] - SF_ID_ISO_Tr[2]);
	else if(IDISOUnc == "Nom")  PUWeight = PUWeight * (SF_ID_ISO_Tr[1]);
      } // if(_idiso_unc)
      
      else if(_tr_unc){
	if     (TrUnc=="Up")   PUWeight=PUWeight*(SF_ID_ISO_Tr[3] + SF_ID_ISO_Tr[4]);
	else if(TrUnc=="Down") PUWeight=PUWeight*(SF_ID_ISO_Tr[3] - SF_ID_ISO_Tr[4]);
	else if(TrUnc=="Nom")  PUWeight=PUWeight*(SF_ID_ISO_Tr[3]);	
      }// if(_tr_unc) 
      
      // Check the SF_e implementation
      // If the electron is in the transition region, SF = 1
      else  if(SF_ID_ISO_Tr[0] != 0.0)  PUWeight=PUWeight*(SF_ID_ISO_Tr[0]); 
      
    }// else(Contain("Data"))
    

    /***************************
            Selection
    ***************************/

    int                            cut = 0; // Single Lepton (from Tree)
    if(NJets > 5)                  cut = 1; // + 6 Jets 
    if(NJets > 5 && NBtagJets > 1) cut = 2; // + 2 b-tag
    if(NJets > 5 && NBtagJets > 3) cut = 3; // + 4 b-tag

    /***************************
        ttbar Categorization
     ***************************/
    if (_ttbar_cat){
      // Categorization using Cone DeltaR
      // Visible Phase Space: 
      int cone_NJets  = (*GenConeCat)[1];
      int cone_NbJets = (*GenConeCat)[2];
      int cone_NcJets = (*GenConeCat)[3];
      
      int cone_NbJetsNoTop = (*GenConeCat)[4];
      
      // Full Phase Space:
      int cone_NaddJets  = (*GenConeCat)[5];
      int cone_NaddbJets = (*GenConeCat)[6];
      
      bool Isttjj = false;
      bool Isttbb = false;
      bool Isttcc = false;
      bool Isttb  = false;
      bool IsttLF = false;
      bool Istt   = false;
      
      if(cone_NbJets > 1 && cone_NJets > 5) Isttjj = true;
      
      if      (cone_NbJets > 3  && cone_NJets > 5) Isttbb = true;
      else if (cone_NbJets > 2  && cone_NJets > 5) Isttb  = true;
      else if (cone_NbJets > 1  && cone_NJets > 5 && cone_NcJets > 1) Isttcc = true;
      else if (cone_NbJets > 1  && cone_NJets > 5) IsttLF = true;
      else Istt = true;
      
      if(ttbar_id == "ttjj" && !Isttjj) cut = -1;
      if(ttbar_id == "ttbb" && !Isttbb) cut = -1;
      if(ttbar_id == "ttb"  && !Isttb ) cut = -1;
      if(ttbar_id == "ttcc" && !Isttcc) cut = -1;
      if(ttbar_id == "ttLF" && !IsttLF) cut = -1;
      if(ttbar_id == "tt"   && !Istt)   cut = -1;

    }// if(_ttbar_cat)
    
    /***************************
        Data and QCD Samples
     ***************************/
    if (fname.Contains("DataSingleMu") && Channel==1)  cut = -1;
    if (fname.Contains("DataSingleEG") && Channel==0)  cut = -1;
    if (fname.Contains("QCD_MuEnr")    && Channel==1)  cut = -1;
    if (fname.Contains("QCD_EGEnr")    && Channel==0)  cut = -1;

    /***************************
          Loop over cuts
    ***************************/
    for(int icut = 0; icut < (cut+1); icut++){
            
      /************************************************************
        pT reweight (Only for ttbar signal)
	It shouldn't be applied to the central value. 
        It is divided over the <SF_pT> (normalization)  
      *************************************************************/
      if(fname.Contains("ttbar")){	
	// // pT Top Reweight
	// TLorentzVector t,tbar;
	// t.SetPxPyPzE(tPx,tPy,tPz,tE);
	// tbar.SetPxPyPzE(tbarPx,tbarPy,tbarPz,tbarE);
	
	// // pT Reweight: Only for Systematic Uncertainty
	// // Definition from pT reweight Twiki
	// // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting and 
	// // https://twiki.cern.ch/twiki/bin/view/CMS/TopSystematicsRun1#pt_top_Reweighting
	// float SF_tPt=sqrt( exp(0.222-0.00197*TMath::Min(t.Pt(), 400.)) * exp(0.222-0.00197*TMath::Min(tbar.Pt(), 400.)) );
	// hSFpT[icut][Channel]->Fill(SF_tPt,PUWeight); 	
	
	// if(fname.Contains("TTJets_Dilep_SYS_pTReweight_Up") || fname.Contains("TTJets_Dilep_SYS_pTReweight_Down")) PUWeight=PUWeight*(SF_tPt/SF_pT_TotalMean); 
	// SF_pTweight[icut][Channel]+=SF_tPt;
	SF_pTweight[icut][Channel] += 1.0;	
      } // if(TTbar) 
      /*************************************************************/
      
      /*******************
        Fill Histograms
      *******************/
      hSFIDISO[icut][Channel]->Fill(SF_ID_ISO_Tr[1],PUWeight);
      hSFIDISOError[icut][Channel]->Fill(SF_ID_ISO_Tr[2],PUWeight);
      hSFTrigger[icut][Channel]->Fill(SF_ID_ISO_Tr[3],PUWeight);
      hSFTriggerError[icut][Channel]->Fill(SF_ID_ISO_Tr[4],PUWeight);
    
      /******************
          Acceptace
      ******************/
      AccEvent[icut][Channel]++;
      EffEvent[icut][Channel]= EffEvent[icut][Channel] + PUWeight;

      /******************
        Kinematic Var.
      ******************/
      
      hPV[icut][Channel]->Fill(GoodPV,PUWeight);

      hMET[icut][Channel]->Fill(MET,PUWeight);
      hMET_Phi[icut][Channel]->Fill(fabs(MET_Phi),PUWeight);
      
      hLepPt[icut][Channel]->Fill(Lep.Pt(),PUWeight);
      hLepEta[icut][Channel]->Fill(fabs(Lep.Eta()),PUWeight);
      hLepPhi[icut][Channel]->Fill(fabs(Lep.Phi()),PUWeight);

      hmT[icut][Channel]->Fill(mT,PUWeight);
      
      /******************
          Jets Var.
      ******************/
      hNJets[icut][Channel]->Fill(NJets,PUWeight); 
      hNBtagJets[icut][Channel]->Fill(NBtagJets,PUWeight);

      // Global btag SF
      h2DSFbtag_Global[icut][Channel]->Fill((*Jet_SF_CSV)[btagUnc::CENTRAL], btagUnc_val, PUWeight);
      hSFbtag_Global[icut][Channel]->Fill((*Jet_SF_CSV)[btagUnc::CENTRAL], PUWeight);
      hSFbtag_Global_var[icut][Channel]->Fill(btagUnc_val, PUWeight);
      
      // CSV discriminant for 3rd and 4th Jet
      if(JetIndex.size() > 3)  h2DCSV_23Jet[icut][Channel]->Fill((*Jet_CSV)[JetIndex[2]], (*Jet_CSV)[JetIndex[3]], PUWeight);


      for(int ijet=0; ijet < JetIndex.size(); ijet++){
	TLorentzVector jet;
	jet.SetPxPyPzE((*Jet_px)[JetIndex[ijet]],(*Jet_py)[JetIndex[ijet]],(*Jet_pz)[JetIndex[ijet]],(*Jet_E)[JetIndex[ijet]]);
	int JetFlav = (*Jet_partonFlavour)[JetIndex[ijet]];
	bool IsbJet = bJet[ijet];

	if (ijet < 4){
	  hCSV  [ijet][icut][Channel]->Fill((*Jet_CSV)[JetIndex[ijet]], PUWeight);
	  hJetPt[ijet][icut][Channel]->Fill(jet.Pt(), PUWeight);
	}

	if(JetFlav == 5){
	  h2DSFbtag_b[icut][Channel]->Fill(jet.Pt(), fabs(jet.Eta()), PUWeight); // b-Flavour
	  if((*Jet_CSV)[JetIndex[ijet]] > CSV_WP) h2DSFbtag_btag_b[icut][Channel]->Fill(jet.Pt(), fabs(jet.Eta()), PUWeight);  
	}	
	else if(JetFlav == 4){
	  h2DSFbtag_c[icut][Channel]->Fill(jet.Pt(), fabs(jet.Eta()), PUWeight); // c-Flavour
	  if((*Jet_CSV)[JetIndex[ijet]] > CSV_WP) h2DSFbtag_btag_c[icut][Channel]->Fill(jet.Pt(), fabs(jet.Eta()), PUWeight);  
	}
	else if(JetFlav == 1 || JetFlav == 2 || JetFlav == 3){
	  h2DSFbtag_l[icut][Channel]->Fill(jet.Pt(), fabs(jet.Eta()), PUWeight); // l-Flavour
	  if((*Jet_CSV)[JetIndex[ijet]] > CSV_WP) h2DSFbtag_btag_l[icut][Channel]->Fill(jet.Pt(), fabs(jet.Eta()), PUWeight);  
	}
      }//for(ijet)     

    }//for(icuts)     
    
    JetIndex.clear();
    bJet.clear();
    
  }//for(events)
  

  /***************************
      TTbar pT Reweight
  ***************************/

  // if(fname.Contains("ttbar")){
    
  //   for(int icut=0; icut<4;icut++){
  //     for(int ich=0; ich<2;ich++){
  // 	// SF_pT_reweight mean value 
  // 	SF_pTweight[icut][ich]=SF_pTweight[icut][ich]/AccEvent[icut][ich]; 
  // 	std::cout << "<SF_pT[cut=" << icut << "][ch=" <<ich << "]>= " << SF_pTweight[icut][ich] << std::endl;	
  //     } // for(icut)
  //   } // for(ich)
  // } // if(TT)


  // Get elapsed time
  sw.Stop();
  std::cout << "==================================================] 100% " << std::endl;
  std::cout << "--- End of event loop: "; sw.Print();
  

  //Acceptance-Efficiency
  std::cout << "--------  Acceptace  --------" << std::endl;
  std::cout << "Number of RAW-mu+Jets events:" << std::endl;
  std::cout << namecut[0] << ": " << AccEvent[0][0] << std::endl;
  std::cout << namecut[1] << ": " << AccEvent[1][0] << std::endl;
  std::cout << namecut[2] << ": " << AccEvent[2][0] << std::endl;
  std::cout << namecut[3] << ": " << AccEvent[3][0] << std::endl;

  std::cout << "--------  Efficiency  --------" << std::endl;
  std::cout << "Number of Weigthed-mu+Jets events:" << std::endl;
  std::cout << namecut[0] << ": " << EffEvent[0][0] << " +/- " << sqrt(AccEvent[0][0])*EffEvent[0][0]/AccEvent[0][0] << std::endl;
  std::cout << namecut[1] << ": " << EffEvent[1][0] << " +/- " << sqrt(AccEvent[1][0])*EffEvent[1][0]/AccEvent[1][0] << std::endl;
  std::cout << namecut[2] << ": " << EffEvent[2][0] << " +/- " << sqrt(AccEvent[2][0])*EffEvent[2][0]/AccEvent[2][0] << std::endl;
  std::cout << namecut[3] << ": " << EffEvent[3][0] << " +/- " << sqrt(AccEvent[3][0])*EffEvent[3][0]/AccEvent[3][0] << std::endl;
  
  std::cout << "--------  Acceptace  --------" << std::endl;
  std::cout << "Number of RAW-e+Jets events:" << std::endl;
  std::cout << namecut[0] << ": " << AccEvent[0][1] << std::endl;
  std::cout << namecut[1] << ": " << AccEvent[1][1] << std::endl;
  std::cout << namecut[2] << ": " << AccEvent[2][1] << std::endl;
  std::cout << namecut[3] << ": " << AccEvent[3][1] << std::endl;

  std::cout << "--------  Efficiency  --------" << std::endl;
  std::cout << "Number of Weigthed-e+Jets events: " << std::endl;
  std::cout << namecut[0] << ": " << EffEvent[0][1] << " +/- " << sqrt(AccEvent[0][1])*EffEvent[0][1]/AccEvent[0][1] << std::endl;
  std::cout << namecut[1] << ": " << EffEvent[1][1] << " +/- " << sqrt(AccEvent[1][1])*EffEvent[1][1]/AccEvent[1][1] << std::endl;
  std::cout << namecut[2] << ": "  << EffEvent[2][1] << " +/- " << sqrt(AccEvent[2][1])*EffEvent[2][1]/AccEvent[2][1] << std::endl;
  std::cout << namecut[3] << ": " << EffEvent[3][1] << " +/- " << sqrt(AccEvent[3][1])*EffEvent[3][1]/AccEvent[3][1] << std::endl;


  //Output Dir
  TString dirname="TopResults";   
  // make a dir if it does not exist!!
  struct stat st;
  if(stat(dirname,&st) != 0) system("mkdir " + dirname);
  

  // Sample name identification
  TString samplename="";
  bool matchsamplename=false;
  
  for(int i=0; i<fname.Sizeof(); i++){
    if (i>3){
      if (fname[i-4]=='b' && 
	  fname[i-3]=='-' && 
	  fname[i-2]=='1' && 
	  fname[i-1]=='_') matchsamplename=true;
    }
    if (matchsamplename) samplename.Append(fname[i]);
  }

  
  if(!_syst){
    // Yields
    // make a dir if it does not exist!!
    TString diryieldsname = dirname + "/Yields_" + hname;
    struct stat st;
    if(stat(diryieldsname,&st) != 0) system("mkdir " + diryieldsname);
    
    TString Yieldfile = diryieldsname + "/" + samplename.Data() + ".h";
    //Yieldfile += ".h";
    // Option a = append: Open file for output at the end of a file.
    // Option w = write: Create an empty file for output operations. 
    FILE* fyields = fopen(Yieldfile, "w"); 

    fprintf(fyields,"\n///////////////////////////////////////////////////////////////////////////////// \n\n");
    fprintf(fyields,"// %s Sample on %s \n", (fname + ".root").Data() , currentDateTime().Data());
    fprintf(fyields,"// %s version \n", hname.Data());
    fprintf(fyields," float  %s[16][3][2][4];//[systematic][variation][channel][Cut] \n",      samplename.Data());
    fprintf(fyields," float  err_%s[16][3][3][4]; //[systematic][variation][channel][Cut] \n", samplename.Data());
    fprintf(fyields,"// Systematic: [0]=Trigger [1]=ID-ISO [2]=LES \n");
    fprintf(fyields,"// Systematic: [3]=JES [4]=JER [5]=b-tag [6]=PileUp \n");
    fprintf(fyields,"// Systematic: [7]=Scale [8]=Matching [9]=pTreweight [9]=PDF \n");
    fprintf(fyields,"// Systematic: [10]=DY-DD [11]=Non-W/Z [15]=Nom \n");
    fprintf(fyields,"// Variation: [0]=Up [1]=Down [2]=Nom \n");
    fprintf(fyields,"// Channel: [0]=mumu [1]=ee [2]=mue \n");
    fprintf(fyields,"// Cut: [0]=Dielpton [1]=Jets+Zveto [2]=MET [3]=btag \n");

    for(int ch=0;ch<2;ch++){
      for(int cut=0;cut<4;cut++){
      fprintf(fyields,"%s[15][2][%i][%i] = %.3f ; \n", samplename.Data(), ch, cut, EffEvent[cut][ch]);
      if(AccEvent[cut][ch]!=0.0) fprintf(fyields,"err_%s[15][2][%i][%i] = %.3f ; \n", samplename.Data(), ch, cut, sqrt(AccEvent[cut][ch])*EffEvent[cut][ch]/AccEvent[cut][ch]);
      else fprintf(fyields,"err_%s[15][2][%i][%i] = 0.0 ; \n", samplename.Data(), ch, cut);
      }
    }
    fclose(fyields);
    
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << "Yields saved into " << Yieldfile << " file" << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;    
  } 
  
  if (_syst){    

    // Name change
    TString processname = samplename;
    TString systname    = "";
    bool matchsystname  = false; 
    
    for(int i=0;i<samplename.Sizeof();i++){
      if (i>3){
	if (samplename[i-4]=='S' && 
	    samplename[i-3]=='Y' && 
	    samplename[i-2]=='S' && 
	    samplename[i-1]=='_'){
	  matchsystname = true;
	}
      }
      if (matchsystname) systname.Append(samplename[i]);
    }
    
    processname.ReplaceAll(("_SYS_" + systname).Data(), "");  
    
    int variation  = -999;
    int systsource = -999;
    
    if     (systname.Contains("Up"))   variation = 0;
    else if(systname.Contains("Down")) variation = 1;
    else if(systname.Contains("Nom"))  variation = 2;

    if(systname.Contains("Trigger"))         systsource = 0;
    else if(systname.Contains("IDISO"))      systsource = 1;
    else if(systname.Contains("LES"))        systsource = 2;
    else if(systname.Contains("JES"))        systsource = 3;
    else if(systname.Contains("JER"))        systsource = 4;
    else if(systname.Contains("PileUp"))     systsource = 6;
    else if(systname.Contains("Scale"))      systsource = 7;
    else if(systname.Contains("Matching"))   systsource = 8;
    else if(systname.Contains("pTReweight")) systsource = 9;
    else if(systname.Contains("Powheg"))     systsource = 10;
    else if(systname.Contains("btagjes"))    systsource = 11;
    else if(systname.Contains("btaglf"))     systsource = 12;
    else if(systname.Contains("btaghf"))     systsource = 13;
    else if(systname.Contains("btaghfsI"))   systsource = 14;
    else if(systname.Contains("btaghfsII"))  systsource = 16;
    else if(systname.Contains("btaglfsI"))   systsource = 17;
    else if(systname.Contains("btaglfsII"))  systsource = 18;
    else if(systname.Contains("btagcfI"))    systsource = 19;
    else if(systname.Contains("btagcfII"))   systsource = 20;

    // Systematic Uncertainty Estimations
    // make a dir if it does not exist!!
    TString dirSysyieldsname = dirname + "/SysYields_" + hname;
    struct stat st;
    if(stat(dirSysyieldsname,&st) != 0) system("mkdir " + dirSysyieldsname);
    
    TString Syst_Yieldfile = dirSysyieldsname + "/" + samplename.Data() + ".h";
    FILE* fSys = fopen(Syst_Yieldfile, "w");        
    
    fprintf(fSys,"\n///////////////////////////////////////////////////////////////////////////////// \n\n");
    fprintf(fSys,"// %s Sample on %s \n", (fname + ".root").Data(), currentDateTime().Data());
    fprintf(fSys,"// %s version \n", hname.Data());
    fprintf(fSys,"// float  %s[16][3][3][4];//[systematic][variation][channel][Cut] \n", processname.Data());
    fprintf(fSys,"// float  err_%s[16][3][3][4]; //[systematic][variation][channel][Cut] \n", processname.Data());
    fprintf(fSys,"// Systematic: [0]=Trigger [1]=ID-ISO [2]=LES \n");
    fprintf(fSys,"// Systematic: [3]=JES [4]=JER [5]=b-tag [6]=PileUp \n");
    fprintf(fSys,"// Systematic: [7]=Scale [8]=Matching [9]=pTreweight [9]=PDF \n");
    fprintf(fSys,"// Systematic: [10]=DY-DD [11]=Non-W/Z [15]=Nom \n");
    fprintf(fSys,"// Variation: [0]=Up [1]=Down [2]=Nom \n");
    fprintf(fSys,"// Channel: [0]=mumu [1]=ee [2]=mue \n");
    fprintf(fSys,"// Cut: [0]=Dielpton [1]=Jets+Zveto [2]=MET [3]=btag \n");

    for(int ch=0;ch<3;ch++){
      for(int cut=0;cut<4;cut++){
	fprintf(fSys,"%s     [%i][%i][%i][%i] = %.3f ; \n", processname.Data(), systsource, variation, ch, cut, EffEvent[cut][ch]);
	if(AccEvent[cut][ch]!=0.0) fprintf(fSys,"err_%s [%i][%i][%i][%i] = %.3f ; \n", processname.Data(), systsource, variation, ch, cut, sqrt(AccEvent[cut][ch])*EffEvent[cut][ch]/AccEvent[cut][ch]);
	else fprintf(fSys,"err_%s [%i][%i][%i][%i] = 0.0 ; \n", processname.Data(), systsource, variation, ch, cut);
      }
    }

    fclose(fSys);
    
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << "Yields saved for Syst. estimation into " << Syst_Yieldfile << " file" << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;    
    
  }
  
  
  // --- Write histograms
  
  TString outfname=dirname + "/hSF-" + hname + "_" + fname + ".root";
  TFile *target  = new TFile(outfname,"RECREATE" );  

  for(int j=0; j<4; j++){
    for(int i=0; i<2; i++){
      
      hEvtCatego[j][i]->Write();
      hPV[j][i]->Write();
      
      hMET[j][i]->Write();
      hMET_Phi[j][i]->Write();

      hmT[j][i]->Write();
      
      hLepPt[j][i]->Write();
      hLepEta[j][i]->Write();
      hLepPhi[j][i]->Write();
      
      hNJets[j][i]->Write();
      hNBtagJets[j][i]->Write();            
      h2DCSV_23Jet[j][i]->Write();
      
      h2DSFbtag_Global[j][i]->Write();
      hSFbtag_Global[j][i]->Write();
      hSFbtag_Global_var[j][i]->Write();

      for(int ij=0; ij<4; ij++){
	hCSV[ij][j][i]->Write();
	hJetPt[ij][j][i]->Write();
      }
      
      hSFpT[j][i]->Write();
      hSFpTError[j][i]->Write();
      
      hSFIDISO[j][i]->Write();
      hSFIDISOError[j][i]->Write();
      hSFTrigger[j][i]->Write();
      hSFTriggerError[j][i]->Write();
      
      h2DSFbtag_b[j][i]->Write();
      h2DSFbtag_c[j][i]->Write();
      h2DSFbtag_l[j][i]->Write();

      h2DSFbtag_btag_b[j][i]->Write();
      h2DSFbtag_btag_c[j][i]->Write();
      h2DSFbtag_btag_l[j][i]->Write();

    }//for(i)

  }//for(j)
  
  std::cout << "File saved as " << outfname << std::endl;

}


#endif

