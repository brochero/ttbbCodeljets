#ifndef __CINT__

#include<string>
#include<iostream>
#include<sstream>
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<set>
#include<vector>
#include <sys/stat.h>

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
#include "TError.h"
// TMVA
#include "TMVA/MethodCategory.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include <TMVAClassification_BDT.class.C>

// TopCode
#include <SFIDISOTrigger.h> // SF_ID-ISO-Trigger
#include <ttbar_category.h> // Event Categorization
#include <SFLumi.h>         // Normalization SF
#include <MVATraining.h>    // Dijet MVA Training

#ifndef __CINT__

// Jet Class
class ComJet: public TLorentzVector{
public:
  float CSV, CvsL, CvsB;
  int Flavour, pTIndex, Mom;
};

void print_progress(int TreeEntries, Long64_t ievt);
const TString currentDateTime();
float DiJetMassCorrection(std::vector<ComJet> &Jets, bool ReArrange);
std::vector<std::vector<float>> DiJetVariableEstimation(std::vector<ComJet>  &Jets, float DRAddJets, bool MVA, IClassifierReader *classReader);


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

// b-tagging SF
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
  const char * _dir      = "../Files_v7-6-4/";
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
  
  //Output Dir
  TString dirname="TopResults";   
  // make a dir if it does not exist!!
  struct stat st;
  if(stat(dirname,&st) != 0) system("mkdir " + dirname);


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
  float Lep_LES=0;

  std::vector<float> *Jet_px=0, *Jet_py=0, *Jet_pz=0, *Jet_E=0;
  std::vector<int>   *Jet_partonFlavour=0;
  std::vector<int>   *Jet_pTIndex=0;
  std::vector<int>   *Jet_GENmatched=0, *GenJet_Mom=0;
  std::vector<float> *Jet_CSV=0;
  std::vector<float> *Jet_SF_CSV=0;
  std::vector<float> *Jet_CSVCvsL=0;
  std::vector<float> *Jet_CvsL=0, *Jet_CvsB=0;
  std::vector<float> *Jet_JER_Up=0, *Jet_JER_Nom=0, *Jet_JER_Down=0;
  std::vector<float> *Jet_JES_Up=0, *Jet_JES_Down=0;

  // GEN Info
  std::vector<int> *GenConeCat=0;
  float DRAddJets;
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
  theTree.SetBranchAddress( "lepton_LES",  &Lep_LES );

  theTree.SetBranchAddress( "jet_px", &Jet_px );
  theTree.SetBranchAddress( "jet_py", &Jet_py );
  theTree.SetBranchAddress( "jet_pz", &Jet_pz );
  theTree.SetBranchAddress( "jet_E",  &Jet_E );

  theTree.SetBranchAddress( "jet_index",  &Jet_pTIndex );

  theTree.SetBranchAddress( "jet_CSV",  &Jet_CSV );
  theTree.SetBranchAddress( "jet_SF_CSV",  &Jet_SF_CSV );
  theTree.SetBranchAddress( "jet_partonFlavour",  &Jet_partonFlavour );

  theTree.SetBranchAddress( "jet_iCSVCvsL", &Jet_CSVCvsL );
  theTree.SetBranchAddress( "jet_CCvsLT",  &Jet_CvsL );
  theTree.SetBranchAddress( "jet_CCvBLT",  &Jet_CvsB );


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
    theTree.SetBranchAddress("draddjets", &DRAddJets);
    theTree.SetBranchAddress("genlepton_pT", &GenLep_pT);
    theTree.SetBranchAddress("jet_MatchedGenJetIndex", &Jet_GENmatched);
    theTree.SetBranchAddress("genjet_mom", &GenJet_Mom);
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
  TH1F *hCSV[6][4][2], *hJetPt[6][4][2], *hJetpTUncVar[6][4][2];
  TH1F *hCvsL[6][4][2], *hCvsB[6][4][2]; 
  TH2F *h2DCSV_23Jet[4][2], *h2DCSV_45Jet[4][2]; 
  TH2F *h2DCSV_24Jet[4][2], *h2DCSV_25Jet[4][2];
  TH2F *h2DCSV_34Jet[4][2], *h2DCSV_35Jet[4][2];
  TH1F *hMassJet[5][6][4][2];
  TH1F *hInvMassjj[4][2], *hMaxMVAjj[4][2];

  TH1F *hSFpT[4][2], *hSFpTError[4][2];
  TH1F *hSFIDISO[4][2], *hSFIDISOError[4][2];
  TH1F *hSFTrigger[4][2], *hSFTriggerError[4][2];
  TH2F *h2DSFbtag_Global[4][2];
  TH1F *hSFbtag_Global[4][2], *hSFbtag_Global_var[4][2];
  TH2F *h2DSFbtag_b[4][2], *h2DSFbtag_c[4][2], *h2DSFbtag_l[4][2], *h2DSFbtag_btag_b[4][2], *h2DSFbtag_btag_c[4][2], *h2DSFbtag_btag_l[4][2]; 

  TH1F *hEvtCatego[4][2];

  TH1F *hTJetPosition, *hWJetPosition, *hOJetPosition;
  TH2F *h2DTJetPosition, *h2DWJetPosition, *h2DOJetPositionMVA;

  TString namech[3];
  namech[0]="mujets";
  namech[1]="ejets";
  namech[2]="lepjets";
  
  TString namecut[4];
  namecut[0]="lepton";
  namecut[1]="6Jets";
  namecut[2]="2btag";
  namecut[3]="3btag";
  
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
      
      h2DCSV_23Jet[j][i] = new TH2F("h2DCSV_23Jet_" + namech[i] + "_" + namecut[j], "CSVv2 Discriminant for 3rd and 4th Jets " + titlenamech[i], 20,0,1,20,0,1);
      h2DCSV_45Jet[j][i] = new TH2F("h2DCSV_45Jet_" + namech[i] + "_" + namecut[j], "CSVv2 Discriminant for 5th and 6th Jets " + titlenamech[i], 20,0,1,20,0,1);
      h2DCSV_24Jet[j][i] = new TH2F("h2DCSV_24Jet_" + namech[i] + "_" + namecut[j], "CSVv2 Discriminant for 3rd and 5th Jets " + titlenamech[i], 20,0,1,20,0,1);
      h2DCSV_25Jet[j][i] = new TH2F("h2DCSV_25Jet_" + namech[i] + "_" + namecut[j], "CSVv2 Discriminant for 3rd and 6th Jets " + titlenamech[i], 20,0,1,20,0,1);
      h2DCSV_34Jet[j][i] = new TH2F("h2DCSV_34Jet_" + namech[i] + "_" + namecut[j], "CSVv2 Discriminant for 4th and 5th Jets " + titlenamech[i], 20,0,1,20,0,1);
      h2DCSV_35Jet[j][i] = new TH2F("h2DCSV_35Jet_" + namech[i] + "_" + namecut[j], "CSVv2 Discriminant for 4th and 6th Jets " + titlenamech[i], 20,0,1,20,0,1);

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
      
      
      TString jetn[6];
      jetn[0]= "Jet-0"; 
      jetn[1]= "Jet-1"; 
      jetn[2]= "Jet-2"; 
      jetn[3]= "Jet-3"; 
      jetn[4]= "Jet-4"; 
      jetn[5]= "Jet-5"; 
      
      for(int ij=0; ij<6; ij++){
	hCSV[ij][j][i] = new TH1F("hCSV_" + jetn[ij] + "_" + namech[i] + "_" + namecut[j],"CSV " + jetn[ij] + " " + titlenamech[i],10,0,1);
	hCSV[ij][j][i]->GetXaxis()->SetTitle("CSVv2");      
	hCvsL[ij][j][i] = new TH1F("hCvsL_" + jetn[ij] + "_" + namech[i] + "_" + namecut[j],"CvsL " + jetn[ij] + " " + titlenamech[i],10,0,1);
	hCvsL[ij][j][i]->GetXaxis()->SetTitle("CvsL");      
	hCvsB[ij][j][i] = new TH1F("hCvsB_" + jetn[ij] + "_" + namech[i] + "_" + namecut[j],"CvsB " + jetn[ij] + " " + titlenamech[i],10,0,1);
	hCvsB[ij][j][i]->GetXaxis()->SetTitle("CvsB");      
	hJetPt[ij][j][i] = new TH1F("hJetPt_" + jetn[ij] + "_" + namech[i] + "_" + namecut[j],"p_{T}^{Jet} " + jetn[ij] + " " + titlenamech[i],10,0,200);
	hJetPt[ij][j][i]->GetXaxis()->SetTitle("p_{T}[GeV]");      

	hJetpTUncVar[ij][j][i] = new TH1F("hJetpTUncVar_" + jetn[ij] + "_" + namech[i] + "_" + namecut[j], "#Delta pT^{Jet} " + jetn[ij] + " " + titlenamech[i], 20.0, 0.0, 2.0);
      }


      TString jetMassn[6][6];
      jetMassn[0][1] = "Jet01";
      jetMassn[0][2] = "Jet02";
      jetMassn[0][3] = "Jet03";
      jetMassn[0][4] = "Jet04";
      jetMassn[0][5] = "Jet05";
      jetMassn[1][2] = "Jet12";
      jetMassn[1][3] = "Jet13";
      jetMassn[1][4] = "Jet14";
      jetMassn[1][5] = "Jet15";
      jetMassn[2][3] = "Jet23";
      jetMassn[2][4] = "Jet24";
      jetMassn[2][5] = "Jet25";
      jetMassn[3][4] = "Jet34";
      jetMassn[3][5] = "Jet35";
      jetMassn[4][5] = "Jet45";

      for(int ja=0; ja<5; ja++){
	for(int jb=ja+1; jb<6; jb++){
	  hMassJet[ja][jb][j][i]    = new TH1F("hMassJet_" + jetMassn[ja][jb] + "_" + namech[i]+"_"+namecut[j],"transverse Mass of Dijets "+ jetMassn[ja][jb] + " " + titlenamech[i],80,0,400);
	}
      }
      hInvMassjj[j][i]  = new TH1F("hInvMassjj_" + namech[i]+"_"+namecut[j],"Compatible Inv. Mass " + titlenamech[i],80,0,400);
      hMaxMVAjj[j][i]   = new TH1F("hMaxMVAjj_"  + namech[i]+"_"+namecut[j],"Max. MVA (BDT) response " + titlenamech[i],20,-0.5,0.5);

      hEvtCatego[j][i]  = new TH1F("hEvtCatego_"+namech[i]+"_"+namecut[j],"ttbar Event Categorization " + titlenamech[i],4,-0.5,3.5);
      hEvtCatego[j][i]->GetXaxis()->SetTitle("ttbar");      

    }//for(i)
  }//for(j)
  

  hTJetPosition = new TH1F("hTJetPosition","CSV Position for jets from Top",  7,0,7);
  hWJetPosition = new TH1F("hWJetPosition","CSV Position for jets from W",    7,0,7);
  hOJetPosition = new TH1F("hOJetPosition","CSV Position for additional jets",7,0,7);
  
  h2DTJetPosition    = new TH2F("h2DTJetPosition","CSV Position for jets from Top Vs Dijet Rank",  7,0,7,7,0,7);
  h2DWJetPosition    = new TH2F("h2DWJetPosition","CSV Position for jets from W Vs Dijet Rank",    7,0,7,7,0,7);
  h2DOJetPositionMVA = new TH2F("h2DOJetPositionMVA","CSV Position for additional jets Vs MVA", 7,0,7,7,0,7);
  

  /************************
     MVA Training Trees
  *************************/
  // --- Write histograms
  TString outTFN = dirname + "/MVATrees" + hname + "_" + fname + ".root";
  TFile *TreeFile = new TFile(outTFN,"RECREATE" );

  std::vector<float> *Tvinput;
  Tvinput = new std::vector<float>;
  for (unsigned int nv = 0; nv < 10; nv++)  Tvinput->push_back(0);

  TTree *TMVASignal;
  TMVASignal = new TTree("TMVASignal","Signal Dijets");
  TMVASignal->Branch("pTjj",   &Tvinput->at(0),   "pTjj/F");
  TMVASignal->Branch("Mjj",    &Tvinput->at(1),    "Mjj/F");
  TMVASignal->Branch("DRjj",   &Tvinput->at(2),   "DRjj/F");
  TMVASignal->Branch("DPhijj", &Tvinput->at(3), "DPhijj/F");
  TMVASignal->Branch("CvsL1",  &Tvinput->at(4),  "CvsL1/F");
  TMVASignal->Branch("CvsL2",  &Tvinput->at(5),  "CvsL2/F");
  TMVASignal->Branch("CvsB1",  &Tvinput->at(6),  "CvsB1/F");
  TMVASignal->Branch("CvsB2",  &Tvinput->at(7),  "CvsB2/F");
  TMVASignal->Branch("CSV1",   &Tvinput->at(8),   "CSV1/F");
  TMVASignal->Branch("CSV2",   &Tvinput->at(9),   "CSV2/F");

  TTree *TMVABkg;
  TMVABkg = new TTree("TMVABkg","Background Dijets");
  TMVABkg->Branch("pTjj",   &Tvinput->at(0),   "pTjj/F");
  TMVABkg->Branch("Mjj",    &Tvinput->at(1),    "Mjj/F");
  TMVABkg->Branch("DRjj",   &Tvinput->at(2),   "DRjj/F");
  TMVABkg->Branch("DPhijj", &Tvinput->at(3), "DPhijj/F");
  TMVABkg->Branch("CvsL1",  &Tvinput->at(4),  "CvsL1/F");
  TMVABkg->Branch("CvsL2",  &Tvinput->at(5),  "CvsL2/F");
  TMVABkg->Branch("CvsB1",  &Tvinput->at(6),  "CvsB1/F");
  TMVABkg->Branch("CvsB2",  &Tvinput->at(7),  "CvsB2/F");
  TMVABkg->Branch("CSV1",   &Tvinput->at(8),   "CSV1/F");
  TMVABkg->Branch("CSV2",   &Tvinput->at(9),   "CSV2/F");


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
  //float CSV_WP = 0.935; // Tight
  int btagSysPar = 0;

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
  NormWeight = SFLumi(fname, 2170, nNorm_Event);  

  std::cout << "-----------------------                                 -------------------------" << std::endl;
  std::cout << "Number of Events     = " << nNorm_Event << std::endl;
  std::cout << "Normalization Factor = " << NormWeight  << std::endl;
  std::cout << "---------------------------------------------------------------------------------" << std::endl;

  /************************
      MVA: BDT branches
  *************************/  
  std::vector<std::string> inputVariableNames;
  inputVariableNames.push_back("pTjj");
  inputVariableNames.push_back("Mjj");
  inputVariableNames.push_back("DRjj");
  inputVariableNames.push_back("DPhijj");    
  inputVariableNames.push_back("CvsL1");    
  inputVariableNames.push_back("CvsL2");    
  inputVariableNames.push_back("CvsB1");    
  inputVariableNames.push_back("CvsB2");    
  inputVariableNames.push_back("CSV1");    
  inputVariableNames.push_back("CSV2");    
  IClassifierReader *classReader = new ReadBDT(inputVariableNames);



  /********************************
             Event Loop
  ********************************/
  std::cout << "--- Processing: " << theTree.GetEntries() << " events" << std::endl;
  
  for (Long64_t ievt=0; ievt<theTree.GetEntries();ievt++) {
    // for (Long64_t ievt=0; ievt<10000;ievt++) {
    
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
    if (_syst && syst_varname.Contains("ScaleR"))
      PUWeight = PUWeight*(*ScaleWeight)[scaleSysPar];
    
    int NJets,NBtagJets;
    
    TLorentzVector Lep;
    
    Lep.SetPxPyPzE(Lep_px,Lep_py,Lep_pz,Lep_E);
    if (_syst && syst_varname.Contains("LES")){
      if (syst_varname.Contains("Up"))     Lep.SetPxPyPzE(Lep_px+(Lep_px*Lep_LES),
							  Lep_py+(Lep_py*Lep_LES),
							  Lep_pz+(Lep_pz*Lep_LES),
							  Lep_E);
      if (syst_varname.Contains("Down"))   Lep.SetPxPyPzE(Lep_px-(Lep_px*Lep_LES),
							  Lep_py-(Lep_py*Lep_LES),
							  Lep_pz-(Lep_pz*Lep_LES),
							  Lep_E);
    }
    if(Lep.Pt() < 30)  continue; // Lep pT >30GeV
    
    // Transverse W Mass
    TLorentzVector METv;
    METv.SetPtEtaPhiE(MET,0.0,MET_Phi,MET);
    
    float mT = sqrt(2*MET*Lep.Pt()*(1.0-cos( Lep.DeltaPhi(METv) )));
    
    // Jets     
    NJets      = 0;
    NBtagJets  = 0;
    
    // Global SF_b-tag
    // From: https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration
    float btagUnc_val = 0.0;

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
    
    std::vector<ComJet> Jets;
    

    for(int ijet=0; ijet < Jet_px->size(); ijet++){

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
      
      ComJet jet;
      jet.SetPxPyPzE(JetSystVar * (*Jet_px)[ijet],
		     JetSystVar * (*Jet_py)[ijet],
		     (*Jet_pz)[ijet],
		     (*Jet_E)[ijet]);
      jet.Flavour = (*Jet_partonFlavour)[ijet];
      jet.CSV = (*Jet_CSV)[ijet];
      jet.CvsL = (*Jet_CvsL)[ijet];
      jet.CvsB = (*Jet_CvsB)[ijet];
      jet.Mom = -1;
      if((fname.Contains("ttbar")  && !fname.Contains("Bkg")) && 
	 (*Jet_GENmatched)[ijet] != -999) jet.Mom = (*GenJet_Mom)[(*Jet_GENmatched)[ijet]];

      if(jet.Pt() > 25){ // Jet pT Cut
	//std::cout << "Evt = " << Event << " RECO-pT = " << jet.Pt() << " RECO-Flv = " << jet.Flavour << " GEN-Match = " << (*Jet_GENmatched)[ijet] << " GEN-Mom = " << jet.Mom << std::endl; 		
	Jets.push_back(jet);

	/*******************************************
                       b-tagging
	*******************************************/    
	// OLD b-tagging Method: 
	// b-tagging WP from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagging
	// Not recommended for analysis with observables that depends ontagging.
	// if (fname.Contains("Data")) btagDisc = fBTagSF->IsTagged((*Jet_CSV)[ijet], -999999, jet.Pt(), jet.Eta());
	// else btagDisc = fBTagSF->IsTagged((*Jet_CSV)[ijet], JetFlav, jet.Pt(), jet.Eta()); 
	// New Method (Event SF from tth group)
	// https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration
	if(jet.CSV > CSV_WP) NBtagJets++; // Number of b-tagged jets
		
      } // if(Jet_pT)
    }// for(jets)
    
    NJets = Jets.size();

    float Mjj    = 0.0; // Dijet Inv. Mass   
    float MaxMVA = 0.0; // Maximum MVA   
    
    // Only Final Cut Level
    if(NJets > 5 && NBtagJets > 1){ 
      // Estimation of the DiJet invariant mass closest to the W mass
      bool ReArrangeMjj = false; // Jets keep the same order
      Mjj = DiJetMassCorrection(Jets, ReArrangeMjj);
      
      /*******************************************
                      TMVA:BDT
      *******************************************/
      std::vector<std::vector<float>> Dijet_Var;
      // 0: Dijet Tagger (To produce the training tree) 
      // 1: Jet0 Index 
      // 2: Jet1 Index 
      // 3: MVA Response 
      // 4: Dijet_M;
      // 5: Dijet_pT;
      // 6: Dijet_DPhi;
      // 7: Dijet_DR;
      bool MVAEstimation = true;
      Dijet_Var = DiJetVariableEstimation(Jets, DRAddJets, MVAEstimation, classReader);

      // Estimation of the BEST Dijet using MVA Response
      std::vector<int> AddDiJetMVAIndex;
      // 0: Jet0  Index
      // 1: Jet1  Index
      // 2: DiJet Index
      AddDiJetMVAIndex  = MVAResponse (Dijet_Var);
      MaxMVA = Dijet_Var.at(3).at(2);

      
      //----------------------------------------
      // MVA Arrange
      if(AddDiJetMVAIndex[0] > 1 && AddDiJetMVAIndex[1] > 1){ // Rearrange jet positions > 1	
      ComJet Jet0 = Jets.at(2);
      ComJet Jet1 = Jets.at(3);
      
      Jets.at(2) = Jets.at(AddDiJetMVAIndex[0]);      
      Jets.at(3) = Jets.at(AddDiJetMVAIndex[1]);      
      Jets.at(AddDiJetMVAIndex[0]) = Jet0;
      Jets.at(AddDiJetMVAIndex[1]) = Jet1;
      
      } // if(AddDiJetMVAIndex[JetPos] > 1)
      
      
      if (fname.Contains("ttbar")  && !fname.Contains("Bkg")){
	/*******************************************
             TEMPORAL -- ttjj Studie: Mom match 
	*******************************************/
	int TJet = 0, WJet = 0, AddJet = 0;
	
	if(NJets == 6){
	  for (unsigned int ij = 0; ij < NJets; ij ++){
	    ComJet jetcand = Jets[ij]; 
	    if (jetcand.Mom == 6   && std::abs(jetcand.Flavour) == 5) TJet++;
	    else if (jetcand.Mom == 24  && std::abs(jetcand.Flavour) != 5) WJet++;
	    else if (jetcand.Mom != 24  && jetcand.Mom != 6) AddJet++;
	  } // for(ij)
	}// if(NJets > 5)
	
	//if (!(TJet == 2 && WJet == 2 && AddJet == 2)) continue; // Jump Event 

	// --WARNING!!
	// W Jets to the end of the Jet Vector    
	// ComJet Wjj_a = Jets[WIndex[0]];
	// ComJet Wjj_b = Jets[WIndex[1]];
	// Jets.erase(Jets.begin() + WIndex[1]);
	// Jets.erase(Jets.begin() + WIndex[0]);
	// Jets.push_back(Wjj_a);
	// Jets.push_back(Wjj_b);
	/*******************************************
                       TEMPORAL END
	*******************************************/
      
	//----------------------------------------
	// Fill Training Trees: Signal (Add DiJet) and Bkg (W Dijet)
	// MVATrainingTree(TMVASignal, TMVABkg, Dijet_Var, Tvinput);
	
	// Estimation of the W and Top Dijet using GEN info
	// 0, 1 and 2: WJet    Index
	// 3, 4 and 5: TopJet  Index
	std::vector<int> WTopDiJetIndex = WTopTagger  (Dijet_Var);
	
	//----------------------------------------
	// Histograms only filled for ttbar samples
	std::vector<int> WIndex;    
	std::vector<int> TopIndex;    
	std::vector<int> AddIndex;    
	
	for (unsigned int ij = 0; ij < NJets; ij ++){
	  ComJet jetcand = Jets[ij]; 
	  if      (jetcand.Mom == 6){
	    hTJetPosition->Fill(ij);
	    TopIndex.push_back(ij);
	  }      
	  else if (jetcand.Mom == 24){
	    hWJetPosition->Fill(ij);
	    WIndex.push_back(ij);
	  }      
	  else{
	    hOJetPosition->Fill(ij);      
	    AddIndex.push_back(ij);
	  }    
	} // for(ij)
	
	//----------------------------------------
	// Re-arrange Jets using GEN info
	// GEN Arrange
	// std::vector<ComJet> arrJets;
	// arrJets.resize(Jets.size());
	
	// arrJets.at(0) = Jets.at(TopIndex[0]);
	// arrJets.at(1) = Jets.at(TopIndex[1]);
	// arrJets.at(2) = Jets.at(AddIndex[0]);
	// arrJets.at(3) = Jets.at(AddIndex[1]);
	// arrJets.at(4) = Jets.at(WIndex[0]);
	// arrJets.at(5) = Jets.at(WIndex[1]);
	
	// Jets.clear();
	// Jets = arrJets;
	
	// Filled ONLY if the event has the 3 components: (2W, 2T 2Add)
	// h2DWJetPosition->Fill(WIndex[0], WTopDiJetIndex[0]);
	// h2DWJetPosition->Fill(WIndex[1], WTopDiJetIndex[1]);
	
	// h2DTJetPosition->Fill(TopIndex[0], WTopDiJetIndex[3]);
	// h2DTJetPosition->Fill(TopIndex[1], WTopDiJetIndex[4]);
	
	// h2DOJetPositionMVA->Fill(AddIndex[0], AddDiJetMVAIndex[0]);
	// h2DOJetPositionMVA->Fill(AddIndex[1], AddDiJetMVAIndex[1]);
	
      }// if (ttbar)
      
    } //if (6jets && 1btag)    
    
    
    /*******************************************
                Number of GENJets
    *******************************************/    
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
      // SFIDISOTrigger(SF_ID_ISO_Tr,
      // 		     Lep, Channel,
      // 		     hmuIDISOSF, hmuTriggerSF,
      // 		     heIDISOSF,  heTriggerSF);
      
      // Easiest adaptation: Create Lep_SF in the same way the Lep_SF vector
      SF_ID_ISO_Tr = (*Lep_SF);
      
      if(_syst && syst_varname.Contains("LepSF")){
	float SFSystUnc = SF_ID_ISO_Tr[0]*0.015;
	if(syst_varname.Contains("Up"))   PUWeight = PUWeight * (SF_ID_ISO_Tr[1] + SFSystUnc) ;
	if(syst_varname.Contains("Down")) PUWeight = PUWeight * (SF_ID_ISO_Tr[2] - SFSystUnc);
      }
      else PUWeight = PUWeight * SF_ID_ISO_Tr[0]; // Central

      // if(_idiso_unc){
      // 	if     (IDISOUnc == "Up")   PUWeight = PUWeight * (SF_ID_ISO_Tr[1] + SF_ID_ISO_Tr[2]);
      // 	else if(IDISOUnc == "Down") PUWeight = PUWeight * (SF_ID_ISO_Tr[1] - SF_ID_ISO_Tr[2]);
      // 	else if(IDISOUnc == "Nom")  PUWeight = PUWeight * (SF_ID_ISO_Tr[1]);
      // } // if(_idiso_unc)
      
      // else if(_tr_unc){
      // 	if     (TrUnc=="Up")   PUWeight=PUWeight*(SF_ID_ISO_Tr[3] + SF_ID_ISO_Tr[4]);
      // 	else if(TrUnc=="Down") PUWeight=PUWeight*(SF_ID_ISO_Tr[3] - SF_ID_ISO_Tr[4]);
      // 	else if(TrUnc=="Nom")  PUWeight=PUWeight*(SF_ID_ISO_Tr[3]);	
      // }// if(_tr_unc) 
      
      // // Check the SF_e implementation
      // // If the electron is in the transition region, SF = 1
      // else  if(SF_ID_ISO_Tr[0] != 0.0)  PUWeight=PUWeight*(SF_ID_ISO_Tr[0]); 
      
    }// else(Contain("Data"))
    

    /***************************
            Selection
    ***************************/
    
    int                            cut = 0; // Single Lepton (from Tree)
    if(NJets > 5)                  cut = 1; // + 6 Jets 
    if(NJets > 5 && NBtagJets > 1) cut = 2; // + 2 b-tag
    if(NJets > 5 && NBtagJets > 2) cut = 3; // + 3 b-tag
    
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
      int cone_NaddcJets = (*GenConeCat)[7];
      
      bool Isttjj = false;
      bool Isttbb = false;
      bool Isttcc = false;
      bool Isttb  = false;
      bool IsttLF = false;
      bool Istt   = false;
      
      if(cone_NbJets > 1 && cone_NJets > 5) Isttjj = true;
      
      // Categorization based in the Visible Ph-Sp
      if      (cone_NbJets > 3  && cone_NJets > 5) Isttbb = true;
      else if (cone_NbJets > 2  && cone_NJets > 5) Isttb  = true;
      else if (cone_NbJets > 1  && cone_NJets > 5 && cone_NcJets > 1) Isttcc = true;
      else if (cone_NbJets > 1  && cone_NJets > 5) IsttLF = true;
      else Istt = true;

      // Categorization based in the Full Ph-Sp
      // if      (cone_NaddJets > 1) Isttjj = true;

      // if      (cone_NaddbJets > 1) Isttbb = true;
      // else if (cone_NaddbJets > 0) Isttb  = true;
      // else if (cone_NaddcJets > 1) Isttcc = true;
      // else if (cone_NaddJets  > 0) IsttLF = true;
      // else Istt = true;


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
      //hSFIDISO[icut][Channel]->Fill(SF_ID_ISO_Tr[1],PUWeight);
      //hSFIDISOError[icut][Channel]->Fill(SF_ID_ISO_Tr[2],PUWeight);
      //hSFTrigger[icut][Channel]->Fill(SF_ID_ISO_Tr[3],PUWeight);
      //hSFTriggerError[icut][Channel]->Fill(SF_ID_ISO_Tr[4],PUWeight);
    
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
      // Dijet InvMass
      hInvMassjj[icut][Channel]->Fill(Mjj, PUWeight);
      // Dijet Max MVA (BDT) Response
      hMaxMVAjj[icut][Channel]->Fill(MaxMVA, PUWeight);
      // Global btag SF
      h2DSFbtag_Global[icut][Channel]->Fill((*Jet_SF_CSV)[btagUnc::CENTRAL], btagUnc_val, PUWeight);
      hSFbtag_Global[icut][Channel]->Fill((*Jet_SF_CSV)[btagUnc::CENTRAL], PUWeight);
      hSFbtag_Global_var[icut][Channel]->Fill(btagUnc_val, PUWeight);
      
      // 2D CSV discriminant plot for all Jets
      if(Jets.size() > 3)  h2DCSV_23Jet[icut][Channel]->Fill(Jets[2].CSV, Jets[3].CSV, PUWeight);
      if(Jets.size() > 4)  h2DCSV_24Jet[icut][Channel]->Fill(Jets[2].CSV, Jets[4].CSV, PUWeight);
      if(Jets.size() > 4)  h2DCSV_34Jet[icut][Channel]->Fill(Jets[3].CSV, Jets[4].CSV, PUWeight);
      if(Jets.size() > 5)  h2DCSV_25Jet[icut][Channel]->Fill(Jets[2].CSV, Jets[5].CSV, PUWeight);
      if(Jets.size() > 5)  h2DCSV_35Jet[icut][Channel]->Fill(Jets[3].CSV, Jets[5].CSV, PUWeight);

      for(int ijet=0; ijet < Jets.size(); ijet++){

	ComJet jet = Jets[ijet];
	
	if (ijet < 6){	  
	  hJetPt[ijet][icut][Channel]->Fill(jet.Pt(), PUWeight);
	  hCSV  [ijet][icut][Channel]->Fill(jet.CSV,  PUWeight);
	  hCvsL [ijet][icut][Channel]->Fill(jet.CvsL, PUWeight);
	  hCvsB [ijet][icut][Channel]->Fill(jet.CvsB, PUWeight);
	}

	if(jet.Flavour == 5){
	  h2DSFbtag_b[icut][Channel]->Fill(jet.Pt(), fabs(jet.Eta()), PUWeight); // b-Flavour
	  if(jet.CSV > CSV_WP) h2DSFbtag_btag_b[icut][Channel]->Fill(jet.Pt(), fabs(jet.Eta()), PUWeight);  
	}	
	else if(jet.Flavour == 4){
	  h2DSFbtag_c[icut][Channel]->Fill(jet.Pt(), fabs(jet.Eta()), PUWeight); // c-Flavour
	  if(jet.CSV > CSV_WP) h2DSFbtag_btag_c[icut][Channel]->Fill(jet.Pt(), fabs(jet.Eta()), PUWeight);  
	}
	else if(jet.Flavour == 1 || jet.Flavour == 2 || jet.Flavour == 3){
	  h2DSFbtag_l[icut][Channel]->Fill(jet.Pt(), fabs(jet.Eta()), PUWeight); // l-Flavour
	  if(jet.CSV > CSV_WP) h2DSFbtag_btag_l[icut][Channel]->Fill(jet.Pt(), fabs(jet.Eta()), PUWeight);  
	}
	
	//Dijet Invariant Mass 
	int jbmax = std::min(6,NJets);
	for(int jjet=ijet+1; jjet < jbmax; jjet++){
	  ComJet jet_ = Jets[jjet];
	  float DijetInvMass = (jet+jet_).M(); 
	  hMassJet[ijet][jjet][icut][Channel]->Fill(DijetInvMass, PUWeight);
	}// for(jjet)
	
      }//for(ijet)     

    }//for(icuts)     
    
    Jets.clear();

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
  std::cout << std::endl;

  //Acceptance-Efficiency
  TH1F *Yields;
  Yields = new TH1F("Yields", "Yields",12,0,12);
  int nbin = 1;
  
  for(int nc = 0; nc < 4; nc++){
    AccEvent[nc][2] = AccEvent[nc][0] + AccEvent[nc][1];
    EffEvent[nc][2] = EffEvent[nc][0] + EffEvent[nc][1];

    for(int nch = 0; nch < 3; nch++){      
      
      float EffError;
      if (AccEvent[nc][nch] != 0.0) EffError = sqrt(AccEvent[nc][nch])*EffEvent[nc][nch]/AccEvent[nc][nch];
      else EffError = 0.0;

      Yields->GetXaxis()->SetBinLabel(nbin, namecut[nc] + " " + namech[nch]);
      Yields->SetBinContent(nbin, EffEvent[nc][nch]);
      Yields->SetBinError  (nbin, EffError);
      nbin++;
      
      std::cout << "-- Acceptace  " << namecut[nc] << " " << namech[nch] << ": " << AccEvent[nc][nch] << std::endl;
      std::cout << "-- Efficiency " << namecut[nc] << " " << namech[nch] << ": " << EffEvent[nc][nch] << " +/- " << EffError << std::endl;
      std::cout << std::endl;
    }
    std::cout << "-----------------------------" << std::endl;
  }
  

  TreeFile->cd();    

  TMVASignal->Write();
  TMVABkg   ->Write();

  TreeFile->Close();    

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
  
  // --- Write histograms
  TString outfname=dirname + "/hSF-" + hname + "_" + fname + ".root";
  TFile *target  = new TFile(outfname,"RECREATE" );  
  
  target->cd();
  
  Yields->Write();

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
      h2DCSV_45Jet[j][i]->Write();
      h2DCSV_24Jet[j][i]->Write();
      h2DCSV_25Jet[j][i]->Write();
      h2DCSV_34Jet[j][i]->Write();
      h2DCSV_35Jet[j][i]->Write();
      
      h2DSFbtag_Global[j][i]->Write();
      hSFbtag_Global[j][i]->Write();
      hSFbtag_Global_var[j][i]->Write();

      for(int ij=0; ij<6; ij++){
	hJetPt[ij][j][i]->Write();
	hCSV  [ij][j][i]->Write();
	hCvsL [ij][j][i]->Write();
	hCvsB [ij][j][i]->Write();
      }

      for(int ja=0; ja<5; ja++){
	for(int jb=ja+1; jb<6; jb++){
	  hMassJet[ja][jb][j][i]->Write();
	}
      }      
      hInvMassjj[j][i]->Write();
      hMaxMVAjj [j][i]->Write();

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

  hTJetPosition->Write();
  hWJetPosition->Write();
  hOJetPosition->Write();

  h2DWJetPosition->Write();
  h2DTJetPosition->Write();
  h2DOJetPositionMVA->Write();
    
  target->Close();

  std::cout << "File saved as " << outfname << std::endl;

  delete classReader;
  delete Tvinput;
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

void print_progress(int TreeEntries, Long64_t ievt){
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


float DiJetMassCorrection(std::vector<ComJet>  &Jets, bool ReArrange = false){

  /*******************************************
            Dijet Invariant Mass 
  *******************************************/
  int InvMassIndex[2];
  float minDeltaMjj = 9999.;
  float Mjj = 0.0;
  
  for(int ijet=0; ijet < Jets.size(); ijet++){
    
    ComJet jet_i = Jets[ijet];
    
    for(int jjet=ijet+1; jjet < Jets.size(); jjet++){
      ComJet jet_j = Jets[jjet];
      
      float DijetInvMass = (jet_i+jet_j).M(); 
      float DeltaMjj = std::abs(DijetInvMass-80.3);
      
      if(minDeltaMjj > DeltaMjj){
	minDeltaMjj = DeltaMjj;
	Mjj = DijetInvMass;
	InvMassIndex[0] = ijet;
	InvMassIndex[1] = jjet;
      }
      
    }// for(jjets)

  }// for(ijet)

  // Change of Jet Ordering
  if(ReArrange && (InvMassIndex[0] > 1 || InvMassIndex[1] > 1)){    
    ComJet Mjj_a = Jets[InvMassIndex[0]];
    ComJet Mjj_b = Jets[InvMassIndex[1]];
  
    // Delete the Mjj Jets: First the last one!
    Jets.erase(Jets.begin()+InvMassIndex[1]);
    Jets.erase(Jets.begin()+InvMassIndex[0]);

    // Move them to the end
    Jets.push_back(Mjj_a);
    Jets.push_back(Mjj_b);
  } // if(ReArrange)

  return Mjj;
}


std::vector<std::vector<float>> DiJetVariableEstimation(std::vector<ComJet>  &Jets, float DRAddJets, bool MVA, IClassifierReader *classReader){

  std::vector<std::vector<float>> DijetOutput;
  std::vector<float> Dijet_Signal; 
  std::vector<float> Dijet_Index0; 
  std::vector<float> Dijet_Index1; 
  std::vector<float> Dijet_MVARes; 
  std::vector<float> Dijet_M; 
  std::vector<float> Dijet_pT; 
  std::vector<float> Dijet_DR; 
  std::vector<float> Dijet_DPhi; 
  std::vector<float> Dijet_CvsL1; 
  std::vector<float> Dijet_CvsL2; 
  std::vector<float> Dijet_CvsLjj; 
  std::vector<float> Dijet_CvsB1; 
  std::vector<float> Dijet_CvsB2; 
  std::vector<float> Dijet_CvsBjj; 
  std::vector<float> Dijet_CSV1; 
  std::vector<float> Dijet_CSV2; 
  
  float MinDelta_DR = 999.;

  int AddDijetIndex     =  0;
  int bestAddDiJetIndex = -1;
  int WDiJetIndex       = -1;
  int TopDiJetIndex     = -1;

  for(int ijet=0; ijet < Jets.size(); ijet++){
    ComJet jet_i = Jets[ijet];
    
    for(int jjet=ijet+1; jjet < Jets.size(); jjet++){
      ComJet jet_j = Jets[jjet];
      
      float vMjj    = (jet_i+jet_j).M();
      float vpTjj   = (jet_i+jet_j).Pt();
      float vDRjj   = jet_i.DeltaR(jet_j);
      float vDPhijj = jet_i.DeltaPhi(jet_j);
      float vCvsL1  = jet_i.CvsL;
      float vCvsL2  = jet_j.CvsL;
      float vCvsLjj = jet_i.CvsL + jet_j.CvsL;
      float vCvsB1  = jet_i.CvsB;
      float vCvsB2  = jet_j.CvsB;
      float vCvsBjj = jet_i.CvsB + jet_j.CvsB;
      float vCSV1   = jet_i.CSV;
      float vCSV2   = jet_j.CSV;

      // Kinematic Variables      
      Dijet_M   .push_back(vMjj); 
      Dijet_pT  .push_back(vpTjj); 
      Dijet_DR  .push_back(vDRjj); 
      Dijet_DPhi.push_back(vDPhijj); 
      Dijet_CvsL1.push_back(vCvsL1); 
      Dijet_CvsL2.push_back(vCvsL2); 
      Dijet_CvsLjj.push_back(vCvsLjj); 
      Dijet_CvsB1.push_back(vCvsB1); 
      Dijet_CvsB2.push_back(vCvsB2); 
      Dijet_CvsBjj.push_back(vCvsBjj); 
      Dijet_CSV1.push_back(vCSV1); 
      Dijet_CSV2.push_back(vCSV2); 
      // Signal/Bkg tagger (Only for trainning)
      Dijet_Signal.push_back(0.0);
      // Jet Index
      Dijet_Index0.push_back(ijet);
      Dijet_Index1.push_back(jjet);

      // Is it a signal (additional) pair?
      if (jet_i.Mom != 6  && jet_j.Mom != 6 && 
	  jet_i.Mom != 24 && jet_j.Mom != 24 ) {
	
	float Delta_DR = std::abs(jet_i.DeltaR(jet_j) - DRAddJets);

	if(Delta_DR < MinDelta_DR){
	  MinDelta_DR = Delta_DR;	  
	  bestAddDiJetIndex = AddDijetIndex;	  
	} // if (Delta_DR)
      
      } // if(Mom)
      
      // Is it a W pair?
      if (jet_i.Mom == 24 && jet_j.Mom == 24 ) WDiJetIndex = AddDijetIndex;
      // Is it a Top pair?
      if (jet_i.Mom == 6 && jet_j.Mom == 6 ) TopDiJetIndex = AddDijetIndex;
      
      // TMVA Response    
      float vMVARes;
      if (MVA){
	std::vector<double> inputVariableValues;
	inputVariableValues.push_back(vMjj);
	inputVariableValues.push_back(vpTjj);
	inputVariableValues.push_back(vDRjj);
	inputVariableValues.push_back(vDPhijj);
	inputVariableValues.push_back(vCvsL1);
	inputVariableValues.push_back(vCvsL2);
	inputVariableValues.push_back(vCvsB1);
	inputVariableValues.push_back(vCvsB2);
	inputVariableValues.push_back(vCSV1);
	inputVariableValues.push_back(vCSV2);

	vMVARes = classReader->GetMvaValue(inputVariableValues);
      }
      else  vMVARes = 0.0;
      Dijet_MVARes.push_back(vMVARes);

      // Global DiJet Index
      AddDijetIndex ++;      

    }// for(jjets)    
  }// for(ijet)
  
  // Signal/Bkg tagger (Only for trainning)
  if(MinDelta_DR < 999.) Dijet_Signal.at(bestAddDiJetIndex) = 1.0;
  if(WDiJetIndex > -1)   Dijet_Signal.at(WDiJetIndex)       = 2.0;
  if(TopDiJetIndex > -1) Dijet_Signal.at(TopDiJetIndex)     = 3.0;
  DijetOutput.push_back(Dijet_Signal);
  // Jet Index
  DijetOutput.push_back(Dijet_Index0);
  DijetOutput.push_back(Dijet_Index1);
  // MVA Response
  DijetOutput.push_back(Dijet_MVARes);
  // Kinematic Variables
  DijetOutput.push_back(Dijet_M);
  DijetOutput.push_back(Dijet_pT);
  DijetOutput.push_back(Dijet_DR);
  DijetOutput.push_back(Dijet_DPhi);
  DijetOutput.push_back(Dijet_CvsL1);
  DijetOutput.push_back(Dijet_CvsL2);
  DijetOutput.push_back(Dijet_CvsLjj);
  DijetOutput.push_back(Dijet_CvsB1);
  DijetOutput.push_back(Dijet_CvsB2);
  DijetOutput.push_back(Dijet_CvsBjj);
  DijetOutput.push_back(Dijet_CSV1);
  DijetOutput.push_back(Dijet_CSV2);

  return DijetOutput;
  // 0: Dijet Tagger (To produce the training tree) 
  // 1: Jet0 Index 
  // 2: Jet1 Index 
  // 3: MVA Response 
  // 4: Dijet_M;
  // 5: Dijet_pT;
  // 6: Dijet_DPhi;
  // 7: Dijet_DR;
  
}

  
#endif

