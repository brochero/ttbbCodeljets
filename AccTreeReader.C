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
#include <ttbar_category.h> // Event Categorization
  
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
  const char * _ttbar_id = "";

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
  std::vector<float> *Jet_px=0, *Jet_py=0, *Jet_pz=0, *Jet_pT=0, *Jet_E=0;
  std::vector<int> *Jet_Mom=0;

  std::vector<int>   *Jet_partonFlavour=0;

  // Categorization
  int  GenCat_ID;
  std::vector<int> *GenConeCat=0; 
 /*********************************
           Tree Branches
  **********************************/
  
  theTree.SetBranchAddress( "genweight",     &GENWeight  );
  theTree.SetBranchAddress( "scaleweight",   &ScaleWeight);
  theTree.SetBranchAddress( "genchannel",    &Channel    );
  theTree.SetBranchAddress( "genhiggscatid", &GenCat_ID  );
  theTree.SetBranchAddress( "genconecatid",  &GenConeCat );


  theTree.SetBranchAddress( "genlepton_pT",  &Lep_pT );
  theTree.SetBranchAddress( "genlepton_eta", &Lep_eta);
  theTree.SetBranchAddress( "genjet_px",     &Jet_px );
  theTree.SetBranchAddress( "genjet_py",     &Jet_py );
  theTree.SetBranchAddress( "genjet_pz",     &Jet_pz );
  theTree.SetBranchAddress( "genjet_E",      &Jet_E  );

  theTree.SetBranchAddress( "genjet_mom",     &Jet_Mom);


  /*********************************
         Output Tree: MVA     
  **********************************/
  float b_pTjj, b_Mjj, b_DRjj, b_DPhijj;    

  TTree *MVASignaltree;
  MVASignaltree = new TTree("MVASignaltree","Dijets from W decay");
  MVASignaltree->Branch("pTjj",   &b_pTjj,   "pTjj/F");
  MVASignaltree->Branch("Mjj",    &b_Mjj,    "Mjj/F");
  MVASignaltree->Branch("DRjj",   &b_DRjj,   "DRjj/F");
  MVASignaltree->Branch("DPhijj", &b_DPhijj, "DPhijj/F");

  TTree *MVAtree;
  MVAtree = new TTree("MVAtree","All Dijets");
  MVAtree->Branch("pTjj",   &b_pTjj,   "pTjj/F");
  MVAtree->Branch("Mjj",    &b_Mjj,    "Mjj/F");
  MVAtree->Branch("DRjj",   &b_DRjj,   "DRjj/F");
  MVAtree->Branch("DPhijj", &b_DPhijj, "DPhijj/F");

  TTree *MVABkgtree;
  MVABkgtree = new TTree("MVABkgtree","All Dijets but W decay");
  MVABkgtree->Branch("pTjj",   &b_pTjj,   "pTjj/F");
  MVABkgtree->Branch("Mjj",    &b_Mjj,    "Mjj/F");
  MVABkgtree->Branch("DRjj",   &b_DRjj,   "DRjj/F");
  MVABkgtree->Branch("DPhijj", &b_DPhijj, "DPhijj/F");

  /*********************************
             Histograms
  **********************************/

  //Correct Statistical Uncertainty Treatment
  //TH1::SetDefaultSumw2(kTRUE);  
  
  TH1F *hNJets[2];
  TH1F *hJetMatch[2];
  TH2F *h2DJetMjjID_I[2], *h2DJetMjjID_II[2];

  TH1F *hJetPt[6][2], *hOJetPt[4][2], *hWJetPt[2][2];

  TH1F *hpTJet [5][6][2], *hMassJet [5][6][2], *hDRJet [5][6][2], *hDPhiJet [5][6][2];
  TH1F *hOpTJet[3][4][2], *hOMassJet[3][4][2], *hODRJet[3][4][2], *hODPhiJet[3][4][2];

  TH1F *hpTWjj[2],*hInvMassWjj[2],*hDRWjj[2], *hDPhiWjj[2];
  TH1F *hpTOjj[2],*hInvMassOjj[2],*hDROjj[2], *hDPhiOjj[2];

  TH1F *hpTminjj [2], *hpTmaxjj [2], *hInvMassjj [2], *hDRminjj [2], *hDRmaxjj [2], *hDPhiminjj [2], *hDPhimaxjj [2];
  TH1F *hOpTminjj[2], *hOpTmaxjj[2], *hOInvMassjj[2], *hODRminjj[2], *hODRmaxjj[2], *hODPhiminjj[2], *hODPhimaxjj[2];

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
        

    hJetMatch[i] = new TH1F("hJetMatch_" + namech[i], "W Jet Match " + titlenamech[i],13,0,13);
    hJetMatch[i]->GetXaxis()->SetBinLabel(1,"Total # Evt");
    hJetMatch[i]->GetXaxis()->SetBinLabel(2,"#Delta Mjj is Wjj");
    hJetMatch[i]->GetXaxis()->SetBinLabel(3,"Min #Delta R_{jj} is Wjj");
    hJetMatch[i]->GetXaxis()->SetBinLabel(4,"Max #Delta R_{jj} is Wjj");
    hJetMatch[i]->GetXaxis()->SetBinLabel(5,"Min #Delta #Phi_{jj} is Wjj");
    hJetMatch[i]->GetXaxis()->SetBinLabel(6,"Max #Delta #Phi_{jj} is Wjj");
    hJetMatch[i]->GetXaxis()->SetBinLabel(7,"At least 1 ->");
    hJetMatch[i]->GetXaxis()->SetBinLabel(8, "#Delta Mjj is Wjj");
    hJetMatch[i]->GetXaxis()->SetBinLabel(9, "Min #Delta R_{jj} is Wjj");
    hJetMatch[i]->GetXaxis()->SetBinLabel(10,"Max #Delta R_{jj} is Wjj");
    hJetMatch[i]->GetXaxis()->SetBinLabel(11,"Min #Delta #Phi_{jj} is Wjj");
    hJetMatch[i]->GetXaxis()->SetBinLabel(12,"Max #Delta #Phi_{jj} is Wjj");
    hJetMatch[i]->GetXaxis()->SetBinLabel(13,"#Delta Mjj && Min #Delta R_{jj} is Wjj");

    h2DJetMjjID_I[i] = new TH2F("h2DJetMjjID1_" + namech[i], "Index ID " + titlenamech[i], 6, 0, 6, 6, 0, 6);
    h2DJetMjjID_I[i]->GetXaxis()->SetTitle("W Index");      
    h2DJetMjjID_I[i]->GetYaxis()->SetTitle("Mjj Index");      

    h2DJetMjjID_II[i] = new TH2F("h2DJetMjjID2_" + namech[i], "Index ID " + titlenamech[i],6, 0, 6, 6, 0, 6);
    h2DJetMjjID_II[i]->GetXaxis()->SetTitle("W Index");      
    h2DJetMjjID_II[i]->GetYaxis()->SetTitle("Mjj Index");      


    TString jetn[6];
    jetn[0]= "Jet-0"; 
    jetn[1]= "Jet-1"; 
    jetn[2]= "Jet-2"; 
    jetn[3]= "Jet-3"; 
    jetn[4]= "Jet-4"; 
    jetn[5]= "Jet-5"; 
    
    for(int ij=0; ij<2; ij++){
      hWJetPt[ij][i] = new TH1F("hWJetPt_" + jetn[ij] + "_" + namech[i], "p_{T}^{Jet} " + jetn[ij] + " " + titlenamech[i],40,0,200);
      hWJetPt[ij][i]->GetXaxis()->SetTitle("p_{T}[GeV]");      
    }
    for(int ij=0; ij<6; ij++){
      hJetPt[ij][i] = new TH1F("hJetPt_" + jetn[ij] + "_" + namech[i], "p_{T}^{Jet} " + jetn[ij] + " " + titlenamech[i],40,0,200);
      hJetPt[ij][i]->GetXaxis()->SetTitle("p_{T}[GeV]");      
    }
    for(int ij=0; ij<4; ij++){
      hOJetPt[ij][i] = new TH1F("hOJetPt_" + jetn[ij] + "_" + namech[i], "p_{T}^{Jet} " + jetn[ij] + " " + titlenamech[i],40,0,200);
      hOJetPt[ij][i]->GetXaxis()->SetTitle("p_{T}[GeV]");      
    }

    
    TString dijetname[5][6];
    dijetname[0][1] = "Jet01";
    dijetname[0][2] = "Jet02";
    dijetname[0][3] = "Jet03";
    dijetname[0][4] = "Jet04";
    dijetname[0][5] = "Jet05";
    dijetname[1][2] = "Jet12";
    dijetname[1][3] = "Jet13";
    dijetname[1][4] = "Jet14";
    dijetname[1][5] = "Jet15";
    dijetname[2][3] = "Jet23";
    dijetname[2][4] = "Jet24";
    dijetname[2][5] = "Jet25";
    dijetname[3][4] = "Jet34";
    dijetname[3][5] = "Jet35";
    dijetname[4][5] = "Jet45";
    
    for(int ja=0; ja<5; ja++){
      for(int jb=ja+1; jb<6; jb++){
	hpTJet[ja][jb][i]   = new TH1F("hpTJet_"   + dijetname[ja][jb] + "_" + namech[i], "transverse pT of Dijets "   + dijetname[ja][jb] + " " + titlenamech[i],80,0,400);
	hMassJet[ja][jb][i] = new TH1F("hMassJet_" + dijetname[ja][jb] + "_" + namech[i], "transverse Mass of Dijets " + dijetname[ja][jb] + " " + titlenamech[i],80,0,400);
	hDRJet[ja][jb][i]   = new TH1F("hDRJet_"   + dijetname[ja][jb] + "_" + namech[i], "#Delta R of Dijets "        + dijetname[ja][jb] + " " + titlenamech[i],50,0,5);
	hDPhiJet[ja][jb][i] = new TH1F("hDPhiJet_" + dijetname[ja][jb] + "_" + namech[i], "#Delta #phi of Dijets "     + dijetname[ja][jb] + " " + titlenamech[i],80,0,4);
      }
    }
    
    for(int ja=0; ja<3; ja++){
      for(int jb=ja+1; jb<4; jb++){
	hOpTJet[ja][jb][i]   = new TH1F("hOpTJet_"   + dijetname[ja][jb] + "_" + namech[i], "transverse pT of Dijets "   + dijetname[ja][jb] + " " + titlenamech[i],80,0,400);
	hOMassJet[ja][jb][i] = new TH1F("hOMassJet_" + dijetname[ja][jb] + "_" + namech[i], "transverse Mass of Dijets " + dijetname[ja][jb] + " " + titlenamech[i],80,0,400);
	hODRJet[ja][jb][i]   = new TH1F("hODRJet_"   + dijetname[ja][jb] + "_" + namech[i], "#Delta R of Dijets "        + dijetname[ja][jb] + " " + titlenamech[i],50,0,5);
	hODPhiJet[ja][jb][i] = new TH1F("hODPhiJet_" + dijetname[ja][jb] + "_" + namech[i], "#Delta #phi of Dijets "     + dijetname[ja][jb] + " " + titlenamech[i],80,0,4);
      }
    }
    
    hpTWjj[i]      = new TH1F("hpTWjj_"      + namech[i] ,"pT jets coming from W "                  + titlenamech[i], 80, 0, 400);
    hInvMassWjj[i] = new TH1F("hInvMassWjj_" + namech[i] ,"Inv. Mass of jets coming from W "        + titlenamech[i], 80, 0, 400);
    hDRWjj[i]      = new TH1F("hDRWjj_"      + namech[i] ,"#Delta R_{jj} of jets coming from W "    + titlenamech[i], 50, 0, 5);
    hDPhiWjj[i]    = new TH1F("hDPhiWjj_"    + namech[i] ,"#Delta #phi_{jj} of jets coming from W " + titlenamech[i], 80, 0, 4);

    hpTOjj[i]      = new TH1F("hpTOjj_"      + namech[i] ,"pT jets coming from W "                  + titlenamech[i], 80, 0, 400);
    hInvMassOjj[i] = new TH1F("hInvMassOjj_" + namech[i] ,"Inv. Mass of jets coming from W "        + titlenamech[i], 80, 0, 400);
    hDROjj[i]      = new TH1F("hDROjj_"      + namech[i] ,"#Delta R_{jj} of jets coming from W "    + titlenamech[i], 50, 0, 5);
    hDPhiOjj[i]    = new TH1F("hDPhiOjj_"    + namech[i] ,"#Delta #phi_{jj} of jets coming from W " + titlenamech[i], 80, 0, 4);
    
    hpTminjj[i]   = new TH1F("hpTminjj_"      + namech[i] ,"Minimum pT^{jj} "          + titlenamech[i], 80, 0, 400);
    hpTmaxjj[i]   = new TH1F("hpTmaxjj_"      + namech[i] ,"Maximum pT^{jj} "          + titlenamech[i], 80, 0, 400);
    hInvMassjj[i] = new TH1F("hInvMassjj_" + namech[i] ,"Compatible Inv. Mass "     + titlenamech[i], 80, 0, 400);
    hDRminjj[i]   = new TH1F("hDRminjj_"   + namech[i] ,"Minimum #Delta R_{jj} "    + titlenamech[i], 50, 0, 5);
    hDRmaxjj[i]   = new TH1F("hDRmaxjj_"   + namech[i] ,"Maximum #Delta R_{jj} "    + titlenamech[i], 50, 0, 5);
    hDPhiminjj[i] = new TH1F("hDPhiminjj_" + namech[i] ,"Minimum #Delta #Phi_{jj} " + titlenamech[i], 80, 0, 4);
    hDPhimaxjj[i] = new TH1F("hDPhimaxjj_" + namech[i] ,"Maximum #Delta #Phi_{jj} " + titlenamech[i], 80, 0, 4);
    
    hOpTminjj[i]   = new TH1F("hOpTminjj_"   + namech[i] ,"Minimum pT^{jj} "          + titlenamech[i], 80, 0, 400);
    hOpTmaxjj[i]   = new TH1F("hOpTmaxjj_"   + namech[i] ,"Maximum pT^{jj} "          + titlenamech[i], 80, 0, 400);
    hOInvMassjj[i] = new TH1F("hOInvMassjj_" + namech[i] ,"Compatible Inv. Mass "     + titlenamech[i], 80, 0, 400);
    hODRminjj[i]   = new TH1F("hODRminjj_"   + namech[i] ,"Minimum #Delta R_{jj} "    + titlenamech[i], 50, 0, 5);
    hODRmaxjj[i]   = new TH1F("hODRmaxjj_"   + namech[i] ,"Maximum #Delta R_{jj} "    + titlenamech[i], 50, 0, 5);
    hODPhiminjj[i] = new TH1F("hODPhiminjj_" + namech[i] ,"Minimum #Delta #Phi_{jj} " + titlenamech[i], 80, 0, 4);
    hODPhimaxjj[i] = new TH1F("hODPhimaxjj_" + namech[i] ,"Maximum #Delta #Phi_{jj} " + titlenamech[i], 80, 0, 4);

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
  float fAccEvent_full_ttbb[2] = {0.0,0.0};
  float fAccEvent_full_ttjj[2] = {0.0,0.0};  
  float fAccEvent_vis_ttbb [2] = {0.0,0.0};
  float fAccEvent_vis_ttjj [2] = {0.0,0.0};

  int AccEvent_full_ttbb[2] = {0,0};
  int AccEvent_full_ttjj[2] = {0,0};
  int AccEvent_vis_ttbb [2] = {0,0};
  int AccEvent_vis_ttjj [2] = {0,0};

  // Uncertainties file name
  if(_syst) fname += "_SYS_" + syst_varname;

  /********************************
             Event Loop
  ********************************/
  std::cout << "--- Processing: " << theTree.GetEntries() << " events" << std::endl;
  
  for (Long64_t ievt=0; ievt<theTree.GetEntries();ievt++) {
    
    theTree.GetEntry(ievt);  
    print_progress(theTree.GetEntries(), ievt);
    
    // Jets 
    int NJets;
    NJets     = 0;    
    std::vector<TLorentzVector> vjets, vOjets, vWjets, vTjets;
    std::vector<int> vIndex, vWIndex, vOIndex, vTIndex;

    for(int ijet=0; ijet < Jet_px->size(); ijet++){
      
      TLorentzVector gjet;
      gjet.SetPxPyPzE((*Jet_px)[ijet], (*Jet_py)[ijet], (*Jet_pz)[ijet], (*Jet_E)[ijet]);

      if(gjet.Pt()>25 && std::abs(gjet.Eta())<2.5){

	vjets.push_back(gjet); // All Jets
	vIndex.push_back(ijet);

	if ((*Jet_Mom)[ijet] == 24){
	  vWjets.push_back(gjet); // Jets coming from W	
	  vWIndex.push_back(ijet);	  
	}	
	else if ((*Jet_Mom)[ijet] == 6){
	  vTjets.push_back(gjet); // Jets coming from W	
	  vTIndex.push_back(ijet);	  
	}	
	else{
	  vOjets.push_back(gjet); // Other Jets
	  vOIndex.push_back(ijet);
	}
      } // if(jet.pT > 25)
    }// for(jets)
    
        
    /***************************
       Categorization GenTop
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

    /*******************************************
          Categorization to use W Jets 
    *******************************************/

  TString CatEvt = "";
  // if (_ttbar_cat){
  if (cone_NaddJets == 0 && ttbar_category("tt", GenCat_ID) && vWjets.size() == 2){
    if      (ttbar_Wjjcategory(GenCat_ID).Contains("bb_") && vTjets.size() == 2) CatEvt = "bbWjj"; // Two b-jets from top inside acceptance
    else if (ttbar_Wjjcategory(GenCat_ID).Contains("b_" ) && vTjets.size() == 1) CatEvt = "bWbjj"; // One b-jet  from top inside acceptance
    else if ((ttbar_Wjjcategory(GenCat_ID) == "_LF" || ttbar_Wjjcategory(GenCat_ID) == "_c") && vTjets.size() == 0) CatEvt = "Wjj";   // No b-jets  from top inside acceptance
	
  }
  // }
  /*******************************************
             Dijet Invariant Mass 
  *******************************************/
    
  if (vWjets.size() == 2 && 
      (CatEvt == "bbWjj" || CatEvt == "bWjj" || CatEvt == "Wjj")){
  // if (vWjets.size() == 2 && vjets.size() == 6 && vTjets.size() == 2 && cone_NaddJets == 2){
      
      NJets = vjets.size();
      
      float pTWjj, MWjj, DRWjj, DPhiWjj;
      
      // Jets from W
      if(vWjets.size() == 2){
	TLorentzVector jet_WI  = vWjets[0];
	TLorentzVector jet_WII = vWjets[1];
	hWJetPt[0][Channel]->Fill(jet_WI.Pt());
	hWJetPt[1][Channel]->Fill(jet_WII.Pt());
	
	pTWjj   = (vWjets.at(0) + vWjets.at(1)).Pt();
	MWjj    = (vWjets.at(0) + vWjets.at(1)).M();
	DRWjj   = vWjets.at(0).DeltaR(vWjets.at(1));
	DPhiWjj = std::abs(vWjets.at(0).DeltaPhi(vWjets.at(1)));
	
	hpTWjj[Channel]     ->Fill(pTWjj); 
	hInvMassWjj[Channel]->Fill(MWjj); 
	hDRWjj[Channel]     ->Fill(DRWjj);   
	hDPhiWjj[Channel]   ->Fill(DPhiWjj);   

 	// MVA Bkg TREE
	b_pTjj   = pTWjj;
	b_Mjj    = MWjj;
	b_DRjj   = DRWjj;
	b_DPhijj = DPhiWjj;
	
	MVASignaltree->Fill();	

      }

    std::vector<TLorentzVector> seljets;
      
      for(int coljet = 0; coljet < 2; coljet++){
	
	if(coljet == 0) seljets = vjets; 
	else seljets = vOjets; 
	
	// Estimation of all the the paramenters of the dijet
	int pTminIndex[2], pTmaxIndex[2], InvMassIndex[2], DRminIndex[2], DRmaxIndex[2], DPhiminIndex[2], DPhimaxIndex[2];
	float minpTjj = 9999., maxpTjj = 0., Mjj = 0., minDeltaMjj = 9999.,  minDeltaRjj = 9999., maxDeltaRjj = 0., minDeltaPhijj = 9999., maxDeltaPhijj = 0.;
	
	// Loop over the first Jet
	for(int ijet=0; ijet < seljets.size(); ijet++){
	  
	  TLorentzVector jet_i = seljets[ijet];
	  
	  // Loop over the second Jet
	  for(int jjet=ijet+1; jjet < seljets.size(); jjet++){
	    
	    TLorentzVector jet_j = seljets[jjet];
	    
	    float DijetpT = (jet_i+jet_j).Pt(); 
	    
	    if(minpTjj > DijetpT){
	      minpTjj = DijetpT;
	      pTminIndex[0] = vIndex[ijet];
	      pTminIndex[1] = vIndex[jjet];
	    }
	    
	    if(maxpTjj < DijetpT){
	      maxpTjj = DijetpT;
	      pTmaxIndex[0] = vIndex[ijet]; pTmaxIndex[1] = vIndex[jjet];
	    }
	    
	    float DijetInvMass = (jet_i+jet_j).M(); 
	    float DeltaMjj = std::abs(DijetInvMass-80.3);
	    
	    if(minDeltaMjj > DeltaMjj){
	      minDeltaMjj = DeltaMjj;
	      Mjj = DijetInvMass;
	      InvMassIndex[0] = vIndex[ijet]; InvMassIndex[1] = vIndex[jjet];
	    }
	    
	    float DijetDeltaR  = jet_i.DeltaR(jet_j); 
	    
	    if(minDeltaRjj > DijetDeltaR){
	      minDeltaRjj = DijetDeltaR;
	      DRminIndex[0] = ijet; DRminIndex[1] = jjet;
	    }
	    
	    if(maxDeltaRjj < DijetDeltaR){
	      maxDeltaRjj = DijetDeltaR;
	      DRmaxIndex[0] = vIndex[ijet]; DRmaxIndex[1] = vIndex[jjet];
	    }
	    
	    float DijetDeltaPhi  = std::abs(jet_i.DeltaPhi(jet_j)); 
	    
	    if(minDeltaPhijj > DijetDeltaPhi){
	      minDeltaPhijj = DijetDeltaPhi;
	      DPhiminIndex[0] = vIndex[ijet]; DPhiminIndex[1] = vIndex[jjet];
	    }
	    
	    if(maxDeltaPhijj < DijetDeltaPhi){
	      maxDeltaPhijj = DijetDeltaPhi;
	      DPhimaxIndex[0] = vIndex[ijet]; DPhimaxIndex[1] = vIndex[jjet];
	    }
	    
	    // Dijet Variables
	    if(coljet == 0){
	      hpTJet   [ijet][jjet][Channel]->Fill(DijetpT);
	      hMassJet [ijet][jjet][Channel]->Fill(DijetInvMass);  
	      hDRJet   [ijet][jjet][Channel]->Fill(DijetDeltaR);
	      hDPhiJet [ijet][jjet][Channel]->Fill(DijetDeltaPhi);

	      // MVA TREE
	      b_pTjj   = DijetpT;
	      b_Mjj    = DijetInvMass;
	      b_DRjj   = DijetDeltaR;
	      b_DPhijj = DijetDeltaPhi;

	      MVAtree->Fill();
	    }
	    else{
	      hOpTJet   [ijet][jjet][Channel]->Fill(DijetpT);
	      hOMassJet [ijet][jjet][Channel]->Fill(DijetInvMass);  
	      hODRJet   [ijet][jjet][Channel]->Fill(DijetDeltaR);
	      hODPhiJet [ijet][jjet][Channel]->Fill(DijetDeltaPhi);	  


	      hpTOjj[Channel]     ->Fill(DijetpT); 
	      hInvMassOjj[Channel]->Fill(DijetInvMass); 
	      hDROjj[Channel]     ->Fill(DijetDeltaR);   
	      hDPhiOjj[Channel]   ->Fill(DijetDeltaPhi);   

	      // MVA Bkg TREE
	      b_pTjj   = DijetpT;
	      b_Mjj    = DijetInvMass;
	      b_DRjj   = DijetDeltaR;
	      b_DPhijj = DijetDeltaPhi;

	      MVABkgtree->Fill();
	    }
	  }// for(jjets)

	  if(coljet == 0) hJetPt[ijet][Channel]->Fill(jet_i.Pt());
	  else hOJetPt[ijet][Channel]->Fill(jet_i.Pt());
	}// for(ijet)
	
	
	if(coljet == 0){
	  hpTminjj[Channel]  ->Fill(minpTjj);   
	  hpTmaxjj[Channel]  ->Fill(maxpTjj);   
	  hInvMassjj[Channel]->Fill(Mjj); 
	  hDRminjj[Channel]  ->Fill(minDeltaRjj);   
	  hDRmaxjj[Channel]  ->Fill(maxDeltaRjj);   
	  hDPhiminjj[Channel]->Fill(minDeltaPhijj);   
	  hDPhimaxjj[Channel]->Fill(maxDeltaPhijj);   
	}
	else{
	  hOpTminjj[Channel]  ->Fill(minpTjj);   
	  hOpTmaxjj[Channel]  ->Fill(maxpTjj);   
	  hOInvMassjj[Channel]->Fill(Mjj); 
	  hODRminjj[Channel]  ->Fill(minDeltaRjj);   
	  hODRmaxjj[Channel]  ->Fill(maxDeltaRjj);   
	  hODPhiminjj[Channel]->Fill(minDeltaPhijj);   
	  hODPhimaxjj[Channel]->Fill(maxDeltaPhijj);   
	}
      }// for(coljet)
      hNJets[Channel]->Fill(NJets); 
    }// if(vWjets.size() == 0)
    
    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------

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
  TString outfname=dirname + "/hAcc-" + hname + "_" + fname  + ".root";
  //TString outfname=dirname + "/hAcc-" + hname + "_" + fname + ttbar_id + ".root";
  TFile *target  = new TFile(outfname,"RECREATE" );  
  
  for(int i=0; i<2; i++){    
    hNJets[i]->Write();

    // pT for EACH jet 
    for(int ij=0; ij<2; ij++) hWJetPt[ij][i]->Write();
    for(int ij=0; ij<6; ij++) hJetPt[ij][i] ->Write();
    for(int ij=0; ij<4; ij++) hOJetPt[ij][i]->Write();
    
    
    // pT for EACH dijet 
    for(int ja=0; ja<5; ja++){
      for(int jb=ja+1; jb<6; jb++) hpTJet[ja][jb][i]  ->Write();
    }      
    for(int ja=0; ja<5; ja++){
      for(int jb=ja+1; jb<6; jb++) hMassJet[ja][jb][i]->Write();
    }      
    for(int ja=0; ja<5; ja++){
      for(int jb=ja+1; jb<6; jb++) hDRJet[ja][jb][i]  ->Write();
    }      
    for(int ja=0; ja<5; ja++){
      for(int jb=ja+1; jb<6; jb++) hDPhiJet[ja][jb][i]->Write();
    }      
    
    for(int ja=0; ja<3; ja++){
      for(int jb=ja+1; jb<4; jb++) hOpTJet[ja][jb][i]  ->Write();
    }      
    for(int ja=0; ja<3; ja++){
      for(int jb=ja+1; jb<4; jb++) hOMassJet[ja][jb][i]->Write();
    }      
    for(int ja=0; ja<3; ja++){
      for(int jb=ja+1; jb<4; jb++) hODRJet[ja][jb][i]  ->Write();
    }      
    for(int ja=0; ja<3; ja++){
      for(int jb=ja+1; jb<4; jb++) hODPhiJet[ja][jb][i]->Write();
    }      
    
    // MAX and MIN parameter per event
    hpTWjj[i]     ->Write(); 
    hInvMassWjj[i]->Write(); 
    hDRWjj[i]     ->Write();   
    hDPhiWjj[i]   ->Write();   
    
    hpTminjj[i]  ->Write();
    hpTmaxjj[i]  ->Write();
    hInvMassjj[i]->Write();
    hDRminjj[i]  ->Write();
    hDRmaxjj[i]  ->Write();
    hDPhiminjj[i]->Write();
    hDPhimaxjj[i]->Write();
    
    hpTOjj[i]     ->Write(); 
    hInvMassOjj[i]->Write(); 
    hDROjj[i]     ->Write();   
    hDPhiOjj[i]   ->Write();   

    hOpTminjj[i]  ->Write();
    hOpTmaxjj[i]  ->Write();
    hOInvMassjj[i]->Write();
    hODRminjj[i]  ->Write();
    hODRmaxjj[i]  ->Write();
    hODPhiminjj[i]->Write();
    hODPhimaxjj[i]->Write();
    
    hJetMatch[i] ->Write();
    
    h2DJetMjjID_I [i] ->Write();
    h2DJetMjjID_II[i] ->Write();
  }//for(i)

  MVASignaltree->Write();
  MVAtree      ->Write();
  MVABkgtree   ->Write();
  
  target->Close();

  std::cout << "File saved as " << outfname << std::endl;
  
}


#endif

