#include "Plots.h"

void IncludeSystPlotsBkg( TString SampleName = "", TString plots="2btag"){
  
  /****************
       Channel
  ****************/
  TString channel[3];
  channel[0] = "mujets";
  channel[1] = "ejets";  
  channel[2] = "lepjet"; 

  TString files  = dirnameIn + fl;
  
  /****************
    ttbar Signal
  ****************/ 
  // ttbar categorization 
  std::vector<histos> Sample;
  Sample = loadhistograms(plots, files + "_" + SampleName);

  // PileUp
  std::vector<histos> Sample_PileUp_Up;
  Sample_PileUp_Up = loadhistograms(plots, files + "_" + SampleName + "_SYS_PileUp_Up");
  std::vector<histos> Sample_PileUp_Down;
  Sample_PileUp_Down = loadhistograms(plots, files + "_" + SampleName + "_SYS_PileUp_Down");
  // JES
  std::vector<histos> Sample_JES_Up;
  Sample_JES_Up = loadhistograms(plots, files + "_" + SampleName + "_SYS_JES_Up");
  std::vector<histos> Sample_JES_Down;
  Sample_JES_Down = loadhistograms(plots, files + "_" + SampleName + "_SYS_JES_Down");
  // JER
  std::vector<histos> Sample_JER_Up;
  Sample_JER_Up = loadhistograms(plots, files + "_" + SampleName + "_SYS_JER_Up");
  std::vector<histos> Sample_JER_Nom;
  Sample_JER_Nom = loadhistograms(plots, files + "_" + SampleName + "_SYS_JER_Nom");
  std::vector<histos> Sample_JER_Down;
  Sample_JER_Down = loadhistograms(plots, files + "_" + SampleName + "_SYS_JER_Down");
  // b-tagging
  std::vector<histos> Sample_btag_Up;
  Sample_btag_Up = loadhistograms(plots, files + "_" + SampleName + "_SYS_btag_Up");
  std::vector<histos> Sample_btag_Down;
  Sample_btag_Down = loadhistograms(plots, files + "_" + SampleName + "_SYS_btag_Down");

  for(unsigned int h = 0; h < Sample.size(); h++){
    for(unsigned int ch=0; ch<2; ch++){// Only mu+Jets and e+Jets  
      for(int ibin = 1; ibin <= Sample[h].hist[ch]->GetNbinsX()+1; ibin++){      
      
	float central = Sample[h].hist[ch]->GetBinContent(ibin);
	// PileUp
	float PileUp_Up   = Sample_PileUp_Up[h].hist[ch]->GetBinContent(ibin);
	float PileUp_Down = Sample_PileUp_Down[h].hist[ch]->GetBinContent(ibin);
	// JES
	float JES_Up   = Sample_JES_Up[h].hist[ch]->GetBinContent(ibin);
	float JES_Down = Sample_JES_Down[h].hist[ch]->GetBinContent(ibin);
	// JER
	float JER_Up   = Sample_JER_Up[h].hist[ch]->GetBinContent(ibin);
	float JER_Nom  = Sample_JER_Nom[h].hist[ch]->GetBinContent(ibin);
	float JER_Down = Sample_JER_Down[h].hist[ch]->GetBinContent(ibin);
	// b-tagging
	float btag_Up   = Sample_btag_Up[h].hist[ch]->GetBinContent(ibin);
	float btag_Down = Sample_btag_Down[h].hist[ch]->GetBinContent(ibin);
	
	float PileUp=0.0;
	if(central != 0.0) PileUp = max(abs(central - PileUp_Up)/central,  abs(central - PileUp_Down)/central);
	float JES=0.0;
	if(central != 0.0) JES = max(abs(central - JES_Up)/central,  abs(central - JES_Down)/central);
	float JER=0.0;
	if(JER_Nom != 0.0)  JER    = max(abs(JER_Nom - JER_Up)/JER_Nom,  abs(JER_Nom - JER_Down)/JER_Nom);
	float btag=0.0;
	if(central != 0.0)  btag   = max(abs(central - btag_Up)/central, abs(central - btag_Down)/central);
	
	float Stat  = Sample[h].hist[ch]->GetBinError(ibin);
	float Syst  = sqrt(PileUp**2 + JES**2 + JER**2 + btag**2)*Sample[h].hist[ch]->GetBinContent(ibin);
	float Total = sqrt(Syst**2 + Stat**2); 
	
	if (Sample[h].hist[ch]->GetBinContent(ibin) == 0.0) Sample[h].hist[ch]->SetBinError(ibin, 0.0);
	else Sample[h].hist[ch]->SetBinError(ibin, Total);

	//cout << "Stat = " << Stat << " Syst = " << Syst << " Total Error = " << Total << endl;	
	
      }// for(bins)
    } // for(channel)    
  }  // for(histograms)  
  
  overwritehistograms(Sample, plots, files + "_" + SampleName);

} //end Plots.C

std::vector<histos> loadhistograms(TString plots, TString namefile){

  TFile *file=NULL; // new TFile(namefile);
  file = TFile::Open(namefile + ".root");
  
  if(!file){
    std::cerr << "ERROR: Could not open " <<  namefile  << " files!!!"  << std::endl;
    std::exit(0);
  }
  
  std::vector<histos> sample;
  histos histofile;
  
  TString channel[3];
  channel[0]="mujets";
  channel[1]="ejets";
  channel[2]="lepjets";
  
  std::vector<TString> histoname;
  
  // Histograms
  histoname.push_back("hPV");
  histoname.push_back("hMET");
  histoname.push_back("hLepPt");
  histoname.push_back("hLepEta");
  histoname.push_back("hLepPhi");
  histoname.push_back("hJetPt_Jet-0");
  histoname.push_back("hJetPt_Jet-1");
  histoname.push_back("hJetPt_Jet-2");
  histoname.push_back("hJetPt_Jet-3");
  histoname.push_back("hNJets");
  histoname.push_back("hNBtagJets");
  histoname.push_back("hCSV_Jet-0");
  histoname.push_back("hCSV_Jet-1");
  histoname.push_back("hCSV_Jet-2");
  histoname.push_back("hCSV_Jet-3");
  histoname.push_back("hmT");

  for(unsigned int h=0; h<histoname.size(); h++){
    for(unsigned int ch=0; ch<2; ch++) histofile.hist[ch] = (TH1F*)file->Get(histoname[h] + "_" + channel[ch] + "_" + plots);    

    sample.push_back(histofile);
    
  }  

  cout << "All the histograms from " << namefile << ".root have been loaded successfully!!!"  << endl; 
  return sample;
}

void overwritehistograms(std::vector<histos> newhistos, TString plots, TString namefile){

  TFile *file=NULL; // new TFile(namefile);
  file = TFile::Open(namefile + "SystError_" + plots + ".root", "RECREATE");
  
  if(!file){
    std::cerr << "ERROR: Could not open " <<  namefile  << " files!!!"  << std::endl;
    std::exit(0);
  }
  
  for(unsigned int h = 0; h < newhistos.size(); h++){
    for(unsigned int ch=0; ch<2; ch++){// Only mu+Jets and e+Jets  
      //newhistos[h].hist[ch]->Write("",TObject::kOverwrite);
      newhistos[h].hist[ch]->Write();
    }
  }
  
  cout << "All the histograms from " << namefile << ".root have been overwrote successfully!!!"  << endl; 
}

