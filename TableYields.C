#include<string>
#include<iostream>
#include<sstream>
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<set>
#include<vector>
#include <sys/stat.h>

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
#include "TError.h"

struct Yields{
  std::vector<float>   Evt;
  std::vector<float>   Error;
  TString              Name;
};

Yields loadhistoYields(int SelCut, TString TName, TString namefile);
void   EntryPrinter   (FILE *file, Yields Sample);

void display_usage()
{
  std::cout << "\033[1;37musage:\033[1;m skimfile cutindex [options]" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "    -i inputfile  Input file without .root" << std::endl;
  std::cout << "    -o name in the output file. Default: Yields" << std::endl;
  std::cout << "    -d Input file directory. Default directory: TopResults" << std::endl;
  std::cout << "    -h displays this help message and exits " << std::endl;
  std::cout << "" << std::endl;
}



int main (int argc, char *argv[]){

  const char * _output   = "Yields";
  const char * _input    = 0;
  const char * _dir      = "TopResults/";
  const char * _cut      = 0;
  
  // Arguments used
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
      if( strcmp(argv[i],"-cut") == 0 ){
	_cut= argv[i+1];
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
  TString outfiledir(_output);
  TString fdir(_dir);
  TString cutname(_cut);

  ///////////////////////////////////////
  // Please, IGNORE. Temporal solution //
  ///////////////////////////////////////
  TCanvas *mydummycanvas=new TCanvas();// 
  ///////////////////////////////////////

  int cut;
  if      (cutname=="lepton")   cut = 0;
  else if (cutname=="6Jets")    cut = 1;
  else if (cutname=="2btag")    cut = 2;
  else if (cutname=="3btag")    cut = 3;

  else{
    std::cerr << "Invalid cut name!!!" << std::endl;
    std::exit(0); 
  }

  std::vector<Yields> Samples;

  Yields Data;
  Data = loadhistoYields(cut, "Data", fdir + fname + "_DataSingleLep");
  Samples.push_back(Data);

  Yields ttbar;
  ttbar = loadhistoYields(cut, "\\ttbar", fdir + fname + "_ttbar_PowhegPythia");
  Samples.push_back(ttbar);
  Yields ttbar_aMCatNLO;
  ttbar_aMCatNLO = loadhistoYields(cut, "ttbar aMCatNLO", fdir + fname + "_ttbar_MCatNLOPythia");
  Samples.push_back(ttbar_aMCatNLO);
  Yields ttbar_Madgraph;
  ttbar_Madgraph = loadhistoYields(cut, "ttbar Madgraph", fdir + fname + "_ttbar_MadgraphPythia");
  Samples.push_back(ttbar_Madgraph);
    
  Yields ttbb;
  ttbb = loadhistoYields(cut, "\\ttbar\\bbbar", fdir + fname + "_ttbar_PowhegPythiattbb");
  Samples.push_back(ttbb);
  Yields ttb;
  ttb = loadhistoYields(cut, "\\ttbar\\qb", fdir + fname + "_ttbar_PowhegPythiattb");
  Samples.push_back(ttb);
  Yields ttcc;
  ttcc = loadhistoYields(cut, "\\ttbar\\ccbar", fdir + fname + "_ttbar_PowhegPythiattcc");
  Samples.push_back(ttcc);
  Yields ttLF;
  ttLF = loadhistoYields(cut, "\\ttbar LF", fdir + fname + "_ttbar_PowhegPythiattLF");
  Samples.push_back(ttLF);
  Yields ttjj;
  ttjj = loadhistoYields(cut, "\\ttbar\\jj", fdir + fname + "_ttbar_PowhegPythiattjj");
  Samples.push_back(ttjj);

  Yields ttOther;
  ttOther = loadhistoYields(cut, "\\ttbar\\ Other", fdir + fname + "_BkgOther");
  Samples.push_back(ttOther);

  Yields tW;
  tW = loadhistoYields(cut, "Single Top", fdir + fname + "_SingleTop");
  Samples.push_back(tW);

  Yields WJets;
  WJets = loadhistoYields(cut, "W+Jets", fdir + fname + "_WJets_MCatNLO");
  Samples.push_back(WJets);

  Yields ZJets;
  ZJets = loadhistoYields(cut, "Z+Jets", fdir + fname + "_ZJets_MCatNLO");
  Samples.push_back(ZJets);

  Yields QCD;
  QCD = loadhistoYields(cut, "QCD", fdir + fname + "_QCD");
  Samples.push_back(QCD);


  // LaTeX table
  TString dirTeXname;
  dirTeXname = "TopResults/TeX_Tables_" + fname + "_" + outfiledir + "/";
  // make a dir if it does not exist!!
  gSystem->mkdir(dirTeXname, kTRUE);

  // Yields
  TString YieldsTeXfile  = dirTeXname + "Yields_" + cutname + ".tex";
  FILE*   Yields         = fopen(YieldsTeXfile, "w");

  fprintf(Yields,"\\begin{tabular}{|l|r@{~$\\pm$~}l|r@{~$\\pm$~}l|r@{~$\\pm$~}l|}\\hline\n");
  fprintf(Yields,"Source & \\multicolumn{2}{c|}{${\\mu}$+Jets} & \\multicolumn{2}{c|}{${\\e}$+Jets} & \\multicolumn{2}{c|}{${\\ell}$+Jets} \\\\\\hline\\hline\n");
  
  for(int ns = 0; ns < Samples.size(); ns++)  EntryPrinter(Yields, Samples.at(ns));

  fprintf(Yields,"\\end{tabular}\n");

  fclose(Yields);

  std::cout << "Table with yields has been done...................." << std::endl;
  std::cout << "File is " << YieldsTeXfile << std::endl;

  return 0;  

}// end Code


Yields loadhistoYields(int SelCut, TString TName,TString namefile){

  TFile *file=NULL; // new TFile(namefile);
  file = TFile::Open(namefile + ".root");

  if(!file){
    std::cerr << "ERROR: Could not open " <<  namefile  << " files!!!"  << std::endl;
    std::exit(0);
  }

  if(SelCut > 3 || SelCut < 0) {
    std::cerr << "ERROR: Incorrect selection cut name!!! "  << std::endl;
    std::exit(0);
  }

  float aYields[4][3];
  float aError [4][3];

  Yields output;

  // Histograms
  TH1F *YieldsSample;
  YieldsSample = (TH1F*)file->Get("Yields");

  int nbin = 1;
  for(int nc = 0; nc < 4; nc++){
    for(int nch = 0; nch < 3; nch++){

      aYields[nc][nch] = YieldsSample->GetBinContent(nbin);
      aError [nc][nch] = YieldsSample->GetBinError  (nbin);

      nbin++;      
    }
  }  

  output.Name = TName;
  std::cout << output.Name << std::endl;

  TString chname[3];
  chname[0] = "mu+Jets ";
  chname[1] = "e+Jets  ";
  chname[2] = "Lep+Jets";

  for(int nch = 0; nch < 3; nch++){
    output.Evt  .push_back(aYields[SelCut][nch]);
    output.Error.push_back(aError [SelCut][nch]);

    std::cout << chname[nch] << ": " << output.Evt[nch] << " +/- "  << output.Error[nch] << std::endl;
  }

  std::cout << "All Yields from " << namefile << ".root have been loaded successfully!!!"  << std::endl;
  return output;
}


void  EntryPrinter(FILE *file, Yields Sample){
  fprintf(file,"%s & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f   \\\\\\hline\n", 
   	  (Sample.Name).Data(), 
	  Sample.Evt[0], Sample.Error[0],
	  Sample.Evt[1], Sample.Error[1],
	  Sample.Evt[2], Sample.Error[2]);
}
