#include <iostream> // Stream declarations
#include <vector>
#include <limits>

#include "TChain.h"
#include "TCut.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH1.h"
#include "TMath.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TROOT.h"

#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"

namespace TMVA {
  
  void Training(){
    std::string factoryOptions( "!V:!Silent:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    
    TString fname = "./tmva_example_multiple_background.root";
    
    TFile *input = NULL;
    input = TFile::Open(fname);

    TTree *signalTree = (TTree*)input->Get("TreeS");
    TTree *bkgTree    = (TTree*)input->Get("TreeB");

    /// global event weights per tree (see below for setting event-wise weights)
    Double_t signalWeight = 1.0;
    Double_t bkgWeight    = 1.0;

    // Create a new root output file.
    TString outfileName( "TMVASignalBackground0.root" );
    TFile*  outputFile = TFile::Open( outfileName, "RECREATE" );

    // ===================== background
    TMVA::Factory *factory = new TMVA::Factory( "TMVAMultiBkg0", outputFile, factoryOptions );
    // Example
    // factory->AddVariable( "varTree", "Variable Name", "units", 'F' );
    factory->AddVariable( "varTree", "Variable Name", "units", 'F' );

    factory->AddSignalTree    ( signalTree, signalWeight );
    factory->AddBackgroundTree( bkgTree,    bkgWeight    );

    //   factory->SetBackgroundWeightExpression("weight");
    TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

    // tell the factory to use all remaining events in the trees after training for testing:
    factory->PrepareTrainingAndTestTree( mycuts, mycutb,
					 "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

    // Boosted Decision Trees
    factory->BookMethod( TMVA::Types::kBDT, "BDTG",
			 "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=5" );
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    outputFile->Close();

    delete factory;
  }
}
