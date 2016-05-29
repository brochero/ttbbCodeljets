#include<vector>
#include "TLorentzVector.h"
#include<iostream>
#include "TTree.h"

void MVATrainingTree(TTree *TreeSignal, TTree *TreeBkg, std::vector<std::vector<float>> variables, std::vector<float> *Tvinput){
  
  int NDiJets = variables.at(4).size();
  int NVar = variables.size();

  for(int ip=0; ip<NDiJets; ip++){
    for(int iv=0; iv<Tvinput->size(); iv++){
      
      (*Tvinput)[iv] = variables.at(iv + 4).at(ip); // 4 -> ignore the first variables      
    }// for(iv)  

    // Is it a signal (add) or W pair?
    if      (variables.at(0).at(ip) == 1.0) TreeSignal->Fill(); // Add dijet
    else if (variables.at(0).at(ip) == 2.0) TreeBkg->Fill();    // Jets from W
    else if (variables.at(0).at(ip) == 3.0) TreeBkg->Fill();    // Jets from Top
  }// for(ip)  
  
}


std::vector<int> MVAResponse (std::vector<std::vector<float>> variables){
 
  std::vector<int> AddDijet;
  AddDijet.push_back(-1); // 0: Add Jet0  Index
  AddDijet.push_back(-1); // 1: Add Jet1  Index
  AddDijet.push_back(-1); // 2: Add DiJet Index

  int NDiJets = variables.at(4).size();

  float MaxRes = -999.;

  for(int ip=0; ip<NDiJets; ip++){        
    float MVARes = variables.at(3).at(ip);      
    if (MVARes > MaxRes){
      MaxRes = MVARes;
      AddDijet.at(0) = variables.at(1).at(ip);
      AddDijet.at(1) = variables.at(2).at(ip);
      AddDijet.at(2) = ip;
    }

  }// for(ip)  
  
  return AddDijet;
}

std::vector<int> WTopTagger (std::vector<std::vector<float>> variables){
 
  std::vector<int> AddDijet;
  AddDijet.push_back(-1); // 0: W Jet0  Index
  AddDijet.push_back(-1); // 1: W Jet1  Index
  AddDijet.push_back(-1); // 2: W DiJet Index
  AddDijet.push_back(-1); // 3: Top Jet0  Index
  AddDijet.push_back(-1); // 4: Top Jet1  Index
  AddDijet.push_back(-1); // 5: Top DiJet Index

  int NDiJets = variables.at(4).size();

  for(int ip=0; ip<NDiJets; ip++){        

    if(variables.at(0).at(ip) == 2.0){
      AddDijet.at(0) = variables.at(1).at(ip);
      AddDijet.at(1) = variables.at(2).at(ip);
      AddDijet.at(2) = ip;
    }      
    
    if(variables.at(0).at(ip) == 3.0){
      AddDijet.at(3) = variables.at(1).at(ip);
      AddDijet.at(4) = variables.at(2).at(ip);
      AddDijet.at(5) = ip;
    }      

  }// for(ip)  
  
  return AddDijet;
}
