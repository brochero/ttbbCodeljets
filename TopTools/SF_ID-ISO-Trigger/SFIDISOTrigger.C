#include<vector>
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include<iostream>


/***************************
  SF: ID, ISO and Trigger
***************************/
      
void SFIDISOTrigger(std::vector<float> &result,
		    TLorentzVector Lep, int channel,
		    TH2F *hmuIDISOSF, TH2F *hmuTriggerSF,
		    TH2F *heIDISOSF,  TH2F *heTriggerSF){

  /***************************
    SF: ID, ISO and Trigger
  ***************************/
  float SF_ID_ISO=1.0;
  float SF_ID_ISO_Error=0.0;
  float SF_Tr=1.0;
  float SF_Tr_Error=0.0;

  float SF_ID_ISO_Tr=1.0;


  /***************************
         Lepton: muon 
  ***************************/
  if(channel == 0){
    // Binning for Muon ID/ISO histo
    int Nbinmueta = hmuIDISOSF->GetNbinsX(); // X-axis (eta)
    int NbinmupT =  hmuIDISOSF->GetNbinsY();  // Y-axis (pT)
    
    std::vector<float> Vbinmueta; // muon eta values
    std::vector<float> VbinmupT;  // muon pT values
    
    // X-axis (eta)
    for(int nbx=1; nbx <= Nbinmueta; nbx++)  Vbinmueta.push_back( hmuIDISOSF->GetXaxis()->GetBinLowEdge(nbx) );
    Vbinmueta.push_back( hmuIDISOSF->GetXaxis()->GetBinLowEdge(Nbinmueta + 1) ); // Add upper edge of the last bin
    
    // Y-axis (pT)
    for(int nby=1; nby <= NbinmupT; nby++)  VbinmupT.push_back( hmuIDISOSF->GetYaxis()->GetBinLowEdge(nby) );
    //VbinmupT.push_back( hmuIDISOSF->GetYaxis()->GetBinLowEdge(NbinmupT + 1) ); // Add upper edge of the last bin
    VbinmupT.push_back(1000.0); // Add Max pT for the last bin
    
    // std::cout << "Number of X-axis bins (eta) = " << Nbinmueta << std::endl; 
    // for(int nbx=0; nbx < Vbinmueta.size(); nbx++) std::cout << "Range of X-axis bins (eta) = " << Vbinmueta[nbx] << std::endl; 
    // std::cout << "Number of Y-axis bins (pT)  = " << NbinmupT << std::endl; 
    // for(int nby=0; nby < VbinmupT.size(); nby++) std::cout << "Range of Y-axis bins (pT) = " << VbinmupT[nby] << std::endl; 

    // Binning for Muon Trigger histo
    int NbinTrmueta = hmuTriggerSF->GetNbinsX(); // X-axis (eta)
    int NbinTrmupT  = hmuTriggerSF->GetNbinsY(); // Y-axis (pT)
    
    std::vector<float> VbinTrmueta; // muon eta values
    std::vector<float> VbinTrmupT;  // muon pT values
    
    // X-axis (eta)
    for(int nbx=1; nbx <= NbinTrmueta; nbx++)  VbinTrmueta.push_back( hmuTriggerSF->GetXaxis()->GetBinLowEdge(nbx) );
    VbinTrmueta.push_back( hmuTriggerSF->GetXaxis()->GetBinLowEdge(NbinTrmueta + 1) ); // Add upper edge of the last bin
    
    // Y-axis (pT)
    for(int nby=1; nby <= NbinTrmupT; nby++)  VbinTrmupT.push_back( hmuTriggerSF->GetYaxis()->GetBinLowEdge(nby) );
    VbinTrmupT.push_back( hmuTriggerSF->GetYaxis()->GetBinLowEdge(NbinTrmupT + 1) ); // Add upper edge of the last bin
    
    //--------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------
    
    // -- ID/ISO
    for(int binmupT=0; binmupT<NbinmupT; binmupT++){
      if(Lep.Pt()>VbinmupT[binmupT] && Lep.Pt()<=VbinmupT[binmupT+1]){ 
	
	for(int binmueta=0; binmueta<Nbinmueta; binmueta++){
	  if(Lep.Eta()>Vbinmueta[binmueta] && Lep.Eta()<=Vbinmueta[binmueta+1]){ 
	    
	    SF_ID_ISO       = hmuIDISOSF->GetBinContent(binmueta+1,binmupT+1);
	    SF_ID_ISO_Error = hmuIDISOSF->GetBinError  (binmueta+1,binmupT+1);
	    
	    break;
	  }// if(Lep.Eta())
	}// for(binmueta)
	
      }// if(Lep.Pt())
    }// for(binmupT)
    
    // -- Trigger
    for(int binTrmupT=0; binTrmupT<NbinTrmupT; binTrmupT++){
      if(Lep.Pt()>VbinTrmupT[binTrmupT] && Lep.Pt()<=VbinTrmupT[binTrmupT+1]){ 
	
	for(int binTrmueta=0; binTrmueta<NbinTrmueta; binTrmueta++){
	  if(Lep.Eta()>VbinTrmueta[binTrmueta] && Lep.Eta()<=VbinTrmueta[binTrmueta+1]){ 
	    
	    SF_Tr       = hmuTriggerSF->GetBinContent(binTrmueta+1,binTrmupT+1);
	    SF_Tr_Error = hmuTriggerSF->GetBinError  (binTrmueta+1,binTrmupT+1);
	    
	    break;
	  }// if(Lep.Eta())
	}// for(binTrmueta)
	
      }// if(Lep.Pt())
    }// for(binTrmupT)
    
  }// if(channel == muon)
  
  /***************************
       Lepton: electron
  ***************************/
  if(channel == 1){
    // Binning for Electron ID/ISO histo
    int Nbineeta = heIDISOSF->GetNbinsX(); // X-axis (eta)
    int NbinepT  = heIDISOSF->GetNbinsY(); // Y-axis (pT)
    
    std::vector<float> Vbineeta; // electron eta values
    std::vector<float> VbinepT;  // electron pT values
    
    // X-axis (eta)
    for(int nbx=1; nbx <= Nbineeta; nbx++)  Vbineeta.push_back( heIDISOSF->GetXaxis()->GetBinLowEdge(nbx) );
    Vbineeta.push_back( heIDISOSF->GetXaxis()->GetBinLowEdge(Nbineeta + 1) ); // Add upper edge of the last bin
    
    // Y-axis (pT)
    for(int nby=1; nby <= NbinepT; nby++)  VbinepT.push_back( heIDISOSF->GetYaxis()->GetBinLowEdge(nby) );
    VbinepT.push_back( heIDISOSF->GetYaxis()->GetBinLowEdge(NbinepT + 1) ); // Add upper edge of the last bin
    
    // Binning for Electron Trigger histo
    int NbinTreeta = heTriggerSF->GetNbinsX(); // X-axis (eta)
    int NbinTrepT  = heTriggerSF->GetNbinsY(); // Y-axis (pT)
    
    std::vector<float> VbinTreeta; // electron eta values
    std::vector<float> VbinTrepT;  // electron pT values
    
    // X-axis (eta)
    for(int nbx=1; nbx <= NbinTreeta; nbx++)  VbinTreeta.push_back( heTriggerSF->GetXaxis()->GetBinLowEdge(nbx) );
    VbinTreeta.push_back( heTriggerSF->GetXaxis()->GetBinLowEdge(NbinTreeta + 1) ); // Add upper edge of the last bin
    
    // Y-axis (pT)
    for(int nby=1; nby <= NbinTrepT; nby++)  VbinTrepT.push_back( heTriggerSF->GetYaxis()->GetBinLowEdge(nby) );
    VbinTrepT.push_back( heTriggerSF->GetYaxis()->GetBinLowEdge(NbinTrepT + 1) ); // Add upper edge of the last bin
    
    //--------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------
    
    // -- ID/ISO  
    for(int binepT=0; binepT<NbinepT; binepT++){
      if(Lep.Pt()>VbinepT[binepT] && Lep.Pt()<=VbinepT[binepT+1]){ 
	
	for(int bineeta=0;bineeta<Nbineeta;bineeta++){
	  if(Lep.Eta()>Vbineeta[bineeta] && Lep.Eta()<=Vbineeta[bineeta+1]){ 
	    
	    SF_ID_ISO       = heIDISOSF->GetBinContent(bineeta+1,binepT+1);
	    SF_ID_ISO_Error = heIDISOSF->GetBinError(bineeta+1,binepT+1);
	    
	    break;
	  }// if(Lep.Eta())
	}// for(bineEta)
	
      }// if(Lep.Pt())
    }// for(binePt)
    
    // -- Trigger
    for(int binTrepT=0; binTrepT<NbinTrepT; binTrepT++){
      if(Lep.Pt()>VbinTrepT[binTrepT] && Lep.Pt()<=VbinTrepT[binTrepT+1]){ 
	
	for(int binTreeta=0;binTreeta<NbinTreeta;binTreeta++){
	  if(Lep.Eta()>VbinTreeta[binTreeta] && Lep.Eta()<=VbinTreeta[binTreeta+1]){ 
	    
	    SF_Tr       = heTriggerSF->GetBinContent(binTreeta+1,binTrepT+1);
	    SF_Tr_Error = heTriggerSF->GetBinError(binTreeta+1,binTrepT+1);
	    
	    break;
	  }// if(Lep.Eta())
	}// for(binTreeta)
	
      }// if(Lep.Pt())
    }// for(binTrePt)
    
  }// if(channel == electron)
  
  /*******************************************
   Trigger,ID & ISO Scale Factors/bin(Pt,Eta)
  *******************************************/
  // Control Information
  //std::cout << "Channel = " << channel << '\n' << "Lep.eta = " << Lep.Eta() << " -- Lep.pT = " << Lep.Pt() << '\n' << "SF = " << SF_ID_ISO << std::endl;
  
  result.push_back(SF_ID_ISO*SF_Tr); //[0]  
  result.push_back(SF_ID_ISO);       //[1]  
  result.push_back(SF_ID_ISO_Error); //[2]  
  result.push_back(SF_Tr);           //[3]  
  result.push_back(SF_Tr_Error);     //[4]  
  
}
 
