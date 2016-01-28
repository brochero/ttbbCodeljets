void Histos(){

  // Arrange the SF in the TH2F
   TFile *MuSF   = TFile::Open("/home/brochero/ttbar/TopTrees_CATuples/ScaleFactors/SF_muon_IDISO_13TeV.root"); // IDISO mu
  //TFile *MuSF   = TFile::Open("/home/brochero/ttbar/TopTrees_CATuples/ScaleFactors/SF_electron_IDISO_13TeV.root"); // IDISO mu

  TH2F *hmuSF;
  hmuSF = (TH2F*) MuSF->Get("GlobalSF")->Clone("hmuSF"); // Trigger
  //hmuSF = (TH2F*) MuSF->Get("GlobalSF")->Clone("hmuSF"); // Trigger

  int nbin_YpT  = hmuSF->GetNbinsY();
  int nbin_Xeta = hmuSF->GetNbinsX();

  cout << "pT bins= "  << nbin_YpT  << endl; 
  cout << "eta bins= " << nbin_Xeta << endl; 

  std::vector<float> vmuXEta; // muon Eta values
  std::vector<float> vmuYPt; // muon Pt values

  cout << " ---- Histo SF_IDISO ----" << endl;

  for(int nby=1; nby <= nbin_YpT; nby++){
    cout <<  "bin " << nby << " width pT= "   << hmuSF->GetYaxis()->GetBinWidth(nby) 
	 << " starting at " << hmuSF->GetYaxis()->GetBinLowEdge(nby) << endl;
  }

  for(int nbx=1; nbx <= nbin_Xeta; nbx++){
    cout <<  "bin " << nbx << " width eta = " << hmuSF->GetXaxis()->GetBinWidth(nbx) 
	 << " starting at " << hmuSF->GetXaxis()->GetBinLowEdge(nbx) << endl;
  }


  //float hbinmuXeta[] = {-2.4 , -1.2 , -0.9 , 0.0 , 0.9 , 1.2 , 2.4}; // muon Eta values
  float hbinmuXeta[] = {-2.4 , 0.0 , 2.4}; // muon Eta values
  float hbinmuYpT [] = {20 , 30 , 40 , 50, 60, 1000}; // muon Pt values

  TH2F *GlobalSF;
  //GlobalSF    = new TH2F("TriggerSF","SF_{Trigger}^{#mu}(#eta,p_{t}) ",2,hbinmuXeta,5,hbinmuYpT);
  GlobalSF    = new TH2F("TriggerSF","SF_{Trigger}^{e}(#eta,p_{t}) ",2,hbinmuXeta,5,hbinmuYpT);
  //GlobalSF    = new TH2F("GlobalSF","SF_{ID,ISO}^{#mu}(#eta,p_{t}) ",2,hbinmuXeta,5,hbinmuYpT);
  //GlobalSF    = new TH2F("GlobalSF","SF_{ID,ISO}^{e}(#eta,p_{t}) ",6,hbinmuXeta,4,hbinmuYpT);

  float SF_T_mu[6] = {0.81,0.81,0.81,0.86,0.87,0.89};
  float SF_T_e [6] = {0.81,0.81,0.81,0.88,0.90,0.92};

  for(int nbx=1; nbx <= 2; nbx++){
    for(int nby=1; nby <= 5; nby++){
      GlobalSF->SetBinContent (nbx,nby, SF_T_e[nby]);
      GlobalSF->SetBinError   (nbx,nby, 0.1); 
    }
  }

  cout << " ---- New Histo ST_IDISO ----" << endl;
  // for(int nbx=1; nbx <= 6; nbx++){
  //   for(int nby=1; nby <= 4; nby++){
   
  //     // GlobalSF->SetBinContent (nbx,nby, 1.0);
  //     // GlobalSF->SetBinError   (nbx,nby, 0.1);
  //     if(nbx < 4){
  //     	GlobalSF->SetBinContent (nbx,nby, hmuSF->GetBinContent (4-nbx,nby));
  //     	GlobalSF->SetBinError   (nbx,nby, hmuSF->GetBinError   (4-nbx,nby));

  //     	cout << "(" << nbx << ","<< nby << ") = " << hmuSF->GetBinContent (4-nbx,nby) << endl;
  //     	cout << "(" << nbx << ","<< nby << ") = " << GlobalSF->GetBinContent (nbx,nby) << endl;
	
  //     	cout << "Error(" << nbx << ","<< nby << ") = " << hmuSF->GetBinError (4-nbx,nby) << endl;
  //     	cout << "Error(" << nbx << ","<< nby << ") = " << GlobalSF->GetBinError (nbx,nby) << endl;
  //     }
  //     else{
  //     	GlobalSF->SetBinContent (nbx,nby, hmuSF->GetBinContent (nbx-3,nby));
  //     	GlobalSF->SetBinError   (nbx,nby, hmuSF->GetBinError   (nbx-3,nby));
  //     	cout << "(" << nbx << ","<< nby << ") = " << hmuSF->GetBinContent (nbx-3,nby) << endl;
  //     	cout << "(" << nbx << ","<< nby << ") = " << GlobalSF->GetBinContent (nbx,nby) << endl;
	
  //     	cout << "Error(" << nbx << ","<< nby << ") = " << hmuSF->GetBinError (nbx-3,nby) << endl;
  //     	cout << "Error(" << nbx << ","<< nby << ") = " << GlobalSF->GetBinError (nbx,nby) << endl;
  //     }

  //   }
  // }



  //TFile *target  = new TFile("Output_SFTrigger.root","RECREATE");
  //TFile *target  = new TFile("Output_SFmu.root","RECREATE");
  TFile *target  = new TFile("Output_SFe.root","RECREATE");
  GlobalSF->Write();

  target->Close();
}
