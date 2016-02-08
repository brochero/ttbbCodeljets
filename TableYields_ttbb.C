#include <algorithm> // max(a,b)

//#include "/home/brochero/ttbar/TopCodeljets/TopResults/Yields_CSVT-v0/Yields.h" //Yields
//#include "/home/brochero/ttbar/TopCodeljets/TopResults/Yields_FullLumi-v3/Yields.h" //Yields
//#include "/home/brochero/ttbar/TopCodeljets/TopResults/Yields_6Jets_2btagM_1btagT/Yields.h" //Yields
#include "/home/brochero/ttbar/ttbbCodeljets/TopResults/Yields_Lumi2260-v3/Yields.h" //Yields

void TableYields_ttbb(bool Isttbb = true){

  TableYields("lepton", Isttbb);
  TableYields("6Jets",  Isttbb);
  TableYields("2btag",  Isttbb);
  TableYields("4btag",  Isttbb);
  
}

void TableYields(TString cutname, bool ttbb = false, TString outfiledir= "Lumi2260-v3"){

  int cut;
  if      (cutname=="lepton")   cut = 0;
  else if (cutname=="6Jets")    cut = 1;
  else if (cutname=="2btag")    cut = 2;
  else if (cutname=="4btag")    cut = 3;

  else{
    cout << "Invalid cut name!!!" << endl;
    break; 
  }

  int index_c = 15;// Index for Central Values

  // [Channels]
  float TTbar_c    [2];
  float err_TTbar_c[2];

  float TTbar_tt_c    [2];
  float err_TTbar_tt_c[2];

  float TTbar_ttLF_c    [2];
  float err_TTbar_ttLF_c[2];

  float TTbar_ttcc_c    [2];
  float err_TTbar_ttcc_c[2];

  float TTbar_ttb_c    [2];
  float err_TTbar_ttb_c[2];

  float TTbar_ttbb_c    [2];
  float err_TTbar_ttbb_c[2];

  float TTbar_ttjj_c    [2];
  float err_TTbar_ttjj_c[2];

  float VV_c    [2];
  float err_VV_c[2];

  float ST_c    [2];
  float err_ST_c[2];

  double ZJets    [2];
  double err_ZJets[2];
  
  double WJets    [2];
  double err_WJets[2];

  double TTbarBkg_c    [2];
  double err_TTbarBkg_c[2];

  double QCD    [2];
  double err_QCD[2];

  float backg_c    [2]; 
  float err_backg_c[2]; 

  float xsec_c [4][2];     //[0=central;1=stat;2=syst;3=lumi][channel]
  float acceptance [4][2]; //[0=central;1=stat;2=syst][channel]

  double data[2];
  
  // Yields
  for(int ch=0; ch<2; ch++){

    // TTbar    
    TTbar_c    [ch]  = ttbar_PowhegPythia    [index_c][2][ch][cut];
    err_TTbar_c[ch]  = err_ttbar_PowhegPythia[index_c][2][ch][cut];
    // TTbar+others    
    TTbar_tt_c    [ch]  = ttbar_PowhegPythiatt    [index_c][2][ch][cut];
    err_TTbar_tt_c[ch]  = err_ttbar_PowhegPythiatt[index_c][2][ch][cut];
    // TTbar+c    
    TTbar_ttLF_c    [ch]  = ttbar_PowhegPythiattLF    [index_c][2][ch][cut];
    err_TTbar_ttLF_c[ch]  = err_ttbar_PowhegPythiattLF[index_c][2][ch][cut];
    // TTbar+cc    
    TTbar_ttcc_c    [ch]  = ttbar_PowhegPythiattcc    [index_c][2][ch][cut];
    err_TTbar_ttcc_c[ch]  = err_ttbar_PowhegPythiattcc[index_c][2][ch][cut];
    // TTbar+b    
    TTbar_ttb_c    [ch]  = ttbar_PowhegPythiattb    [index_c][2][ch][cut];
    err_TTbar_ttb_c[ch]  = err_ttbar_PowhegPythiattb[index_c][2][ch][cut];
    // TTbar+ttbb    
    TTbar_ttbb_c    [ch]  = ttbar_PowhegPythiattbb    [index_c][2][ch][cut];
    err_TTbar_ttbb_c[ch]  = err_ttbar_PowhegPythiattbb[index_c][2][ch][cut];
    // TTbar+ttjj    
    TTbar_ttjj_c    [ch]  = ttbar_PowhegPythiattjj    [index_c][2][ch][cut];
    err_TTbar_ttjj_c[ch]  = err_ttbar_PowhegPythiattjj[index_c][2][ch][cut];

    // VV
    VV_c    [ch]     = WW[index_c][2][ch][cut] + WZ[index_c][2][ch][cut] + ZZ[index_c][2][ch][cut];
    err_VV_c[ch]     = sqrt(err_WW[index_c][2][ch][cut]**2 + err_WZ[index_c][2][ch][cut]**2 + err_ZZ[index_c][2][ch][cut]**2);
    // ST
    ST_c     [ch]    = tW[index_c][2][ch][cut] + tbarW[index_c][2][ch][cut] + tbar_tchannel[index_c][2][ch][cut] + t_tchannel[index_c][2][ch][cut];
    err_ST_c [ch]    = sqrt(err_tW[index_c][2][ch][cut]**2 + err_tbarW[index_c][2][ch][cut]**2 + err_tbar_tchannel[index_c][2][ch][cut]**2 + err_t_tchannel[index_c][2][ch][cut]**2);
    // ZJets
    ZJets    [ch]    = ZJets_M10to50_MCatNLO[index_c][2][ch][cut] + ZJets_M50_MCatNLO[index_c][2][ch][cut];
    err_ZJets[ch]    = sqrt(err_ZJets_M10to50_MCatNLO[index_c][2][ch][cut]**2 + err_ZJets_M50_MCatNLO[index_c][2][ch][cut]**2); 
    // WJets
    WJets    [ch]    = WJets_MCatNLO[index_c][2][ch][cut];
    err_WJets[ch]    = err_WJets_MCatNLO[index_c][2][ch][cut];
    // TTbar Background
    TTbarBkg_c    [ch]    = ttbar_PowhegPythiaBkg[index_c][2][ch][cut];
    err_TTbarBkg_c[ch]    = err_ttbar_PowhegPythiaBkg[index_c][2][ch][cut];
  }
  
  // QCD Muon
  QCD    [0]    = (QCD_MuEnr_20to30[index_c][2][0][cut] + 
		   QCD_MuEnr_30to50[index_c][2][0][cut] + 
		   //QCD_MuEnr_50to80[index_c][2][0][cut] + 
		   QCD_MuEnr_80to120[index_c][2][0][cut] + 
		   QCD_MuEnr_120to170[index_c][2][0][cut] + 
		   QCD_MuEnr_170to300[index_c][2][0][cut] + 
		   QCD_MuEnr_300to470[index_c][2][0][cut] + 
		   QCD_MuEnr_470to600[index_c][2][0][cut] + 
		   QCD_MuEnr_800to1000[index_c][2][0][cut]); //+ 
		   //QCD_MuEnr_1000toInf[index_c][2][0][cut]);
  
  err_QCD[0]    = sqrt(err_QCD_MuEnr_20to30[index_c][2][0][cut]**2 + 
		       err_QCD_MuEnr_30to50[index_c][2][0][cut]**2 + 
		       //err_QCD_MuEnr_50to80[index_c][2][0][cut]**2 + 
		       err_QCD_MuEnr_80to120[index_c][2][0][cut]**2 + 
		       err_QCD_MuEnr_120to170[index_c][2][0][cut]**2 + 
		       err_QCD_MuEnr_170to300[index_c][2][0][cut]**2 + 
		       err_QCD_MuEnr_300to470[index_c][2][0][cut]**2 + 
		       err_QCD_MuEnr_470to600[index_c][2][0][cut]**2 + 
		       err_QCD_MuEnr_800to1000[index_c][2][0][cut]**2);// + 
		       //err_QCD_MuEnr_1000toInf[index_c][2][0][cut]**2);
  // QCD Electron
  QCD    [1]    = (QCD_EGEnr_15to20[index_c][2][1][cut] + 
		   QCD_EGEnr_20to30[index_c][2][1][cut] + 
		   QCD_EGEnr_30to50[index_c][2][1][cut] + 
		   QCD_EGEnr_50to80[index_c][2][1][cut] + 
		   QCD_EGEnr_80to120[index_c][2][1][cut] + 
		   //QCD_EGEnr_120to170[index_c][2][1][cut] + 
		   QCD_EGEnr_170to300[index_c][2][1][cut] + 
		   QCD_EGEnr_300toInf[index_c][2][1][cut]);
  
  err_QCD[1]    = sqrt(err_QCD_EGEnr_15to20[index_c][2][1][cut]**2 + 
		       err_QCD_EGEnr_20to30[index_c][2][1][cut]**2 + 
		       err_QCD_EGEnr_30to50[index_c][2][1][cut]**2 + 
		       //err_QCD_EGEnr_50to80[index_c][2][1][cut]**2 + 
		       err_QCD_EGEnr_80to120[index_c][2][1][cut]**2 + 
		       //err_QCD_EGEnr_120to170[index_c][2][1][cut]**2 + 
		       err_QCD_EGEnr_170to300[index_c][2][1][cut]**2 + 
		       err_QCD_EGEnr_300toInf[index_c][2][1][cut]**2);
  
  // Data
  data[0] = DataSingleMu  [index_c][2][0][cut];
  data[1] = DataSingleEG  [index_c][2][1][cut];
  
  
  cout << "Cross Section extraction at 13TeV...................." << endl;

  float luminosity   = 2100.0; // [pb-1]
  float SysLumi      = 0.022;   // 2.2%
  double xsec_theory = 831.76;  // [pb] (13TeV NNLO+NNLL)
  
  

  // Table I: Cross Section Variations
  TString dirTeXname;
  if (ttbb) dirTeXname = "TopResults/TeX_Tables_ttbb_" + outfiledir + "/";
  else      dirTeXname = "TopResults/TeX_Tables_" + outfiledir + "/";
  // make a dir if it does not exist!!
  gSystem->mkdir(dirTeXname, kTRUE);


  for(int ch=0; ch<2; ch++){               // Channel
   
    /*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
                   Cross Section Central Value     
    *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/

    // Total background: VV + TW + ZJets + Fakes  
    backg_c[ch]     = VV_c[ch] + ST_c[ch] + ZJets[ch] + WJets[ch] + TTbarBkg_c[ch] + QCD[ch]; // Nom
    err_backg_c[ch] = sqrt(err_VV_c[ch]**2 + err_ST_c[ch]**2 + err_ZJets[ch]**2 + err_WJets[ch]**2 + err_TTbarBkg_c[ch]**2 + err_QCD[ch]**2); // Nom
    // Cross Section
    xsec_c[0][ch] = (xsec_theory/TTbar_c[ch]) * (data[ch] - backg_c[ch]);  // Nom
    // Statistical Uncertainty
    xsec_c[1][ch] = sqrt(data[ch]) * (xsec_c[0][ch]/(data[ch] - backg_c[ch]));  // stat
    // Systematic Uncertainty
    xsec_c[2][ch] = 0.0;
    // Luminosity Uncertainty
    xsec_c[3][ch] = 0.0;

    acceptance[0][ch] = TTbar_c[ch] / (luminosity * xsec_theory);
    acceptance[1][ch] = err_TTbar_c[ch] / (luminosity * xsec_theory);
    acceptance[2][ch] = (TTbar_c[ch]*0.0) / (luminosity * xsec_theory);
    
  }// for(ch) 
  

  TString XsecTeXSummary = dirTeXname + "Xsec_" + cutname + ".tex";
  FILE*   fxsecsummary = fopen(XsecTeXSummary, "w");

  fprintf(fxsecsummary,"\\begin\{tabular\}\{\|l\|c\|c\|\}\\hline\n");
  fprintf(fxsecsummary,"Channel \& \$\\sigma_{\\ttbar}\$[pb] \& Acceptance [\\%] \\\\\\hline\\hline\n");
  fprintf(fxsecsummary,"\$\{\\mu}\$+Jets \& \$ %.1f \\pm \%.1f \\pm %.1f \\pm %.1f \$ \& \$ %.3f \\pm %.3f \\pm %.3f \$ \\\\\\hline\n", xsec_c[0][0], xsec_c[1][0], xsec_c[2][0], xsec_c[3][0], acceptance[0][0]*100., acceptance[1][0]*100., acceptance[2][0]*100. );
  fprintf(fxsecsummary,"\$\{\\e}\$+Jets \& \$ %.1f \\pm \%.1f \\pm %.1f \\pm %.1f \$ \& \$ %.3f \\pm %.3f \\pm %.3f \$ \\\\\\hline\n",     xsec_c[0][1], xsec_c[1][1], xsec_c[2][1], xsec_c[3][1], acceptance[0][1]*100., acceptance[1][1]*100., acceptance[2][1]*100. );
  fprintf(fxsecsummary,"\\end\{tabular\}\n");

  fclose(fxsecsummary);

  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  cout << "---------------------------------------------------------" << endl;
  cout << cutname << endl;
  cout << "Xsec_mu+Jets = " << xsec_c[0][0] << " +/- " << xsec_c[1][0] << " +/- " << xsec_c[2][0] << " +/- " << xsec_c[3][0] <<endl;
  cout << "Xsec_e+Jets  = " << xsec_c[0][1] << " +/- " << xsec_c[1][1] << " +/- " << xsec_c[2][1] << " +/- " << xsec_c[3][1] <<endl;
  cout << "---------------------------------------------------------" << endl;
  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;


  // Yields

  // LaTeX table
  TString YieldsTeXfile  = dirTeXname + "Yields_" + cutname + ".tex";
  FILE*   Yields         = fopen(YieldsTeXfile, "w");

  fprintf(Yields,"\\begin\{tabular\}\{\|l\|c\|c\|\}\\hline\n");
  fprintf(Yields,"Source \& \$\{\\mu}+Jets\$ \& \$\{\\e}+Jets\$ \\\\\\hline\\hline\n");

  fprintf(Yields,"Single Top \&\$ %.1f \\pm %.1f \$ \&\$ %.1f \\pm %.1f \$   \\\\\\hline\n",    ST_c[0], err_ST_c[0], ST_c[1], err_ST_c[1]);
  fprintf(Yields,"Boson-Boson \&\$ %.1f \\pm %.1f \$  \&\$ %.1f \\pm %.1f \$   \\\\\\hline\n",    VV_c[0], err_VV_c[0], VV_c[1], err_VV_c[1]);

  fprintf(Yields,"Drell-Yan \&\$ %.1f \\pm %.1f \$  \&\$ %.1f \\pm %.1f \$   \\\\\\hline\n",    ZJets[0], err_ZJets[0], ZJets[1], err_ZJets[1] );
  fprintf(Yields,"W+Jets \&\$ %.1f \\pm %.1f \$  \&\$ %.1f \\pm %.1f \$    \\\\\\hline\n",    WJets[0],err_WJets[0], WJets[1], err_WJets[1] );
  fprintf(Yields,"\$\\ttbar\$ (Bkg) \&\$ %.1f \\pm %.1f \$  \&\$ %.1f \\pm %.1f \$    \\\\\\hline\n",    TTbarBkg_c[0],err_TTbarBkg_c[0], TTbarBkg_c[1], err_TTbarBkg_c[1] );
  fprintf(Yields,"QCD \&\$ %.1f \\pm %.1f \$  \&\$ %.1f \\pm %.1f \$    \\\\\\hline\\hline\n",    QCD[0],err_QCD[0], QCD[1], err_QCD[1] );

  fprintf(Yields,"Total Bkg \&\$ %.1f \\pm %.1f \$  \&\$ %.1f \\pm %.1f \$  \\\\\\hline\\hline\n",    backg_c[0], err_backg_c[0], backg_c[1],err_backg_c[1] );

  if(ttbb){
    fprintf(Yields,"\$\\ttbar\$+Others \&\$ %.1f \\pm %.1f \$  \&\$ %.1f \\pm %.1f \$  \\\\\\hline\n",   TTbar_tt_c[0], err_TTbar_tt_c[0], TTbar_tt_c[1],err_TTbar_tt_c[1] );
    fprintf(Yields,"\$\\ttbar\$+LF \&\$ %.1f \\pm %.1f \$  \&\$ %.1f \\pm %.1f \$  \\\\\\hline\n",   TTbar_ttLF_c[0], err_TTbar_ttLF_c[0], TTbar_ttLF_c[1],err_TTbar_ttLF_c[1] );
    fprintf(Yields,"\$\\ttbar\$+cc \&\$ %.1f \\pm %.1f \$  \&\$ %.1f \\pm %.1f \$  \\\\\\hline\n",   TTbar_ttcc_c[0], err_TTbar_ttcc_c[0], TTbar_ttcc_c[1],err_TTbar_ttcc_c[1] );
    fprintf(Yields,"\$\\ttbar\$+b \&\$ %.1f \\pm %.1f \$  \&\$ %.1f \\pm %.1f \$  \\\\\\hline\n",   TTbar_ttb_c[0], err_TTbar_ttb_c[0], TTbar_ttb_c[1],err_TTbar_ttb_c[1] );
    fprintf(Yields,"\$\\ttbar\$+bb \&\$ %.1f \\pm %.1f \$  \&\$ %.1f \\pm %.1f \$  \\\\\\hline\\hline\n",   TTbar_ttbb_c[0], err_TTbar_ttbb_c[0], TTbar_ttbb_c[1],err_TTbar_ttbb_c[1] );
    fprintf(Yields,"\$\\ttbar\$+jj \&\$ %.1f \\pm %.1f \$  \&\$ %.1f \\pm %.1f \$  \\\\\\hline\\hline\n",   TTbar_ttjj_c[0], err_TTbar_ttjj_c[0], TTbar_ttjj_c[1],err_TTbar_ttjj_c[1] );
  }

  else fprintf(Yields,"\$\\ttbar\$ (Signal) \&\$ %.1f \\pm %.1f \$  \&\$ %.1f \\pm %.1f \$    \\\\\\hline\\hline\n",   TTbar_c[0], err_TTbar_c[0], TTbar_c[1],err_TTbar_c[1] );

  fprintf(Yields,"Total MC \&\$ %.1f \\pm %.1f \$  \&\$ %.1f \\pm %.1f \$  \\\\\\hline\\hline\n",    backg_c[0] + TTbar_c[0], sqrt(err_backg_c[0]**2 + err_TTbar_c[0]**2), backg_c[1] + TTbar_c[1], sqrt(err_backg_c[1]**2 + err_TTbar_c[1]**2) );
  fprintf(Yields,"Data  \& %.1f \& %.1f \\\\\\hline\n",    data[0], data[1] );
  fprintf(Yields,"\\end\{tabular\}\n");


  fclose(Yields);
  
}// end Code
