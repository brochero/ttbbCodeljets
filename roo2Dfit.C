#include "roo2Dfit.h"

void roo2Dfit(TString SystVar = "", TString nModel = "RttbCon"){

  setTDRStyle();

  gSystem->Load("libRooFit") ;
  using namespace RooFit;

  HistoFit InFile[15];
  float ratio[15];
  int colors[15];

  InFile[data] = LoadSample("_DataSingleLep.root");

  // InFile[ttbb]   = LoadSample("_ttbar_PowhegPythiattbb"   + SystVar + ".root");
  // InFile[ttb]    = LoadSample("_ttbar_PowhegPythiattb"    + SystVar + ".root");
  // InFile[ttcc]   = LoadSample("_ttbar_PowhegPythiattcc"   + SystVar + ".root");
  // InFile[ttLF]   = LoadSample("_ttbar_PowhegPythiattLF"   + SystVar + ".root"); // Includes ttc
  // InFile[ttccLF] = LoadSample("_ttbar_PowhegPythiattccLF" + SystVar + ".root");

  InFile[ttbb]   = LoadSample("_ttbar_PowhegPythiattbb" + SystVar + ".root");
  InFile[ttb]    = LoadSample("_ttbar_PowhegPythiattb"  + SystVar + ".root");
  InFile[ttcc]   = LoadSample("_ttbar_PowhegPythiattcc" + SystVar + ".root");
  InFile[ttLF]   = LoadSample("_ttbar_PowhegPythiattLF" + SystVar + ".root"); // Includes ttc
  InFile[ttccLF] = LoadSample("_ttbar_PowhegPythiattccLF" + SystVar + ".root");

  InFile[WJets]     = LoadSample("_WJets_MCatNLO.root");
  InFile[ZJets]     = LoadSample("_ZJets_MCatNLO.root");
  //InFile[SingleTop] = LoadSample("_SingleTop" + SystVar + ".root");
  InFile[SingleTop] = LoadSample("_SingleTop" + SystVar + ".root");
  InFile[VV]        = LoadSample("_VV.root");
  InFile[QCD]       = LoadSample("_QCD.root");

  // Add Backgrounds
  InFile[Bkgtt]    = LoadSample("_ttbar_PowhegPythia" + SystVar + "Bkgtt.root"); // ttbarBkg + tt
  InFile[BkgOther] = LoadSample("_BkgOther" + SystVar + ".root"); 
  InFile[BkgFull]  = LoadSample("_BkgFull" + SystVar + ".root"); // BkgOther + Bkgtt


  TString name_ch[3];
  name_ch[0] = "muJets";
  name_ch[1] = "eJets";
  name_ch[2] = "LepJets";

  // Plots: Color
  colors[ttbb]     = TColor::GetColor("#660000");
  colors[ttb]      = TColor::GetColor("#ffcc00");
  colors[ttcc]     = TColor::GetColor("#cc6600");
  colors[ttLF]     = TColor::GetColor("#ff0000");
  colors[ttccLF]   = TColor::GetColor("#cc6600");
  colors[ttjj]     = TColor::GetColor("#3366ff");
  colors[Bkgtt]    = TColor::GetColor("#ff6565");
  colors[BkgOther] = TColor::GetColor("#8000ff");
  colors[BkgFull]  = TColor::GetColor("#33cc33");
  
  // Plots: Canvas
  TCanvas *canvas_comp = new TCanvas("canvas_comp",   "Component Plots");
  canvas_comp->Divide(2,2);
  TCanvas *canvas_data = new TCanvas("canvas_data",   "Data Plots");
  canvas_data->Divide(2,2);
  TCanvas *canvas_ratio_k = new TCanvas("canvas_ratio_k", "Parameters");
  canvas_ratio_k->Divide(2,1);
  
  
  TString dirfigname_pdf;
  dirfigname_pdf = dirnameIn + "FIT-" + nModel + "_figures_" + fl + SystVar + "/pdf/";
  // make a dir if it does not exist!!
  gSystem->mkdir(dirfigname_pdf,       kTRUE);
  TString Rfile = dirnameIn + "FIT-" + nModel + "_figures_" + fl + SystVar + "/FitResult.log";
  FILE* fResult = fopen(Rfile, "w");
  
  
  for (unsigned int ch=0; ch<3; ch++){
    
    float n_ttjj = InFile[ttbb].events[ch] + 
      InFile[ttcc].events[ch] + 
      InFile[ttb].events[ch] + 
      InFile[ttLF].events[ch];    
    
    ratio[ttbb]    = InFile[ttbb].events[ch]/n_ttjj; 
    ratio[ttcc]    = InFile[ttcc].events[ch]/n_ttjj; 
    ratio[ttb]     = InFile[ttb].events[ch]/n_ttjj; 
    ratio[ttccLF]  = InFile[ttccLF].events[ch]/n_ttjj; 
    
    // Variables
    RooRealVar *CSV2 = new RooRealVar("CSV2", "CSV for Jet 3", InFile[ttbb].hist1D[0][0]->GetXaxis()->GetXmin(), InFile[ttbb].hist1D[0][0]->GetXaxis()->GetXmax()); 
    RooRealVar *CSV3 = new RooRealVar("CSV3", "CSV for Jet 4", InFile[ttbb].hist1D[0][0]->GetXaxis()->GetXmin(), InFile[ttbb].hist1D[0][0]->GetXaxis()->GetXmax());  
    
    // Ratio at RECO level
    RooRealVar *RECO_ratio_ttbb   = new RooRealVar("RECO_ratio_ttbb",   "RECO ratio ttbb/ttjj",   ratio[ttbb],   ratio[ttbb],    ratio[ttbb]);
    RooRealVar *RECO_ratio_ttb    = new RooRealVar("RECO_ratio_ttb",    "RECO ratio ttb/ttjj",    ratio[ttb],    ratio[ttb],     ratio[ttb]);
    RooRealVar *RECO_ratio_ttcc   = new RooRealVar("RECO_ratio_ttcc",   "RECO ratio ttcc/ttjj",   ratio[ttcc],   ratio[ttcc],    ratio[ttcc]);
    RooRealVar *RECO_ratio_ttccLF = new RooRealVar("RECO_ratio_ttccLF", "RECO ratio ttccLF/ttjj", ratio[ttccLF], ratio[ttccLF] , ratio[ttccLF]);

    // Signal: Fitted ratios
    RooRealVar  *fit_ratio_ttbb   = new RooRealVar("fit_ratio_ttbb",   "FITTED ratio ttbb/ttjj",   ratio[ttbb],   0.0, 0.1); 
    RooRealVar  *fit_ratio_ttb    = new RooRealVar("fit_ratio_ttb",    "FITTED ratio ttb/ttjj",    ratio[ttb],    0.0, 1.0); 
    RooRealVar  *fit_ratio_ttcc   = new RooRealVar("fit_ratio_ttcc",   "FITTED ratio ttcc/ttjj",   ratio[ttcc],   0.0, 1.0); 
    RooRealVar  *fit_ratio_ttccLF = new RooRealVar("fit_ratio_ttccLF", "FITTED ratio ttccLF/ttjj", ratio[ttccLF], 0.0, 1.0); 
    RooRealVar  *fit_ratio_ttjj   = new RooRealVar("fit_ratio_ttjj",   "FITTED ratio ttjj/Total",  0.8,           0.0, 1.0); 

    // Constrained
    RooFormulaVar *fit_ratio_ttbb_con = new RooFormulaVar("fit_ratio_ttbb_con", "FITTED ratio ttbb/ttjj contrained", "fit_ratio_ttbb*(RECO_ratio_ttb/RECO_ratio_ttbb)", RooArgList(*fit_ratio_ttbb, *RECO_ratio_ttbb, *RECO_ratio_ttb));

    // Normalization Constant
    RooRealVar *k = new RooRealVar("k", "Normalization factor", 1.0, 0.50, 3.50);

    // Background
    RooRealVar  *fit_ratio_Bkgtt    = new RooRealVar("fit_ratio_Bkgtt",    "FITTED ratio bkgtt/FullBkg",    0.4, 0.0, 1.0); 
    RooRealVar  *fit_ratio_BkgOther = new RooRealVar("fit_ratio_BkgOther", "FITTED ratio bkgOther/FullBkg", 0.4, 0.0, 1.0); 

    // Signal
    RooRealVar *n_ttbb_var   = new RooRealVar("n_ttbb_var",   "number of ttbb events",   InFile[ttbb].events[ch],   InFile[ttbb].events[ch],   InFile[ttbb].events[ch]);
    RooRealVar *n_ttb_var    = new RooRealVar("n_ttb_var",    "number of ttb events",    InFile[ttb].events[ch],    InFile[ttb].events[ch],    InFile[ttb].events[ch]);
    RooRealVar *n_ttccLF_var = new RooRealVar("n_ttccLF_var", "number of ttccLF events", InFile[ttccLF].events[ch], InFile[ttccLF].events[ch], InFile[ttccLF].events[ch]);
    RooRealVar *n_ttjj_var   = new RooRealVar("n_ttjj_var",   "number of ttjj events",   n_ttjj, n_ttjj, n_ttjj);

    RooFormulaVar *kn_ttjj_var = new RooFormulaVar("kn_ttjj_var", "Normalized number of ttjj events", 
						  "k*n_ttjj_var", RooArgList(*k, *n_ttjj_var));
    
    // Background
    RooRealVar *n_Bkgtt_var    = new RooRealVar("n_Bkgtt_var",    "number of tt background events",    InFile[Bkgtt].events[ch],    InFile[Bkgtt].events[ch],    InFile[Bkgtt].events[ch]);
    RooRealVar *n_BkgOther_var = new RooRealVar("n_BkgOther_var", "number of Other background events", InFile[BkgOther].events[ch], InFile[BkgOther].events[ch], InFile[BkgOther].events[ch]);
    RooRealVar *n_BkgFull_var  = new RooRealVar("n_BkgFull_var",  "number of Total background events", InFile[BkgFull].events[ch],  InFile[BkgFull].events[ch],  InFile[BkgFull].events[ch]);

    RooFormulaVar *kn_BkgFull_var  = new RooFormulaVar("kn_BkgFull",  "Normalized number of background events", 
				   "k*n_BkgFull_var",  RooArgList(*k, *n_BkgFull_var));
    RooFormulaVar *kn_Bkgtt_var    = new RooFormulaVar("kn_Bkgtt",    "Normalized number of tt background events", 
				   "k*n_Bkgtt_var",    RooArgList(*k, *n_Bkgtt_var));
    RooFormulaVar *kn_BkgOther_var = new RooFormulaVar("kn_BkgOther", "Normalized number of Other background events", 
				   "k*n_BkgOther_var", RooArgList(*k, *n_BkgOther_var));
    

    RooArgList *arg_CSV = new RooArgList(*CSV2, *CSV3);
    //histograms    
    RooDataHist *data_his   = new RooDataHist("data_his",   "data Histogram",   
    			    *arg_CSV, InFile[data].hist2D[ch]);
    RooDataHist *ttbb_his   = new RooDataHist("ttbb_his",   "ttbb Histogram",   
			    *arg_CSV, InFile[ttbb].hist2D[ch]);
    RooDataHist *ttcc_his   = new RooDataHist("ttcc_his",   "ttbb Histogram",   
			    *arg_CSV, InFile[ttcc].hist2D[ch]);
    RooDataHist *ttb_his    = new RooDataHist("ttb_his",    "ttb Histogram",    
			    *arg_CSV, InFile[ttb].hist2D[ch]);
    RooDataHist *ttLF_his   = new RooDataHist("ttLF_his", "ttLF Histogram", 
			    *arg_CSV, InFile[ttLF].hist2D[ch]);
    RooDataHist *ttccLF_his = new RooDataHist("ttccLF_his", "ttccLF Histogram", 
			    *arg_CSV, InFile[ttccLF].hist2D[ch]);

    RooDataHist *Bkgtt_his    = new RooDataHist("Bkgtt_his",    "Bkgtt Histogram",     
			      *arg_CSV, InFile[Bkgtt].hist2D[ch]);
    RooDataHist *BkgOther_his = new RooDataHist("BkgOther_his", "BkgOther Histogram",    
			      *arg_CSV, InFile[BkgOther].hist2D[ch]);
    RooDataHist *BkgFull_his  = new RooDataHist("BkgFull_his",  "BkgFull Histogram",    
			      *arg_CSV, InFile[BkgFull].hist2D[ch]);

    RooDataHist *ttbbC_his   = new RooDataHist("ttbbC_his",   "ttbb Histogram",   
					       *arg_CSV, InFile[ttbb].hist2D[ch]);
       
    //pdf 
    RooHistPdf *ttbb_hispdf   = new RooHistPdf("ttbb_hispdf",   "PDF for ttbb",     RooArgSet(*arg_CSV), *ttbb_his);
    RooHistPdf *ttcc_hispdf   = new RooHistPdf("ttcc_hispdf",   "PDF for ttcc",     RooArgSet(*arg_CSV), *ttcc_his);
    RooHistPdf *ttb_hispdf    = new RooHistPdf("ttb_hispdf",    "PDF for ttb",      RooArgSet(*arg_CSV), *ttb_his);
    RooHistPdf *ttLF_hispdf   = new RooHistPdf("ttLF_hispdf",   "PDF for ttccLF",   RooArgSet(*arg_CSV), *ttLF_his);
    RooHistPdf *ttccLF_hispdf = new RooHistPdf("ttccLF_hispdf", "PDF for ttccLF",   RooArgSet(*arg_CSV), *ttccLF_his);

    RooHistPdf *Bkgtt_hispdf    = new RooHistPdf("Bkgtt_hispdf",    "PDF for ttbar Bkg", RooArgSet(*arg_CSV), *Bkgtt_his);
    RooHistPdf *BkgOther_hispdf = new RooHistPdf("BkgOther_hispdf", "PDF for Other Bkg", RooArgSet(*arg_CSV), *BkgOther_his);
    RooHistPdf *BkgFull_hispdf  = new RooHistPdf("BkgFull_hispdf",  "PDF for all Bkg",   RooArgSet(*arg_CSV), *BkgFull_his);

    /////////////////
    // Model
    /////////////////
    RooAddPdf *model_ttjj;
    
    if (nModel == "RttbCon")   model_ttjj = new RooAddPdf("model_ttjj", "Model For Signal (ttjj)", 
							  RooArgList(*ttbb_hispdf, *ttb_hispdf, *ttccLF_hispdf), 
							  RooArgList(*fit_ratio_ttbb, *fit_ratio_ttbb_con));
    if (nModel == "RttbFree")  model_ttjj = new RooAddPdf("model_ttjj", "Model For Signal (ttjj)", 
							  RooArgList(*ttbb_hispdf, *ttb_hispdf, *ttccLF_hispdf), 
							  RooArgList(*fit_ratio_ttbb, *fit_ratio_ttb));
    if (nModel == "RttccFree") model_ttjj = new RooAddPdf("model_ttjj", "Model For Signal (ttjj)", 
							  RooArgList(*ttbb_hispdf, *ttb_hispdf, *ttcc_hispdf, *ttLF_hispdf), 
							  RooArgList(*fit_ratio_ttbb, *fit_ratio_ttbb_con, *fit_ratio_ttcc));
    
    RooAddPdf *model = new RooAddPdf("model", "Model For Signal + Background", 
				     RooArgList(*model_ttjj, *Bkgtt_hispdf, *BkgOther_hispdf), 
				     RooArgList(*kn_ttjj_var, *n_Bkgtt_var, *n_BkgOther_var)); 
    
    
    // Plot I and II
    RooPlot *CSV2Tot_f = CSV2->frame();
    RooPlot *CSV3Tot_f = CSV3->frame();
    
    ttbbC_his->plotOn(CSV2Tot_f,LineColor(0),MarkerSize(0)); //Scale PDF to the HISTOS
    ttbbC_his->plotOn(CSV3Tot_f,LineColor(0),MarkerSize(0)); //Scale PDF to the HISTOS
    if (nModel == "RttccFree"){
      model->plotOn (CSV2Tot_f, Components(RooArgSet(*ttbb_hispdf, *ttb_hispdf, *ttcc_hispdf, *ttLF_hispdf)), LineColor(colors[ttjj]), RooFit::Name("CSV2Tot_f_ttjj"));
      model->plotOn (CSV3Tot_f, Components(RooArgSet(*ttbb_hispdf, *ttb_hispdf, *ttcc_hispdf, *ttLF_hispdf)), LineColor(colors[ttjj]));
    }
    else{
      model->plotOn (CSV2Tot_f, Components(RooArgSet(*ttbb_hispdf, *ttb_hispdf, *ttccLF_hispdf)), LineColor(colors[ttjj]), RooFit::Name("CSV2Tot_f_ttjj"));
      model->plotOn (CSV3Tot_f, Components(RooArgSet(*ttbb_hispdf, *ttb_hispdf, *ttccLF_hispdf)), LineColor(colors[ttjj]));
    }
    model->plotOn (CSV2Tot_f, Components(RooArgSet(*Bkgtt_hispdf, *BkgOther_hispdf)), LineColor(colors[BkgFull]), LineStyle(kDashed), RooFit::Name("CSV2Tot_f_BkgFull"));
    model->plotOn (CSV3Tot_f, Components(RooArgSet(*Bkgtt_hispdf, *BkgOther_hispdf)), LineColor(colors[BkgFull]), LineStyle(kDashed));
    
    // Plot III and IV
    RooPlot *CSV2_f = CSV2->frame();
    RooPlot *CSV3_f = CSV3->frame();

    ttbb_his->plotOn(CSV2_f,LineColor(0),MarkerSize(0)); //Scale PDF to the HISTOS
    ttbb_his->plotOn(CSV3_f,LineColor(0),MarkerSize(0)); //Scale PDF to the HISTOS

    model->plotOn(CSV2_f, Components(RooArgSet(*ttbb_hispdf)),   LineColor(colors[ttbb]),   RooFit::Name("CSV2_f_ttbb"));
    model->plotOn(CSV2_f, Components(RooArgSet(*ttb_hispdf)),    LineColor(colors[ttb]),    RooFit::Name("CSV2_f_ttb"));
    model->plotOn(CSV3_f, Components(RooArgSet(*ttbb_hispdf)),   LineColor(colors[ttbb]));
    model->plotOn(CSV3_f, Components(RooArgSet(*ttb_hispdf)),    LineColor(colors[ttb]));

    if (nModel == "RttccFree"){
      model->plotOn(CSV2_f, Components(RooArgSet(*ttcc_hispdf)), LineColor(colors[ttcc]), RooFit::Name("CSV2_f_ttcc"));
      model->plotOn(CSV2_f, Components(RooArgSet(*ttLF_hispdf)), LineColor(colors[ttLF]), RooFit::Name("CSV2_f_ttLF"));
      model->plotOn(CSV3_f, Components(RooArgSet(*ttcc_hispdf)), LineColor(colors[ttcc]));
      model->plotOn(CSV3_f, Components(RooArgSet(*ttLF_hispdf)), LineColor(colors[ttLF]));
    }
    else{
      model->plotOn(CSV2_f, Components(RooArgSet(*ttccLF_hispdf)), LineColor(colors[ttccLF]), RooFit::Name("CSV2_f_ttccLF"));
      model->plotOn(CSV3_f, Components(RooArgSet(*ttccLF_hispdf)), LineColor(colors[ttccLF]));
    }


    model->plotOn(CSV2_f, Components(RooArgSet(*Bkgtt_hispdf)),    LineColor(colors[Bkgtt]),    LineStyle(kDashed), RooFit::Name("CSV2_f_Bkgtt"));
    model->plotOn(CSV2_f, Components(RooArgSet(*BkgOther_hispdf)), LineColor(colors[BkgOther]), LineStyle(kDashed), RooFit::Name("CSV2_f_BkgOther"));
    model->plotOn(CSV3_f, Components(RooArgSet(*Bkgtt_hispdf)),    LineColor(colors[Bkgtt]),    LineStyle(kDashed));
    model->plotOn(CSV3_f, Components(RooArgSet(*BkgOther_hispdf)), LineColor(colors[BkgOther]), LineStyle(kDashed));


    model->fitTo(*data_his);
    
    // Plot I and II
    RooPlot *CSV2Data_f = CSV2->frame();
    RooPlot *CSV3Data_f = CSV3->frame();

    data_his->plotOn(CSV2Data_f, RooFit::Name("CSV2Data_f_Data"));
    data_his->plotOn(CSV3Data_f);

    model->plotOn(CSV2Data_f, LineColor(kRed), RooFit::Name("CSV2Data_f_Fit"));
    model->plotOn(CSV3Data_f, LineColor(kRed));

    model->plotOn(CSV2Data_f, Components(RooArgSet(*model_ttjj)), LineColor(colors[ttjj]), RooFit::Name("CSV2Data_f_Fitttjj"));
    model->plotOn(CSV3Data_f, Components(RooArgSet(*model_ttjj)), LineColor(colors[ttjj]));

    model->plotOn(CSV2Data_f, Components(RooArgSet(*Bkgtt_hispdf, *BkgOther_hispdf)), LineColor(colors[BkgFull]), LineStyle(kDashed), RooFit::Name("CSV2Data_f_FitBkg"));
    model->plotOn(CSV3Data_f, Components(RooArgSet(*Bkgtt_hispdf, *BkgOther_hispdf)), LineColor(colors[BkgFull]), LineStyle(kDashed));

    // Plot III and IV
    RooPlot *CSV2Datattbb_f = CSV2->frame();
    RooPlot *CSV3Datattbb_f = CSV3->frame();

    data_his->plotOn(CSV2Datattbb_f, RooFit::Name("CSV2Datattbb_f_Data"));
    data_his->plotOn(CSV3Datattbb_f);

    model->plotOn(CSV2Datattbb_f, LineColor(kRed), RooFit::Name("CSV2Datattbb_f_Fit"));
    model->plotOn(CSV3Datattbb_f, LineColor(kRed));

    model->plotOn(CSV2Datattbb_f, Components(RooArgSet(*ttbb_hispdf)), LineColor(colors[ttbb]), RooFit::Name("CSV2Datattbb_f_Fitttbb"));
    model->plotOn(CSV3Datattbb_f, Components(RooArgSet(*ttbb_hispdf)), LineColor(colors[ttbb]));
    
    if (nModel == "RttccFree"){
      model->plotOn(CSV2Datattbb_f, Components(RooArgSet(*ttcc_hispdf)), LineColor(colors[ttcc]), RooFit::Name("CSV2Datattbb_f_Fitttcc"));
      model->plotOn(CSV3Datattbb_f, Components(RooArgSet(*ttcc_hispdf)), LineColor(colors[ttcc]));
    }
    // Parameters
    RooArgSet *params = model->getVariables();
    params->Print("v");

    // double ttb       = fit_ratio_ttb->getVal();
    // cout << ttb << "------------------------"<< endl;

    double k_pf       = k->getVal();
    double k_pf_error = k->getError();

    double ratio_ttbb_pf       = fit_ratio_ttbb->getVal();
    double ratio_ttbb_pf_error = fit_ratio_ttbb->getError();

    // Plot: Ratio and Normalization Constant
    RooAbsReal *nll_ratio = model->createNLL(*data_his);

    RooPlot *ratio_ttbb_f = fit_ratio_ttbb->frame();
    nll_ratio->plotOn(ratio_ttbb_f, ShiftToZero());

    RooPlot *k_f = k->frame();
    nll_ratio->plotOn(k_f, ShiftToZero());

    // Efficiency Ratios
    float eff_Ratiobbjj[3]; // Eff_ttjj/Eff_ttbb    
    eff_Ratiobbjj[0] = 0.1885/0.4359; 
    eff_Ratiobbjj[1] = 0.1583/0.3710;
    eff_Ratiobbjj[2] = 0.1734/0.4034;
    // Acceptance Ratios
    float acc_Ratiobbjj[3]; // Acc_ttjj/Eff_ttbb    
    acc_Ratiobbjj[0] = 0.276/0.322; 
    acc_Ratiobbjj[1] = 0.275/0.320; 
    //acc_Ratiobbjj[2] = 0.275/0.320; 
    acc_Ratiobbjj[2] = 0.094428/0.109935; // Includes ttjj events in all channels

    float ratio_ttbb_Vis        = ratio_ttbb_pf  * eff_Ratiobbjj[ch];
    float ratio_ttbb_Vis_error  = ratio_ttbb_pf  * eff_Ratiobbjj[ch] * (ratio_ttbb_pf_error/ratio_ttbb_pf);
    float ratio_ttbb_Full       = ratio_ttbb_Vis * acc_Ratiobbjj[ch];
    float ratio_ttbb_Full_error = ratio_ttbb_Vis * acc_Ratiobbjj[ch] * (ratio_ttbb_Vis_error/ratio_ttbb_Vis);

    float N_ttjj_pf = k_pf * n_ttjj_var->getVal();
    float N_ttbb_pf = ratio_ttbb_pf * N_ttjj_pf;


    /***********************
            Plots
    ***********************/    
    TLegend *legCSV2Tot_f = new TLegend(0.2,0.7,0.48,0.9);
    legCSV2Tot_f->AddEntry(CSV2Tot_f->findObject("CSV2Tot_f_ttjj"),"t#bar{t}jj","l");
    legCSV2Tot_f->AddEntry(CSV2Tot_f->findObject("CSV2Tot_f_BkgFull"),"Bkg","l");
    
    canvas_comp->cd(1);
    CSV2Tot_f->Draw();
    CSV2Tot_f->SetMaximum(70);
    legCSV2Tot_f->Draw("SAME");

    canvas_comp->cd(2);
    CSV3Tot_f->Draw();

    TLegend *legCSV2_f = new TLegend(0.2,0.7,0.48,0.9);
    legCSV2_f->AddEntry(CSV2_f->findObject("CSV2_f_ttbb"),    "t#bar{t}b#bar{b}","l");
    legCSV2_f->AddEntry(CSV2_f->findObject("CSV2_f_ttb"),     "t#bar{t}bj","l");
    if(nModel == "RttccFree"){
      legCSV2_f->AddEntry(CSV2_f->findObject("CSV2_f_ttcc"),  "t#bar{t}c#bar{c}","l");
      legCSV2_f->AddEntry(CSV2_f->findObject("CSV2_f_ttLF"),  "t#bar{t}LF","l");
    }
    else legCSV2_f->AddEntry(CSV2_f->findObject("CSV2_f_ttccLF"),  "t#bar{t}c#bar{c} + t#bar{t}LF","l");
    legCSV2_f->AddEntry(CSV2_f->findObject("CSV2_f_Bkgtt"),   "t#bar{t} Bkg","l");
    legCSV2_f->AddEntry(CSV2_f->findObject("CSV2_f_BkgOther"),"Other Bkg","l");

    canvas_comp->cd(3);
    canvas_comp->cd(3)->SetLogy();
    CSV2_f->SetMinimum(1.e-3);
    CSV2_f->SetMaximum(10000);
    CSV2_f->Draw();
    legCSV2_f->Draw("SAME");

    canvas_comp->cd(4);
    canvas_comp->cd(4)->SetLogy();
    CSV3_f->SetMinimum(1.e-3);
    CSV3_f->SetMaximum(500);
    CSV3_f->Draw();


    TLegend *legCSV2Data_f = new TLegend(0.2,0.7,0.48,0.9);
    legCSV2Data_f->AddEntry(CSV2Data_f->findObject("CSV2Data_f_Data"),"Data","lp");
    legCSV2Data_f->AddEntry((TObject*)0,"FIT:","");
    legCSV2Data_f->AddEntry(CSV2Data_f->findObject("CSV2Data_f_Fit"),"Total","l");
    legCSV2Data_f->AddEntry(CSV2Data_f->findObject("CSV2Data_f_Fitttjj"),"t#bar{t}jj","l");
    legCSV2Data_f->AddEntry(CSV2Data_f->findObject("CSV2Data_f_FitBkg"),"Bkg","l");

    canvas_data->cd(1);
    float maxh = CSV2Data_f->GetMaximum();
    CSV2Data_f->SetMaximum(1.5*maxh);
    CSV2Data_f->Draw();
    legCSV2Data_f->Draw("SAME");

    canvas_data->cd(2);
    CSV3Data_f->Draw();
    
    TLegend *legCSV2Datattbb_f = new TLegend(0.2,0.7,0.48,0.9);
    legCSV2Datattbb_f->AddEntry(CSV2Datattbb_f->findObject("CSV2Datattbb_f_Data"),"Data","lp");
    legCSV2Datattbb_f->AddEntry((TObject*)0,"FIT:","");
    legCSV2Datattbb_f->AddEntry(CSV2Datattbb_f->findObject("CSV2Datattbb_f_Fit"),"Total","l");
    legCSV2Datattbb_f->AddEntry(CSV2Datattbb_f->findObject("CSV2Datattbb_f_Fitttbb"),"t#bar{t}b#bar{b}","l");
    if (nModel == "RttccFree") legCSV2Datattbb_f->AddEntry(CSV2Datattbb_f->findObject("CSV2Datattbb_f_Fitttcc"),"t#bar{t}c#bar{c}","l");

    canvas_data->cd(3);
    canvas_data->cd(3)->SetLogy();
    CSV2Datattbb_f->SetMaximum((2*CSV2Datattbb_f->GetMaximum()/10. * CSV2Datattbb_f->GetMaximum()));
    CSV2Datattbb_f->SetMinimum(0.01);
    CSV2Datattbb_f->Draw();
    legCSV2Datattbb_f->Draw("SAME");

    canvas_data->cd(4);
    canvas_data->cd(4)->SetLogy();
    CSV3Datattbb_f->SetMaximum((2*CSV3Datattbb_f->GetMaximum()/10. * CSV3Datattbb_f->GetMaximum()));
    CSV3Datattbb_f->SetMinimum(0.01);
    CSV3Datattbb_f->Draw();


    double DeltaSigma = 0.5; //one sigma variation. for two sigma, it should be 2.0 

    canvas_ratio_k->cd(1);
    ratio_ttbb_f->SetMaximum(0.6);
    ratio_ttbb_f->SetMinimum(0.0); 
    ratio_ttbb_f->GetYaxis()->SetNdivisions(607);
    ratio_ttbb_f->GetXaxis()->SetNdivisions(509); //(402)
    ratio_ttbb_f->Draw();

    TLine *line_ratio = new TLine(ratio_ttbb_f->GetXaxis()->GetXmin(), DeltaSigma, ratio_ttbb_f->GetXaxis()->GetXmax(), DeltaSigma);
    line_ratio->SetLineColor(kRed);
    line_ratio->Draw("SAME");
      
    canvas_ratio_k->cd(2);
    k_f->SetMaximum(0.6);
    k_f->SetMinimum(0.0);
    k_f->GetYaxis()->SetNdivisions(607);
    k_f->GetXaxis()->SetNdivisions(509); //(402)
    k_f->Draw();

    TLine *line_k = new TLine(k_f->GetXaxis()->GetXmin(), DeltaSigma, k_f->GetXaxis()->GetXmax(), DeltaSigma);
    line_k->SetLineColor(kRed);
    line_k->Draw("SAME");


    /***********************
         Save Histos
    ***********************/    
    canvas_comp->    SaveAs(dirfigname_pdf + "Comp_"   + name_ch[ch] + ".pdf");
    canvas_data->    SaveAs(dirfigname_pdf + "Data_"   + name_ch[ch] + ".pdf");
    canvas_ratio_k-> SaveAs(dirfigname_pdf + "Para_"   + name_ch[ch] + ".pdf");

    
    fprintf(fResult,"--------------------------------------------------------------\n" );
    fprintf(fResult,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" );
    fprintf(fResult,"--------------------------------------------------------------\n" );
    fprintf(fResult,"FINAL: RESULT %s \n", name_ch[ch].Data()      );
    fprintf(fResult,"FINAL: R_ttbb/ttjj (RECO) = %.4f +- %.4f \n", ratio_ttbb_pf, ratio_ttbb_pf_error );
    fprintf(fResult,"FINAL: R_ttbb/ttjj (Vis)  = %.4f +- %.4f \n", ratio_ttbb_Vis, ratio_ttbb_Vis_error  ); 
    fprintf(fResult,"FINAL: R_ttbb/ttjj (Full) = %.4f +- %.4f \n", ratio_ttbb_Full, ratio_ttbb_Full_error  );
    fprintf(fResult,"FINAL: k                  = %.4f +- %.4f \n", k_pf,  k_pf_error );

  } // for(ch)
  
  fclose(fResult);
  
}


HistoFit LoadSample(TString FileName){

  
  TFile* fInput  = new TFile(dirnameIn + fl +  FileName);
  
  TString Cut = "2btag";   
  
  HistoFit Output;
  
  // 3rd Jet
  Output.hist1D[0][0] = (TH1F*)fInput->Get("hCSV_Jet-2_mujets_" + Cut);
  Output.hist1D[0][1] = (TH1F*)fInput->Get("hCSV_Jet-2_ejets_"  + Cut);
  Output.hist1D[0][2] = (TH1F*)Output.hist1D[0][0]->Clone();
  Output.hist1D[0][2]->Add(Output.hist1D[0][0], Output.hist1D[0][1]);
  // 4th Jet
  Output.hist1D[1][0] = (TH1F*)fInput->Get("hCSV_Jet-3_mujets_" + Cut);
  Output.hist1D[1][1] = (TH1F*)fInput->Get("hCSV_Jet-3_ejets_"  + Cut);
  Output.hist1D[1][2] = (TH1F*)Output.hist1D[1][0]->Clone();
  Output.hist1D[1][2]->Add(Output.hist1D[1][0], Output.hist1D[1][1]);
  // 3rd and 4th Jet
  Output.hist2D[0] = (TH2F*)fInput->Get("h2DCSV_23Jet_mujets_" + Cut);
  Output.hist2D[1] = (TH2F*)fInput->Get("h2DCSV_23Jet_ejets_"  + Cut);
  Output.hist2D[2] = (TH2F*)Output.hist2D[0]->Clone();
  Output.hist2D[2]->Add(Output.hist2D[0], Output.hist2D[1]);
  // Numbewr of Events
  for(unsigned int i=0; i<3; i++) Output.events[i] = Output.hist1D[0][i]->Integral(0,1000);
  
  cout << "All the histograms from " << dirnameIn + fl +  FileName << " have been loaded successfully!!!" << endl;

  return Output;
}

