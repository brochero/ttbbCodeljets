#include "Plots.h"

void PlotsSyst(TString plots="2btag", bool LogScale=false, TString nSyst = ""){
  
  /****************
        Style
  ****************/
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1000);
  gStyle->SetOptStat("emruo");
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  
  gROOT->ProcessLine(".L tdrStyle.C");
  setTDRStyle();
  
  Int_t chatch = 1756;
  TColor *color = new TColor(chatch, 0.3, 0.5, 0.5, "", 0.45); // alpha = 0.5
  TString files  = dirnameIn + fl;  
  
  int col_ttbb = TColor::GetColor("#660000");
  int col_ttb  = TColor::GetColor("#ffcc00");
  int col_ttcc = TColor::GetColor("#cc6600");
  int col_ttc  = TColor::GetColor("#cc6600");
  int col_ttLF = TColor::GetColor("#ff0000");
  int col_tt   = TColor::GetColor("#FF7F7F");

  int col_ttbarBkg  = TColor::GetColor("#ff6565");
  int col_SingleTop = TColor::GetColor("#ff00ff");
  int col_WJets     = TColor::GetColor("#33cc33");
  int col_ZJets     = TColor::GetColor("#3366ff");
  int col_QCD       = TColor::GetColor("#ffff00");
  int col_VV        = TColor::GetColor("#ffffff");
  
  /****************
       Channel
  ****************/
  TString channel[3];
  channel[0] = "mujets";
  channel[1] = "ejets";  
  channel[2] = "lepjet"; 
 
  /****************
    ttbar Signal
  ****************/ 
  // ttbar categorization 
  std::vector<histos> ttbar_0;
  ttbar_0 = loadhistograms(plots, files + "_ttbar_PowhegPythia" + nSyst + "tt");
  setuphistograms(ttbar_0, col_tt);
  std::vector<histos> ttbar_0_ttbb;
  ttbar_0_ttbb = loadhistograms(plots, files + "_ttbar_PowhegPythia" + nSyst + "ttbb");
  setuphistograms(ttbar_0_ttbb, col_ttbb);
  std::vector<histos> ttbar_0_ttb;
  ttbar_0_ttb = loadhistograms(plots, files + "_ttbar_PowhegPythia" + nSyst + "ttb");
  setuphistograms(ttbar_0_ttb, col_ttb);
  std::vector<histos> ttbar_0_ttcc;
  ttbar_0_ttcc = loadhistograms(plots, files + "_ttbar_PowhegPythia" + nSyst + "ttcc");
  setuphistograms(ttbar_0_ttcc, col_ttcc);
  std::vector<histos> ttbar_0_ttLF;
  ttbar_0_ttLF = loadhistograms(plots, files + "_ttbar_PowhegPythia" + nSyst + "ttLF");
  setuphistograms(ttbar_0_ttLF, col_ttLF);

  /****************
     ttbar Bkg
  ****************/ 
  std::vector<histos> ttbar_bkg;
  ttbar_bkg = loadhistograms(plots, files + "_ttbar_PowhegPythiaBkg" + nSyst);
  setuphistograms(ttbar_bkg, col_ttbarBkg);
  std::vector<histos> ttOther;
  ttOther = addhistograms(ttbar_bkg, ttbar_0);
  setuphistograms(ttOther, col_ttbarBkg);
  
  /****************
     Draw Histos
  ****************/ 
  TCanvas *histocanvas;
  histocanvas = new TCanvas("plots", "Plots");

  for(unsigned int h=0; h<ttbar_0_ttbb.size(); h++){
    for(int ch=0; ch<3; ch++){
      
      histocanvas->cd();
      if(LogScale) histocanvas->cd()->SetLogy();

      ttbar_0_ttbb[h].hist[ch]->GetXaxis()->SetRange(ttbar_0_ttbb[h].hist[ch]->GetXaxis()->GetFirst(), 
						ttbar_0_ttbb[h].hist[ch]->GetXaxis()->GetLast());
      //ttbar_0_ttbb[h].hist[ch]->GetYaxis()->SetTitle("Events");
      ttbar_0_ttbb[h].hist[ch]->GetYaxis()->SetTitleOffset(1.2);
      ttbar_0_ttbb[h].hist[ch]->GetYaxis()->SetTitleSize(0.07);
      ttbar_0_ttbb[h].hist[ch]->GetYaxis()->SetLabelSize(0.055);
      ttbar_0_ttbb[h].hist[ch]->GetYaxis()->SetNdivisions(607);
      //ttbar_0_ttbb[h].hist[ch]->GetYaxis()->SetLabelSize(0.05);
      TGaxis *hYaxis = (TGaxis*)ttbar_0_ttbb[h].hist[ch]->GetYaxis();
      //hYaxis->SetMaxDigits(3);
      //ttbar_0_ttbb[h].hist[ch]->GetXaxis()->SetLabelSize(0.0);
      //ttbar_0_ttbb[h].hist[ch]->GetXaxis()->SetTitle("");


      // Produce enough vertical space for the legend 
      float MaxHisto;

      if(LogScale) MaxHisto = 8;
      else MaxHisto = 1.4; 
      float maxh;
      maxh = ttOther[h].hist[ch]->GetMaximum();
      if (maxh < ttbar_0_ttbb[h].hist[ch]->GetMaximum()) maxh = ttbar_0_ttbb[h].hist[ch]->GetMaximum();
      if (maxh < ttbar_0_ttb[h].hist[ch]->GetMaximum()) maxh = ttbar_0_ttb[h].hist[ch]->GetMaximum();
      if (maxh < ttbar_0_ttcc[h].hist[ch]->GetMaximum()) maxh = ttbar_0_ttcc[h].hist[ch]->GetMaximum();
      if (maxh < ttbar_0_ttLF[h].hist[ch]->GetMaximum()) maxh = ttbar_0_ttLF[h].hist[ch]->GetMaximum();

      ttbar_0_ttbb[h].hist[ch]->SetMaximum(MaxHisto*maxh);
      //if(LogScale) ttbar_0_ttbb[h].hist[ch]->SetMinimum(0.1);

      ttbar_0_ttbb[h].hist[ch]->Draw("hist");
      ttbar_0_ttb[h].hist[ch]->Draw("histSAME");
      ttbar_0_ttcc[h].hist[ch]->Draw("histSAME");
      ttbar_0_ttLF[h].hist[ch]->Draw("histSAME");
      ttOther[h].hist[ch]->SetLineStyle(4);
      ttOther[h].hist[ch]->Draw("histSAME");
      

      //ttbar_1[h].hist[ch]->Draw("histoSAME");

      /***********************
             Legends
      ***********************/
      TLegend *leg;
      float legx1=0.70;
      float legy1=0.80;
      float legx2=0.93;
      float legy2=0.90;
      leg = new TLegend(legx1,legy1,legx2,legy2);
      leg->SetFillColor(0);
      leg->SetLineColor(1);
      leg->SetTextFont(62);
      leg->SetTextSize(0.03);
      leg->SetNColumns(2);

      leg->AddEntry(ttOther[h].hist[ch],"t#bar{t}+other","l");
      leg->AddEntry(ttbar_0_ttLF[h].hist[ch],"t#bar{t}+LF","l");
      leg->AddEntry(ttbar_0_ttcc[h].hist[ch],"t#bar{t}+cc","l");
      leg->AddEntry(ttbar_0_ttb[h].hist[ch],"t#bar{t}+b","l");
      leg->AddEntry(ttbar_0_ttbb[h].hist[ch],"t#bar{t}+bb","l");
      //leg->AddEntry((TObject*)0,"","");
      leg->Draw("SAME");

      //-------------------------------------------------------
      // CMS Legend
      //-------------------------------------------------------
      TString htitleCMSChannel[3];
      htitleCMSChannel[0] = "#mu^{#pm}+jets channel";
      htitleCMSChannel[1] = "e^{#pm}+jets channel";
      htitleCMSChannel[2] = "l^{#pm}+jets channel";
      
      titlePr  = new TLatex(-20.,50.,"Preliminary");
      titlePr->SetNDC();
      titlePr->SetTextAlign(12);
      titlePr->SetX(0.25);
      titlePr->SetY(0.97);
      titlePr->SetTextColor(2);
      titlePr->SetTextFont(42);
      titlePr->SetTextSize(0.05);
      titlePr->SetTextSizePixels(24);
      titlePr->Draw("SAME");
      
      title  = new TLatex(-20.,50.,"CMS #sqrt{s} = 13TeV, L = 2.26 fb^{-1}");
      title->SetNDC();
      title->SetTextAlign(12);
      title->SetX(0.20);
      title->SetY(0.90);
      title->SetTextFont(42);
      title->SetTextSize(0.03);
      title->SetTextSizePixels(24);
      title->Draw("SAME");
      
      chtitle  = new TLatex(-20.,50.,htitleCMSChannel[ch]);
      chtitle->SetNDC();
      chtitle->SetTextAlign(12);
      chtitle->SetX(0.20);
      chtitle->SetY(0.86);
      chtitle->SetTextFont(42);
      chtitle->SetTextSize(0.03);
      chtitle->SetTextSizePixels(24);
      chtitle->Draw("SAME");

            
      /***********************
            Save Histos
      ***********************/    
      TString dirfigname_log;
      if(LogScale) dirfigname_log = "_log";
      else dirfigname_log = "";
      TString dirfigname_pdf;
      TString dirfigname_png;
      dirfigname_pdf = dirnameIn + "figuresSystComp_" + fl + "/ttbb" + nSyst + "/pdf" + dirfigname_log + "/";
      //dirfigname_png = dirnameIn + "figuresSystComp_" + fl + "/ttbb/png" + dirfigname_log + "/";
      // make a dir if it does not exist!!
      gSystem->mkdir(dirfigname_pdf,       kTRUE);
      histocanvas->SaveAs(dirfigname_pdf + ttbar_0[h].hist[ch]->GetName() + ".pdf");
      //gSystem->mkdir(dirfigname_png,       kTRUE);
      //histocanvas->SaveAs(dirfigname_png + WJets[h].hist[ch]->GetName() + ".png");
      
      // clear Canvas
      histocanvas->Clear();    

    }
  }
  
} //end Plots.C

void setuphistograms(std::vector<histos> histos, int color){
  
  for(unsigned int h = 0; h < histos.size(); h++){
    for(unsigned int ch=0; ch<3; ch++){  
      
      // Overflows in the last bin
      histos[h].hist[ch]->SetBinContent(histos[h].hist[ch]->GetNbinsX(),
					(histos[h].hist[ch]->GetBinContent(histos[h].hist[ch]->GetNbinsX()+1) + 
					 histos[h].hist[ch]->GetBinContent(histos[h].hist[ch]->GetNbinsX())
					 )
					);
      // !!!!!!!!!!!!!!!!!!!!!
      // Temporal Rebinning
      // !!!!!!!!!!!!!!!!!!!!!
      //histos[h].hist[ch]->Rebin(2);
      // Set overflows to 0
      histos[h].hist[ch]->SetBinContent(histos[h].hist[ch]->GetNbinsX()+1, 0);
      float Norm = histos[h].hist[ch]->Integral();
      if (Norm! = 0.0) Norm = 1.0/Norm;
      else Norm = 0.0;
      histos[h].hist[ch]->Scale(Norm);

      // Style
      // if(color == kRed+1)  histos[h].hist[ch]->SetLineColor(1); // Except ttbar
      // else 
      histos[h].hist[ch]->SetLineColor(color);
      //histos[h].hist[ch]->SetLineColor(1); 
      histos[h].hist[ch]->SetLineWidth(2); 
      //histos[h].hist[ch]->SetFillColor(color);
      //histos[h].hist[ch]->SetFillStyle(1001);

    }    
  }    
}

std::vector<histos> addhistograms(std::vector<histos> histos_0, std::vector<histos> histos_1){

  std::vector<histos> sum;
  histos histos;

  for(unsigned int h = 0; h < histos_0.size(); h++){
    for(unsigned int ch=0; ch<3; ch++){  
      histos.hist[ch] = (TH1F*)histos_0[h].hist[ch]->Clone();
      histos.hist[ch]->Add(histos_0[h].hist[ch], histos_1[h].hist[ch]);
    }
    sum.push_back(histos);
  }

  return sum;
}

std::vector<histos> addstack(std::vector<histos> stack_0, std::vector<histos> histos_0){

  for(unsigned int h = 0; h < histos_0.size(); h++){
    for(unsigned int ch=0; ch<3; ch++){
      stack_0[h].mc[ch]->Add(histos_0[h].hist[ch]); // add histo to stack
      stack_0[h].hist[ch] = (TH1F*) stack_0[h].mc[ch]->GetStack()->Last()->Clone(); // create histo with the final stack
    }
  }
  
  return stack_0;
}

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
    // lep + jets
    histofile.hist[2] = (TH1F*)histofile.hist[0]->Clone();
    histofile.hist[2]->Add(histofile.hist[0], histofile.hist[1]);
    histofile.hist[2]->SetName(histoname[h] + "_" + channel[2] + "_" + plots);
    
    sample.push_back(histofile);
    
  }  

  cout << "All the histograms from " << namefile << ".root have been loaded successfully!!!"  << endl; 
  return sample;
}

