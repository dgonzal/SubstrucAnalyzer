//bool print_ratio(const char *dirname, TString file_one, TString file_two, TString tag, TString hist);
#include <vector>
#include <iostream>
#include <string>

//ROOT includes
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TLine.h"


using namespace std;

bool plot_effi_and_ratio(vector<TGraphAsymmErrors*> histo, vector<TH1F*> ratio_histo, TString epsname, TString* type);
TH1F* plain_hist(const char *dirname, TString file, TString tag_one, TString hist1);
TH1F* scaled(const char *dirname, TString file, TString tag_one, TString h_entries, TString hist);
TH1F* ratio(const char *dirname, TString file, TString tag_one, TString tag_two, TString hist1, TString hist2);
TH1F* ratio_files(const char *dirname, TString file_one, TString file_two, TString tag_one, TString hist1);
TGraphAsymmErrors* efficiency(const char *dirname, TString file, TString tag_one, TString tag_two, TString hist1, TString hist2);
bool plot_nhists(vector<TH1F*> histo, TString epsname, TString* type);
bool plot_effi(vector<TGraphAsymmErrors*> histo, TString epsname, TString* type);
void set_style();

int main()
{
  //TString Folder[] = {"SVecQ_Selection_0Lepton/","SVecQ_Selection_1Lepton/","SVecQ_Selection_2Lepton/","SVecQ_Selection_3Lepton/"};

  const char *dirname="~/git_project/SubstrucAnalyzer/SubstrucAnalyzer/python/";

  //TString filename[] = {"fast_sim_flat","full_sim_flat"};
  TString filename[] = {"fastsim400","fullsim400"};

  TString legensname[] = {"Fast Sim","Full Sim"};
  
  TString tagnames[] = {"cmsTopTagPFJets","ak8PFJets"};//,"NSubjettinesPFJets"}; 
  TString hist[] = {"jet_pT","jet_phi","jet_eta","jet_mass"};
  TString hist_jetmatch[] = {"jet_matched_pT","jet_matched_phi","jet_matched_eta","jet_matched_mass"};
  TString hist_tagged[] = {"Topjet_pT","Topjet_phi","Topjet_eta","Topjet_mass"};
  TString hist_matched[] = {"Topjet_matched_pT","Topjet_matched_phi","Topjet_matched_eta","Topjet_matched_mass"};
  TString subjet_hists[] = {"subjet_pT","subjet_phi","subjet_eta","subjet_mass"};

  for(unsigned int i =0; i<sizeof(tagnames)/sizeof(TString); ++i){
    vector< TGraphAsymmErrors* > histo_effi;
    vector< TGraphAsymmErrors* > histo_effi_matched;

    //vector< TGraphAsymmErrors* > subjet_effi;

    vector< TH1F* > histo_ratio;
    vector< TH1F* > histo_ratio_matched;

    vector< TH1F* > histo_jet_N;
    vector< TH1F* > histo_jet_norm;
    //vector< TH1F* > histo_jet_ratio;

    vector< TH1F* > histo_subjet_norm;
    //vector< TH1F* > histo_subjet_ratio;

    vector< TH1F* > histo_jet_plain;
    vector< TH1F* > histo_topjet_plain;
    vector< TH1F* > histo_subjet_plain;

    for(unsigned int m =0; m<sizeof(hist)/sizeof(TString); ++m){
      for(unsigned int j =0; j<sizeof(filename)/sizeof(TString); ++j){

	histo_jet_plain.push_back(plain_hist(dirname,filename[j],tagnames[i],hist[m]));
	histo_topjet_plain.push_back(plain_hist(dirname,filename[j],tagnames[i],hist_tagged[m]));
	histo_subjet_plain.push_back(plain_hist(dirname,filename[j],tagnames[i],subjet_hists[m]));

	//histo_effi_matched.push_back(efficiency(dirname,filename[j],tagnames[i],tagnames[i],hist_matched[m], hist_jetmatch[m]));
	//histo_effi.push_back(efficiency(dirname,filename[j],tagnames[i],tagnames[i],hist_tagged[m],hist[m]));

	//histo_ratio_matched.push_back(ratio(dirname,filename[j],tagnames[i],tagnames[i],hist_matched[m], hist_jetmatch[m]));
	//histo_ratio.push_back(ratio(dirname,filename[j],tagnames[i],tagnames[i],hist_tagged[m],hist[m]));

	histo_jet_N.push_back(scaled(dirname,filename[j],tagnames[i],"jet_N",hist[m]));
	histo_jet_norm.push_back(scaled(dirname,filename[j],tagnames[i],hist[m],hist[m]));

	histo_subjet_norm.push_back(scaled(dirname,filename[j],tagnames[i],subjet_hists[m],subjet_hists[m]));
      }
      //histo_subjet_ratio.push_back(ratio_files(dirname,filename[0],filename[1],tagnames[i],subjet_hists[m]));
      //histo_jet_ratio.push_back(ratio_files(dirname,filename[0],filename[1],tagnames[i],hist[m]));
    }

    if(!plot_nhists(histo_jet_plain,"plots/jet_plain_"+tagnames[i]+".eps",legensname)) return 1;  

    if(tagnames[i]!="cmsTopTagPFJets") continue;

    if(!plot_nhists(histo_topjet_plain,"plots/topjet_plain_"+tagnames[i]+".eps",legensname)) return 1;  
    if(!plot_nhists(histo_subjet_plain,"plots/subjet_plain_"+tagnames[i]+".eps",legensname)) return 1;  

    if(!plot_effi(histo_effi_matched,"plots/effi_matched_"+tagnames[i]+".eps",legensname)) return 1;   
    if(!plot_effi(histo_effi,"plots/effi_"+tagnames[i]+".eps",legensname)) return 1;   
    if(!plot_effi_and_ratio(histo_effi, histo_ratio,"plots/effi_"+tagnames[i]+".eps", legensname))return 1;

    if(!plot_nhists(histo_ratio_matched,"plots/ratio_matched_"+tagnames[i]+".eps",legensname)) return 1;   
    if(!plot_nhists(histo_ratio,"plots/ratio_"+tagnames[i]+".eps",legensname)) return 1;  

    if(!plot_nhists(histo_jet_N,"plots/jet_N_"+tagnames[i]+".eps",legensname)) return 1;  
    if(!plot_nhists(histo_jet_norm,"plots/jet_norm_"+tagnames[i]+".eps",legensname)) return 1;  
    //if(!plot_nhists(histo_jet_ratio,"plots/jet_ratio_"+tagnames[i]+".eps",legensname)) return 1;  
    
    if(!plot_nhists(histo_subjet_norm,"plots/subjet_norm_"+tagnames[i]+".eps",legensname)) return 1;  
    //if(!plot_nhists(histo_subjet_ratio,"plots/subjet_ratio_"+tagnames[i]+".eps",legensname)) return 1;  

  }

  return 0;
}

TH1F* scaled(const char *dirname, TString file, TString tag_one, TString h_entries, TString hist)
{
  TFile* f1  = new TFile(dirname+file+".root","READ");

  TH1F* h1 = (TH1F*)f1->Get(tag_one+"/"+tag_one+"_"+hist);
  TH1F* entries = (TH1F*)f1->Get(tag_one+"/"+tag_one+"_"+h_entries);

  //cout<<entries->GetEntries()<<endl;
  h1->Scale(1/entries->GetEntries());

  return h1;
}


TGraphAsymmErrors* efficiency(const char *dirname, TString file, TString tag_one, TString tag_two, TString hist1, TString hist2)
{
  TFile* f1  = new TFile(dirname+file+".root","READ");

  TH1F* h1 = (TH1F*)f1->Get(tag_one+"/"+tag_one+"_"+hist1);
  TH1F* h2 = (TH1F*)f1->Get(tag_two+"/"+tag_two+"_"+hist2);
  

  h1->Rebin(5);
  h2->Rebin(5);

  //cout<<h1->GetName()<<"/"<<h2->GetName()<<dirname+file+".root" <<endl;

  TGraphAsymmErrors* gr = new TGraphAsymmErrors();
  gr->BayesDivide(h1,h2);
  gr->SetTitle(h1->GetTitle());

  return gr;
}

TH1F* plain_hist(const char *dirname, TString file, TString tag_one, TString hist1)
{
  TFile* f1  = new TFile(dirname+file+".root","READ");
  TH1F* h1 = (TH1F*)f1->Get(tag_one+"/"+tag_one+"_"+hist1);
  return h1;
}

TH1F* ratio_files(const char *dirname, TString file_one, TString file_two, TString tag_one, TString hist1)
{
  TFile* f1  = new TFile(dirname+file_one+".root","READ");
  TFile* f2  = new TFile(dirname+file_two+".root","READ");

  TH1F* h1 = (TH1F*)f1->Get(tag_one+"/"+tag_one+"_"+hist1);
  TH1F* h2 = (TH1F*)f2->Get(tag_one+"/"+tag_one+"_"+hist1);

  //h1->Rebin();
  //h2->Rebin();

  //cout<<h1->GetName()<<"/"<<h2->GetName()<<endl;
  h1->Divide(h2);
  h1->SetMarkerStyle(21);

  return h1;
}

TH1F* ratio(const char *dirname, TString file, TString tag_one, TString tag_two, TString hist1, TString hist2)
{

  //cout<<tag_one<<"/"<<tag_one<<"_"<<hist1<<"/"<<tag_two<<"/"<<tag_two<<"_"+hist2<<endl;
  
  TFile* f1  = new TFile(dirname+file+".root","READ");
  //TFile* f2  = new TFile(dirname+fiel_two,"READ");

  TH1F* h1 = (TH1F*)f1->Get(tag_one+"/"+tag_one+"_"+hist1);
  TH1F* h2 = (TH1F*)f1->Get(tag_two+"/"+tag_two+"_"+hist2);

  //cout<<h1->GetName()<<"/"<<h2->GetName()<<endl;
  h1->Divide(h2);

  return h1;
}

bool plot_effi(vector<TGraphAsymmErrors*> histo, TString epsname, TString* type)
{
  set_style();

  TCanvas* can = new TCanvas("can", "can", 600, 700); 
  can->cd();

  can->Print(epsname+"[");
  //can->SetLogy();
  
  for(unsigned int n = 0; n<histo.size()/2; ++n){ 
    TMultiGraph *mg = new TMultiGraph();
    int jumper = n *2;
    //int b = 0;
    //TString  histo_titel = histo[jumper]->GetTitle();
    
    TLegend * legend = new TLegend(0.7,0.8,.92,0.99);
    legend->SetTextFont(72);
    legend->SetTextSize(0.04);
    legend->SetFillColor(kWhite);

    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0.05);
    //pad1->SetLogy();
    pad1->Draw();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0.05);
    
    pad2->Draw();
    
    pad1->cd();

    //histo[jumper]->SetMaximum(0.5);
    histo[jumper]->SetLineColor(2);
    legend->AddEntry(histo[jumper],type[0],"l");
    mg->Add(histo[jumper]);
    //histo[jumper+1]->SetMaximum(0.5);
    histo[jumper+1]->SetLineColor(3);
    legend->AddEntry(histo[jumper+1],type[1],"l");
    mg->Add(histo[jumper+1]);

    /*
    for(unsigned int i = 0; i<histo.size(); ++i){ 
      if(histo_titel.EqualTo(histo[i]->GetTitle())){
	histo[i]->SetMaximum(histo[i]->GetMaximum()<1 ? 1. : 1.4*histo[i]->GetMaximum());
	histo[i]->SetLineColor(b+2);
	legend->AddEntry(histo[i],type[b],"l");
	mg->Add(histo[i]);
	++b;
      }
    }
    */
    //mg->SetLableSize(0.05);
    mg->Draw("ap");
    legend->Draw();
    /*
    pad2->cd();

    TGraphAsymmErrors *hnew = (TGraphAsymmErrors*)histo[jumper]->Clone("help_hist");
    //hnew->SetLabelSize(0.07,"Y");
    //hnew->SetLabelSize(0.00,"X");
    //hnew->Sumw2();
    hnew->SetMaximum(1.5);
    hnew->SetMinimum(.5);
    //hnew->SetStats(0);
    hnew->Divide(hnew,histo[jumper+1]);
    hnew->SetMarkerStyle(21);

    hnew->SetTitle("sdsd");

    //TLine* line = new TLine(histo[i]->GetXaxis()->GetXmin,0,histo[i]->GetXaxis()->GetXmax,0); 
    TLine* line = new TLine(hnew->GetXaxis()->GetXmin(),1,hnew->GetXaxis()->GetXmax(),1); 

    //cout<<hnew->GetXaxis()->GetXmax()<<endl;

    hnew->Draw("HISTp");
    line->Draw();
    */
    can->cd();


    can->Print(epsname);
  }
  can->Print(epsname+"]");
  can->Close();
  
  return true;
}


bool plot_nhists(vector<TH1F*> histo, TString epsname, TString* type)
{

  set_style();

  TCanvas* can = new TCanvas("can", "can", 600, 700); 
  can->cd();
  can->Print(epsname+"[");
  //can->SetLogy();
  



  for(unsigned int i = 0; i<histo.size(); i+=2){
    histo[i]->Rebin(2);
    histo[i+1]->Rebin(2);
    

    TLegend * legend = new TLegend(0.7,0.8,.92,0.99);
    legend->SetTextFont(72);
    legend->SetTextSize(0.04);
    legend->SetFillColor(kWhite);

    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0.05);
    pad1->SetLogy();
    pad1->Draw();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0.05);
   
    pad2->Draw();

    pad1->cd();
    /*
    histo[i]->GetXaxis()->SetLabelFont(63); //font in pixels
    histo[i]->GetXaxis()->SetLabelSize(16); //in pixels
    histo[i]->GetYaxis()->SetLabelFont(63); //font in pixels
    histo[i]->GetYaxis()->SetLabelSize(16); //in pixels
    */

    histo[i]->SetMaximum(1.4*histo[i]->GetMaximum());
    histo[i]->SetLineColor(kRed);
    histo[i+1]->SetLineColor(kGreen);
    histo[i]->Draw("E0");
    histo[i+1]->Draw("E0 same");

    legend->AddEntry(histo[i],type[0]);
    legend->AddEntry(histo[i+1],type[1]);
    legend->Draw();
    //can->cd();

    
    
    pad2->cd();

    //cout<<histo[i]->GetNbinsX()<<" "<<histo[i+1]->GetNbinsX()<<endl;

    TH1F *hnew = (TH1F*)histo[i]->Clone("help_hist");
    TH1F *h_denom = (TH1F*)histo[i+1]->Clone("help_denom");

    hnew->SetLabelSize(0.07,"Y");
    hnew->SetLabelSize(0.00,"X");
    
    //hnew->Sumw2();
    //h_denom->Sumw2();

    hnew->SetStats(0);
    hnew->Divide(hnew,h_denom,1,1,"B");
    hnew->SetMarkerStyle(21);

    hnew->SetMaximum(1.5);
    hnew->SetMinimum(.0);

    TGraphAsymmErrors* gr = new TGraphAsymmErrors();
    //histo[i]->Sumw2();
    //histo[i+1]->Sumw2();

    TH1F *h_nom = (TH1F*)histo[i]->Clone("help_nom");
    //TH1F *h_denom = (TH1F*)histo[i+1]->Clone("help_denom");
    //h_nom->Sumw2();
   
    gr->SetTitle(h_nom->GetTitle());

    //for(int ibin=0; ibin < h_nom->GetNbinsX(); ibin++)
    //  cout<<h_nom->GetBinContent(ibin)<<"/"<<h_denom->GetBinContent(ibin)<<endl;

    //gr->Divide(h_denom,h_nom);
    //gr->SetTitle(h1->GetTitle());
    //gr->SetLabelSize(0.07,"Y");
    //gr->SetLabelSize(0.00,"X");
    //gr->SetMaximum(1.5);
    //gr->SetMinimum(.5);
    
    TMultiGraph *mg = new TMultiGraph();
  
    mg->Add(gr);
    //hnew->SetTitle("sdsd");

    //TLine* line = new TLine(histo[i]->GetXaxis()->GetXmin,0,histo[i]->GetXaxis()->GetXmax,0); 
    TLine* line = new TLine(hnew->GetXaxis()->GetXmin(),1,hnew->GetXaxis()->GetXmax(),1); 

    //cout<<hnew->GetXaxis()->GetXmax()<<endl;
    //gStyle->SetErrorX(0);

    //mg->Draw("ap");
    hnew->Draw("p");
    line->Draw();
    can->cd();

    //gStyle->SetErrorX(1);
    
    can->Print(epsname);
    //can->Flush();
    //pad1->Close();
    //pad2->Close();

 }

  can->Print(epsname+"]");
  can->Close();
  
  return true;
}



void set_style()
{
  // general appearance and style
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0); 

  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
    
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);
    
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
    
  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  
  gStyle->UseCurrentStyle();

}
/*
TH1F *h1 = new TH1F("h1","test1",100,-3,3);
h1->FillRandom("gaus",200000);
h1->GetXaxis()->SetLabelFont(63); //font in pixels
h1->GetXaxis()->SetLabelSize(16); //in pixels
h1->GetYaxis()->SetLabelFont(63); //font in pixels
h1->GetYaxis()->SetLabelSize(16); //in pixels
TH1F *h2 = new TH1F("h2","test2",100,-3,3);
h2->FillRandom("gaus",100000);
TCanvas *c1 = new TCanvas("c1","example",600,700);
TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
pad1->SetBottomMargin(0);
pad1->Draw();
pad1->cd();
h1->DrawCopy();
h2->Draw("same");
c1->cd();
TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
pad2->SetTopMargin(0);
pad2->Draw();
pad2->cd();
h1->Sumw2();
h1->SetStats(0);
h1->Divide(h2);
h1->SetMarkerStyle(21);
h1->Draw("ep");
c1->cd();
*/
bool plot_effi_and_ratio(vector<TGraphAsymmErrors*> histo, vector<TH1F*> ratio_histo, TString epsname, TString* type)
{
  set_style();

  TCanvas* can = new TCanvas("can", "can", 600, 700); 
  can->cd();

  can->Print(epsname+"[");
  //can->SetLogy();
  
  for(unsigned int n = 0; n<histo.size()/2; ++n){ 
    TMultiGraph *mg = new TMultiGraph();
    int jumper = n *2;
    //int b = 0;
    //TString  histo_titel = histo[jumper]->GetTitle();
    
    TLegend * legend = new TLegend(0.7,0.8,.92,0.99);
    legend->SetTextFont(72);
    legend->SetTextSize(0.04);
    legend->SetFillColor(kWhite);

    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0.05);
    //pad1->SetLogy();
    pad1->Draw();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0.05);
    
    pad2->Draw();
    
    pad1->cd();

    //histo[jumper]->SetMaximum(0.5);
    histo[jumper]->SetLineColor(2);
    legend->AddEntry(histo[jumper],type[0],"l");
    mg->Add(histo[jumper]);
    //histo[jumper+1]->SetMaximum(0.5);
    histo[jumper+1]->SetLineColor(3);
    legend->AddEntry(histo[jumper+1],type[1],"l");
    mg->Add(histo[jumper+1]);

    //mg->SetLableSize(0.05);
    mg->Draw("ap");
    legend->Draw();
    
    pad2->cd();

    TH1F *hnew = (TH1F*)ratio_histo[jumper]->Clone("help_hist");
    hnew->Rebin(5);
    ratio_histo[jumper+1]->Rebin(5);
    hnew->SetLabelSize(0.07,"Y");
    hnew->SetLabelSize(0.00,"X");
    //hnew->Sumw2();
    hnew->Divide(hnew,ratio_histo[jumper+1]);
    hnew->SetMarkerStyle(21);

    hnew->SetMaximum(1.5);
    hnew->SetMinimum(.5);
    hnew->SetStats(0);

    //hnew->SetTitle("sdsd");

    //TLine* line = new TLine(histo[i]->GetXaxis()->GetXmin,0,histo[i]->GetXaxis()->GetXmax,0); 
    TLine* line = new TLine(hnew->GetXaxis()->GetXmin(),1,hnew->GetXaxis()->GetXmax(),1); 

    //cout<<hnew->GetXaxis()->GetXmax()<<endl;

    hnew->Draw("HISTp");
    line->Draw();
    
    can->cd();


    can->Print(epsname);
  }
  can->Print(epsname+"]");
  can->Close();
  
  return true;
}
