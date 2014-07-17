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


using namespace std;

TH1F* scaled(const char *dirname, TString file, TString tag_one, TString h_entries, TString hist);
TH1F* ratio(const char *dirname, TString file, TString tag_one, TString tag_two, TString hist1, TString hist2);
bool plot_nhists(vector<TH1F*> histo, TString epsname, TString* type);
void set_style();

int main()
{
  //TString Folder[] = {"SVecQ_Selection_0Lepton/","SVecQ_Selection_1Lepton/","SVecQ_Selection_2Lepton/","SVecQ_Selection_3Lepton/"};

  const char *dirname="/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_7_1_0_pre9/src/Substructure/SubstrucAnalyzer/python/";

  TString filename[] = {"ttbar_fs2","ttbar2"};
  
  TString tagnames[] = {"cmsTopTagPFJets"};//,"NSubjettinesPFJets"}; 
  TString hist[] = {"jet_pT","jet_phi","jet_eta","jet_mass"};
  TString hist_jetmatch[] = {"jet_matched_pT","jet_matched_phi","jet_matched_eta","jet_matched_mass"};
  TString hist_tagged[] = {"Topjet_pT","Topjet_phi","Topjet_eta","Topjet_mass"};
  TString hist_matched[] = {"Topjet_matched_pT","Topjet_matched_phi","Topjet_matched_eta","Topjet_matched_mass"};
  TString subjet_hists[] = {"subjet_pT","subjet_phi","subjet_eta","subjet_mass"};

  for(unsigned int i =0; i<sizeof(tagnames)/sizeof(TString); ++i){
    vector< TH1F* > histo_ratio;
    vector< TH1F* > histo_ratio_matched;
    vector< TH1F* > histo_jet_N;
    vector< TH1F* > histo_jet_norm;
    vector< TH1F* > histo_subjet_norm;
    for(unsigned int m =0; m<sizeof(hist)/sizeof(TString); ++m){
      for(unsigned int j =0; j<sizeof(filename)/sizeof(TString); ++j){
	histo_ratio_matched.push_back(ratio(dirname,filename[j],tagnames[i],tagnames[i],hist_matched[m], hist_jetmatch[m]));
	histo_ratio.push_back(ratio(dirname,filename[j],tagnames[i],tagnames[i],hist_tagged[m],hist[m]));
	histo_jet_N.push_back(scaled(dirname,filename[j],tagnames[i],"jet_N",hist[m]));
	histo_jet_norm.push_back(scaled(dirname,filename[j],tagnames[i],hist[m],hist[m]));
	histo_subjet_norm.push_back(scaled(dirname,filename[j],tagnames[i],subjet_hists[m],subjet_hists[m]));
      }
    }
    if(!plot_nhists(histo_ratio_matched,"plots/ratio_matched_"+tagnames[i]+".eps",filename)) return 1;   
    if(!plot_nhists(histo_ratio,"plots/ratio_"+tagnames[i]+".eps",filename)) return 1;   
    if(!plot_nhists(histo_jet_N,"plots/jet_N_"+tagnames[i]+".eps",filename)) return 1;  
    if(!plot_nhists(histo_jet_norm,"plots/jet_norm_"+tagnames[i]+".eps",filename)) return 1;  
    if(!plot_nhists(histo_subjet_norm,"plots/subjet_norm_"+tagnames[i]+".eps",filename)) return 1;  
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


bool plot_nhists(vector<TH1F*> histo, TString epsname, TString* type)
{

  set_style();

  TCanvas* can = new TCanvas("can", "can", 600, 600); 
  can->cd();

  can->Print(epsname+"[");
  can->SetLogy();
  
  TLegend * legend = new TLegend(0.7,0.8,.92,0.99);
  legend->SetTextFont(72);
  legend->SetTextSize(0.04);
  legend->SetFillColor(kWhite);
  
  TString  histo_titel = histo[0]->GetTitle();
  int b =0;

  for(unsigned int i = 0; i<histo.size(); ++i){ 
    if(histo_titel==histo[i]->GetTitle()){
      histo[i]->SetMaximum(histo[i]->GetMaximum()<1 ? 1.2 : 1.4*histo[i]->GetMaximum());
      histo[i]->Draw("same");
      histo[i]->SetLineColor(b+2);
      legend->AddEntry(histo[i],type[b]);
      b+=1;
    }
    else{
      
      legend->Draw();
      can->Print(epsname);
      legend->Clear();
      b=0;
      histo[i]->SetMaximum(histo[i]->GetMaximum()<1 ? 1.2 : 1.4*histo[i]->GetMaximum());
      histo[i]->Draw();
      histo[i]->SetLineColor(b+2);
      legend->AddEntry(histo[i],type[b]);
      b+=1;
      histo_titel=histo[i]->GetTitle();
    }
  }

  legend->Draw();
  can->Print(epsname);

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
