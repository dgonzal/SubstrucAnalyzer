#include "TSystemDirectory.h"
#include "TList.h"
#include "TSystemFile.h"
#include "TString.h"
#include "TF1.h"
#include "TFile.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include <TROOT.h>
#include <TStyle.h>
#include "TLine.h"

#include <iostream>
#include <vector>
#include <sstream>

using namespace std;

int main()
{

  const char *dirname="../";
  TMultiGraph *mg = new TMultiGraph();
  TMultiGraph *rmsmg = new TMultiGraph();
  TMultiGraph *fitmeanmg = new TMultiGraph();
  TMultiGraph *fitrmsmg = new TMultiGraph();

  int count_files= -1;


  TString epsname = "graphs.eps";
  TCanvas* can = new TCanvas("can", "can", 600, 700); 
  can->cd();
  can->Print(epsname+"[");
  
  vector<TH1F> hcal_detid_energies;
  vector<TH1F> ecal_detid_energies;

  vector<TH1D> slices_0T;
  vector<TH1D> slices_4T;
  vector<TH1D> slices_0T_fast;
  vector<TH1D> slices_4T_fast;

   TSystemDirectory dir(dirname, dirname);
   TList *fileslist = dir.GetListOfFiles();
   if (fileslist) {
      TSystemFile *systemfile;
      TString fname;
      TIter next(fileslist);
      while ((systemfile=(TSystemFile*)next())) {
         fname = systemfile->GetName();
         if (!systemfile->IsDirectory() && (fname.EndsWith("T.root") || fname.EndsWith("T_fast.root"))) {
           cout<<fname<<endl;
	   count_files+=1;

	   TFile* file = new TFile(dirname+fname,"READ");

	   //TH2F* h_fit  = (TH2F*) file->Get("genenergy_energy"); 	   
	   //TH2F* h_fit  = (TH2F*) file->Get("genenergy_relativCalenergy"); 	   
           //TH2F* h_fit  = (TH2F*) file->Get("ecal_detId_genenergy");
           TH2F* h_fit  = (TH2F*) file->Get("hcal_detId_genenergy");

	   //TH1F* ecal = 	   

	   ecal_detid_energies.push_back(*((TH1F*)file->Get("ecal_detId_energy_response")));	     
	   hcal_detid_energies.push_back(*((TH1F*)file->Get("hcal_detId_energy_response")));	

	   fname.Prepend("Fit_");
	   TFile* new_file = new TFile(fname,"RECREATE");

	   TObjArray* aSlice     = new TObjArray;
	   h_fit->FitSlicesX(0,0,-1,50,"QNG5",aSlice);

	   TH1D* hSlice        = (TH1D*) aSlice->At(0);
	   TH1D* hSlice_mean   = (TH1D*) aSlice->At(1);
	   TH1D* hSlice_rms    = (TH1D*) aSlice->At(2);
	   TH1D* hSlice_chifit = (TH1D*) aSlice->At(3);

	   hSlice->Write("Fit",TObject::kOverwrite);
	   hSlice_mean->Write("FitMeans",TObject::kOverwrite);
	   hSlice_rms->Write("FitRMS",TObject::kOverwrite);
	   hSlice_chifit->Write("FitChi",TObject::kOverwrite);
	   
	   new_file->mkdir("single_projections");
	   new_file->cd("single_projections");

	   double multi = 10;


	   //cout<<h_fit->GetXaxis()->GetNbins()/multi<<" "<<h_fit->GetXaxis()->GetBinCenter(1)<<" "<<h_fit->GetXaxis()->GetBinCenter(h_fit->GetXaxis()->GetNbins()+1)<<endl;;

	   TH2D* mean = new TH2D("mean","Mean",h_fit->GetXaxis()->GetNbins()/multi,h_fit->GetXaxis()->GetBinCenter(1),h_fit->GetXaxis()->GetBinCenter(h_fit->GetXaxis()->GetNbins()+1),100,-1,1);
	   TH2D* rms  = new TH2D("rms","RMS",h_fit->GetYaxis()->GetNbins()/multi,0.5,h_fit->GetYaxis()->GetNbins()/multi+0.5,100,0,1);

	   TH1D* mean_values = new TH1D("mean_values","Mean Values",100,-1,1);


	  
	   Double_t x[int(h_fit->GetXaxis()->GetNbins()/multi)], y[int(h_fit->GetXaxis()->GetNbins()/multi)], rmsy[int(h_fit->GetXaxis()->GetNbins()/multi)];
	   Double_t fit_mean[int(h_fit->GetXaxis()->GetNbins()/multi)], fit_rms[int(h_fit->GetXaxis()->GetNbins()/multi)];
	   //Double_t x[], y[], ey[];


	   //cout<<int(h_fit->GetXaxis()->GetNbins()/multi)<<" "<<sizeof(x)/sizeof(Double_t) <<endl;

	   for(double i = 1; i< h_fit->GetXaxis()->GetNbins(); i+=multi){
	     TH1D* single_projection = h_fit-> ProjectionY("projection",i,i+multi);
	     
       	     //if(single_projection->GetEntries()<100) cout<<i<<endl;
	     
	     //cout<<" Entries " <<single_projection->GetEntries()<<endl;
	     //if(single_projection->GetEntries()<100)continue;
	     TF1* f1 = new TF1("gaus","gaus",-0.8,1);

	     single_projection->Fit(f1,"WMQER");
	     

	     //cout<<f1->GetParameter(0)<<" "<<f1->GetParameter(1)<<" "<<f1->GetParameter(2)<<endl;
	     
	     //single_projection->Draw();
	     //can->Print(epsname);


	     single_projection->Write();
	    
	     double bincenter = (h_fit->GetXaxis()->GetBinCenter(i)+h_fit->GetXaxis()->GetBinCenter(i+multi))/2;
	    

	     cout<<h_fit->GetXaxis()->GetBinCenter(i)<<endl;

	     //cout<<i<<" "<<int(i/multi)<<" "<<bincenter<<" "<< single_projection->GetMean()<<" "<<single_projection->GetRMS()<<endl;
 
	     mean->Fill(bincenter,single_projection->GetMean());
	     //mean->SetBinError(i/multi,single_projection->GetRMS());
	     rms->Fill(i/multi,single_projection->GetRMS());
	    
	     //cout<<int(i/multi)<<"/"<<int((h_fit->GetYaxis()->GetNbins()+1-multi)/multi)<<endl;
	     //TString title = h_fit->GetXaxis()->GetBinCenter(i)+"<"+bincenter+"<"+h_fit->GetXaxis()->GetBinCenter(i+multi)


	     stringstream title;
	     title<< h_fit->GetXaxis()->GetBinCenter(i);
	     title<<"<";
	     title<<bincenter;
	     title<<"<";
	     title<<h_fit->GetXaxis()->GetBinCenter(i+multi);

	     cout<<title.str().c_str()<<endl;

	     single_projection->SetTitle(title.str().c_str());

	     x[int(i/multi)]=bincenter;
	     y[int(i/multi)]= single_projection->GetMean();
	     rmsy[int(i/multi)]=single_projection->GetRMS();
	     fit_mean[int(i/multi)]=f1->GetParameter(1);
	     fit_rms[int(i/multi)]=f1->GetParameter(2);

	     mean_values->Fill(single_projection->GetMean());

	     if(fname.EndsWith("0T.root"))slices_0T.push_back(*single_projection);
	     if(fname.EndsWith("4T.root"))slices_4T.push_back(*single_projection);
	     if(fname.EndsWith("0T_fast.root"))slices_0T_fast.push_back(*single_projection);
	     if(fname.EndsWith("4T_fast.root"))slices_4T_fast.push_back(*single_projection);



	     //mean_values->SetBinError(i/multi,single_projection->GetRMS());
	   }
	   

	   //cout<<h_fit->GetXaxis()->GetNbins()/multi-2<<endl;

	   TGraphErrors* fitmeanGraph = new TGraphErrors(h_fit->GetXaxis()->GetNbins()/multi,x,fit_mean,0,0);
	   fitmeanGraph->SetTitle(fname+ "Fitted Mean");
	   if(fname.Contains("fast"))fitmeanGraph->SetMarkerColor(kRed);
	   if(fname.Contains("0"))fitmeanGraph->SetMarkerStyle(21);
	   if(fname.Contains("4"))fitmeanGraph->SetMarkerStyle(24);
	   fitmeanmg->Add(fitmeanGraph);

	   TGraphErrors* fitrmsGraph = new TGraphErrors(h_fit->GetXaxis()->GetNbins()/multi,x,fit_rms,0,0);
	   fitrmsGraph->SetTitle(fname+" Sigma");
	   if(fname.Contains("fast"))fitrmsGraph->SetMarkerColor(kRed);
	   if(fname.Contains("0"))fitrmsGraph->SetMarkerStyle(21);
	   if(fname.Contains("4"))fitrmsGraph->SetMarkerStyle(24);
	   fitrmsmg->Add(fitrmsGraph);


	   TGraphErrors* errorGraph = new TGraphErrors(h_fit->GetXaxis()->GetNbins()/multi,x,y,0,0);
	   errorGraph->SetTitle(fname+" Mean");
	   if(fname.Contains("fast"))errorGraph->SetMarkerColor(kRed);
	   if(fname.Contains("0"))errorGraph->SetMarkerStyle(21);
	   if(fname.Contains("4"))errorGraph->SetMarkerStyle(24);
	   mg->Add(errorGraph);

	   TGraphErrors* rmsGraph = new TGraphErrors(h_fit->GetXaxis()->GetNbins()/multi,x,rmsy,0,0);
	   rmsGraph->SetTitle(fname+" RMS");
	   if(fname.Contains("fast"))rmsGraph->SetMarkerColor(kRed);
	   if(fname.Contains("0"))rmsGraph->SetMarkerStyle(21);
	   if(fname.Contains("4"))rmsGraph->SetMarkerStyle(24);
	   rmsmg->Add(rmsGraph);
	   

	   //errorGraph->SetMarkerStyle(20);
	   errorGraph->Draw("AP");
	   can->Print(epsname);
	   //mean->Draw();
	   //can->Print(epsname);
	   //rms->Draw();
	   //can->Print(epsname);
	   mean_values->Draw();
	   can->Print(epsname);


	   errorGraph->Write();

	   mean->Write();
	   rms->Write();
	   mean_values->Write();
	 
	   /*
	   delete hSlice;
	   delete hSlice_mean;
	   delete hSlice_rms;
	   delete hSlice_chifit;
	   delete mean;
	   delete rms;
	   delete mean_values;
	   //delete errorGraph;
     
	   delete aSlice;
 	   delete file;
	   delete new_file;
	   */
         }
      }
   }

 
   //cout<<slices_0T.size()<<" "<<slices_4T.size()<<" "<<slices_0T_fast.size()<<" "<<slices_4T_fast.size()<<endl;

   //can->SetLogy();
   for(unsigned int i =0;i<slices_0T.size(); ++i){



     TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
     pad1->SetBottomMargin(0.05);
     //pad1->SetLogy();
     pad1->Draw();
     TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
     pad2->SetTopMargin(0.05);
   
     pad2->Draw();
     
     pad1->cd();
     

     gStyle->SetOptStat(kFALSE);
     int divide = 25;

     slices_0T[i].Rebin(divide);
     slices_4T[i].Rebin(divide);
     slices_0T_fast[i].Rebin(divide);
     slices_4T_fast[i].Rebin(divide);

     double maximum = slices_0T[i].GetMaximum();
     if(maximum<slices_4T[i].GetMaximum()) maximum=slices_4T[i].GetMaximum();
     if(maximum<slices_0T_fast[i].GetMaximum()) maximum=slices_0T_fast[i].GetMaximum();
     if(maximum<slices_4T_fast[i].GetMaximum()) maximum=slices_4T_fast[i].GetMaximum();

     slices_0T[i].SetMarkerColor(kBlack);
     slices_0T[i].SetMarkerStyle(21);
     slices_4T[i].SetMarkerColor(kBlack);
     slices_4T[i].SetMarkerStyle(24);

     slices_0T_fast[i].SetMarkerColor(kRed);
     slices_0T_fast[i].SetMarkerStyle(21);
     slices_4T_fast[i].SetMarkerColor(kRed);
     slices_4T_fast[i].SetMarkerStyle(24);
     
     slices_0T[i].SetMaximum(maximum+20);
     
     slices_0T[i].Draw("P");
     slices_4T[i].Draw("PSAME");
     slices_0T_fast[i].Draw("PSAME");
     slices_4T_fast[i].Draw("PSAME");

     //can->cd();
     pad2->cd();

      
     TH1F* h_0T_help = (TH1F*)slices_0T[i].Clone("0T_hist");
     TH1F* h_4T_help = (TH1F*)slices_0T[i].Clone("4T_hist");


     h_0T_help->Divide(&slices_0T_fast[i]);
     h_4T_help->Divide(&slices_4T_fast[i]);
     h_0T_help->SetMarkerColor(kRed);
     h_4T_help->SetMarkerColor(kRed);
     h_0T_help->SetMarkerStyle(21);
     h_4T_help->SetMarkerStyle(24);


     TLine* line = new TLine( slices_0T[i].GetXaxis()->GetXmin(),1,slices_0T[i].GetXaxis()->GetXmax(),1); 

     //maximum = slices_0T[i].GetMaximum();
     //if(maximum<slices_4T[i].GetMaximum()) maximum=slices_4T[i].GetMaximum();

     //cout<<maximum<<endl;

     h_0T_help->SetMaximum(8.1);

     h_0T_help->Draw("P");
     h_4T_help->Draw("PSAME");

     line->Draw();
     can->cd();

     can->Print(epsname);
   }
   //can->SetLogy(0);



   mg->SetTitle("Mean");
   mg->Draw("ap");
   can->Print(epsname);
   rmsmg->SetTitle("RMS");
   rmsmg->Draw("ap");
   can->Print(epsname);
   fitmeanmg->SetTitle("Fitted Mean");
   fitmeanmg->Draw("ap");
   can->Print(epsname);
   fitrmsmg->SetTitle("Sigma");
   fitrmsmg->Draw("ap");
   can->Print(epsname);

   
   //cout<<hcal_detid_energies.size()<<endl;

   can->SetLogy();

   TH1F sum_detid_0T = hcal_detid_energies[0];
   sum_detid_0T.Add(&ecal_detid_energies[0]);
   
   TH1F sum_detid_4T = hcal_detid_energies[1];
   sum_detid_4T.Add(&ecal_detid_energies[1]);

   TH1F sum_detid_4T_fast = hcal_detid_energies[2];
   sum_detid_4T_fast.Add(&ecal_detid_energies[2]);

   TH1F sum_detid_0T_fast = hcal_detid_energies[3];
   sum_detid_0T_fast.Add(&ecal_detid_energies[3]);

   sum_detid_0T.SetTitle("Sum DetId Energy");


   int rebin_detid =25;

   sum_detid_0T.Rebin(rebin_detid);
   sum_detid_4T.Rebin(rebin_detid);
   sum_detid_4T_fast.Rebin(rebin_detid);
   sum_detid_0T_fast.Rebin(rebin_detid);


   sum_detid_0T.SetMarkerStyle(21);
   sum_detid_4T.SetMarkerStyle(24);
   sum_detid_0T_fast.SetMarkerStyle(21);
   sum_detid_4T_fast.SetMarkerStyle(24);
   sum_detid_0T_fast.SetMarkerColor(kRed);
   sum_detid_4T_fast.SetMarkerColor(kRed);

   sum_detid_0T.Draw("P"); 
   sum_detid_4T.Draw("PSAME");
   sum_detid_4T_fast.Draw("PSAME");
   sum_detid_0T_fast.Draw("PSAME");
   can->Print(epsname);
 
   can->SetLogy(0);
  
   sum_detid_0T.SetTitle("Ratio DetId Energy");

   sum_detid_0T.Divide(&sum_detid_0T_fast); 
   sum_detid_4T.Divide(&sum_detid_4T_fast);

   
   sum_detid_4T.Draw("P");
   sum_detid_0T.Draw("PSAME"); 
   
   can->Print(epsname);




   /*
   for(unsigned int i=0; i< hcal_detid_energies.size();++i){
     
     hcal_detid_energies;
     ecal_detid_energies;
   }
   */

   
   


   can->Print(epsname+"]");


   return 0;
}
