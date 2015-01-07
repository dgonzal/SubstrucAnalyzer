//ROOT includes
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
#include "TTree.h"

//RooFit includes
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGamma.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooHist.h"
#include "RooMCStudy.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "TMatrixDSym.h"
#include "RooGenericPdf.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooFitResult.h"
//#include "RooDataSet.h"

//MyCLass
#include "RooDoubleSidedCball.h"



//c++ libs
#include <algorithm>    
#include <iomanip>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <math.h>  


using namespace std;
using namespace RooFit;


struct fitParams{
  int status;
  double energy;
  double eta;
  double chi;// chi2/ndof
  double mu;
  double sigma;
  double rms;
  double aL;
  double nL;
  double aR;
  double nR;
};




bool eta_energy_sort(fitParams x1, fitParams x2){
  if(round(x1.eta*10)<round(x2.eta*10))return true;
  else if(x1.energy < x2.energy && round(x1.eta*10)==round(x2.eta*10)) return true;

  return false;
};
bool energy_eta_sort(fitParams x1, fitParams x2){
     
  if(x1.energy+ x1.eta/100 <x2.energy+x2.eta/100) return true;
 
   
  return false;
};


int main()
{

  const char *dirname="/nfs/dust/cms/user/gonvaq/fastsim/fullsim/no_magneticField/";
  
  //TFile* fit_results = new TFile("FitResults","RECREATE");

  vector<fitParams> parameter_list;

  TString epsname = "ResponseFitResults.eps";
  TCanvas* can = new TCanvas("can", "can", 600, 700); 
  can->cd();
  can->Print(epsname+"[");
  
   /*
  RooRealVar test_range("test_range","test_range",-10,4);
  RooRealVar test_mu("test_mu","test_mu",0) ;
  RooRealVar test_sigma("test_sigma","test_sigma",1);
  RooRealVar test_aL("test_aL","test_aL",1);
  RooRealVar test_aR("test_aR","test_aR",1);
  RooRealVar test_nR("test_nR","test_nR",1.01);
  RooRealVar test_nL("test_nL","test_nL",1.01);
 
  RooDoubleSidedCball crystalball_check("crystalball_check","crystal ball check",test_range,test_mu,test_sigma,test_aL,test_nL,test_aR,test_nR);  
 
  // Activate debug-level messages for topic integration to be able to follow actions below
  RooMsgService::instance().addStream(DEBUG,Topic(Integration)) ;

  RooAbsReal* int_crystalball = crystalball_check.createIntegral(test_range) ;
  Double_t val = int_crystalball->getVal() ;
  cout << " [1] int_dx crystalball(x) = " << setprecision(15) << val << endl ;

  
  RooPlot* testframe = test_range.frame();
  crystalball_check.plotOn(testframe) ;
  testframe->Draw();
  can->Print(epsname);
  */

  TSystemDirectory dir(dirname, dirname);
  TList *fileslist = dir.GetListOfFiles();
  if (fileslist) {
    TSystemFile *systemfile;
    TString fname;
    TIter next(fileslist);
    while ((systemfile=(TSystemFile*)next())) {
      fname = systemfile->GetName();
      if (!systemfile->IsDirectory() && fname.EndsWith(".root") && fname.Contains("response_")&& fname.Contains("")){//eta_2.8_2.9_E_9.0 eta_0.5_0.6_E_1.0
	
	//cout<<fname<<endl;
	
	
	TFile* file = new TFile(dirname+fname,"READ");
	TH1F* energy_response  = (TH1F*) file->Get("energy_response");
	TH1F* gen_eta  = (TH1F*) file->Get("gen_eta");
	TH1F* gen_energy  = (TH1F*) file->Get("gen_energy");
	
	TTree* tree = (TTree*) file->Get("Total");

	

	double mu_estimate = 0;
	double sigma_estimate =0;
	double peak_bin =energy_response->GetBinContent(1);
	int peak_position = -1;

	
	cout<<"----------------------"<<endl;
	cout<<"----------------------"<<endl;


	for(unsigned int i = 1; i<energy_response->GetNbinsX(); ++i){
	  //cout<<energy_response->GetBinContent(i)<<" "<<energy_response->GetXaxis()->GetBinCenter(i)<<endl;
	  
	  if(peak_bin<=energy_response->GetBinContent(i)){
	    //mu_estimate=energy_response->GetXaxis()->GetBinCenter(i);
	    peak_bin =energy_response->GetBinContent(i);
	    peak_position =i;
	  }
	}
	
	if(gen_energy->GetMean()<=energy_response->GetMean()){
	  int zerosAfterPeak =0;
	  int endposition = -1;
	  
	  for (unsigned int i = peak_position;i<energy_response->GetNbinsX(); ++i){
	    if(energy_response->GetBinContent(i) < 5){
	      zerosAfterPeak+=1;
	      endposition = i;
	    }
	    if(zerosAfterPeak>100) break;
	  }


	  cout <<endposition<<endl;
	  if(endposition==-1) endposition = energy_response->GetNbinsX()-2;
	
	  energy_response->GetXaxis()->SetRange(1,endposition);

	  mu_estimate = energy_response->GetMean();
	  sigma_estimate = energy_response->GetRMS();



	}
	else{
	  int minBin = energy_response->FindBin(energy_response->GetMean()- energy_response->GetRMS());
	  if(minBin<1) minBin =1;
	  int maxBin = energy_response->FindBin(energy_response->GetMean()+ energy_response->GetRMS());
	  if(maxBin<0) maxBin = energy_response->GetNbinsX()-2;

	  energy_response->GetXaxis()->SetRange(minBin,maxBin);

	  mu_estimate = energy_response->GetMean();
	  sigma_estimate = energy_response->GetRMS();


	  cout<< mu_estimate<<" "<< sigma_estimate<<" "<< minBin << " "<<maxBin<<endl;

	  if(sigma_estimate == 0) sigma_estimate = 0.1;
	}
	
	
	if(sigma_estimate!=sigma_estimate) assert(0==1);


	double xlow,xhigh;
	
	xlow = energy_response->GetXaxis()->GetXmin();
	xhigh =  9;

	

	if(gen_energy->GetMean()>=5 && gen_energy->GetMean()<30){
	  //energy_response->Rebin(10);
	  xhigh = 55;
	}
	//else if(gen_energy->GetMean()>=1 && gen_energy->GetMean()<5)
	  //energy_response->Rebin(5);
	else if(gen_energy->GetMean()>=30&&gen_energy->GetMean()<90){
	  //energy_response->Rebin(20);
	  xlow = 15;
	  xhigh = 120;
	}
	else if(gen_energy->GetMean()>90){
	  //energy_response->Rebin(50);
	  xlow = 60;
	  xhigh = energy_response->GetXaxis()->GetXmax();
	}


	cout<<"----------------------"<<endl;
	cout<<"----------------------"<<endl;

	//firstfit
	//RooRealVar first_x("first_x","first_x",mu_estimate,mu_estimate-sigma_estimate>0?mu_estimate-sigma_estimate:0,mu_estimate+sigma_estimate<energy_response->GetXaxis()->GetXmax()?mu_estimate+sigma_estimate:energy_response->GetXaxis()->GetXmax()));
      /*
	RooRealVar first_mu("mu","mu",mu_estimate,mu_estimate-sigma_estimate>0?mu_estimate-sigma_estimate:0,mu_estimate+sigma_estimate<energy_response->GetXaxis()->GetXmax()?mu_estimate+sigma_estimate:energy_response->GetXaxis()->GetXmax());
	RooRealVar first_sigma("sigma","sigma",sigma_estimate,0,3000);

	RooGaussian gaus("gaus","gaus",first_x,first_mu,first_sigma);
	gaus.fitTo(roodatahist, Strategy(2));
      */


	//assert( energy_response->GetMean()<9);
	
	RooRealVar e_calo("e_calo","e_calo",0,1000);
	e_calo.removeMax();

	RooRealVar x("x","x",0,5000);
	RooRealVar mu("mu","mu",mu_estimate,energy_response->GetXaxis()->GetXmin(),energy_response->GetXaxis()->GetXmax()) ;
     	RooRealVar sigma("sigma","sigma",sigma_estimate,0.01,1500);

	double aL_in = 2.6;
	double aR_in = 1.;
	double nR_in = 5.;
	double nL_in = 1.1;

	if( fname.Contains("eta_2.4_2.5_E_30_")||fname.Contains("eta_2.3_2.4_E_50_")){
	  aL_in = 1.1;
	  aR_in = 1.1;
	  nL_in = 1.1;
	  nR_in = 1.1;
	}
	
	RooRealVar aL("aL","aL",aL_in,0.001,50.) ;
	RooRealVar aR("aR","aR",aR_in,0.001,50.) ;
	RooRealVar nR("nR","nR",nR_in,1.001,300.) ;
	RooRealVar nL("nL","nL",nL_in,1.001,250.) ;
	

	RooDataHist roodatahist("roodatahist","roodatahist",x,RooFit::Import(*energy_response));
	RooDataSet unbinneddata("roodataset","energy calorimeter",e_calo,Import(*tree)) ; 


	// PDF
	RooDoubleSidedCball crystalball("crystalball","crystal ball",e_calo,mu,sigma,aL,nL,aR,nR); 

	RooFitResult* fitResult = crystalball.fitTo(unbinneddata, Strategy(2), Minimizer("MINUIT","MIGRADIMPROVED"),Save(true) );
      
	//if(fitResult->status()=!0)
	//  int a=1;

	//Plot PDF and toy data overlaid
	
	
	RooPlot* single_fit = e_calo.frame(xlow,xhigh);
	single_fit->SetTitle(fname);
	//roodatahist.plotOn(single_fit) ;
	unbinneddata.plotOn(single_fit) ;
	crystalball.plotOn(single_fit) ;
	//gaus.plotOn(single_fit) ;
	single_fit->Draw();


	struct fitParams params = {fitResult->status(),gen_energy->GetMean(),gen_eta->GetMean(),single_fit->chiSquare(5),mu.getVal(),sigma.getVal(),energy_response->GetRMS(),aL.getVal(),nL.getVal(),aR.getVal(),nR.getVal()};
	parameter_list.push_back(params);
	
	
	//energy_response->SetTitle(fname);

	//energy_response->Draw();
	can->Print(epsname);
     
      
      }
    }
  }

  can->Print(epsname+"]");


  sort(parameter_list.begin(),parameter_list.end(),energy_eta_sort);
  
  ofstream paramFile;
  paramFile.open("paramFile.txt");

  for(unsigned int i =0; i< parameter_list.size();++i){
    paramFile<<setprecision(2)<<fixed<<"E "<<parameter_list[i].energy<<" fit "<<parameter_list[i].status<<" eta "<<parameter_list[i].eta<<scientific<<" chi2 "<<parameter_list[i].chi<<" mu "<<parameter_list[i].mu<<" RMS "<<parameter_list[i].rms<<" sigma "<<parameter_list[i].sigma<<" aL " <<parameter_list[i].aL<<" nL " <<parameter_list[i].nL<<" aR " <<parameter_list[i].aR<<" nR "<<parameter_list[i].nR<<endl;
  }

  paramFile.close();

  
  //sort(parameter_list.begin(),parameter_list.end(),eta_energy_sort);

  
  ofstream resultFile;
  resultFile.open("resultFile.txt");

  char* variable_array[] ={"mu" ,"sigma","aL","nL","aR","nR"}; 

  
  
  char* detector_part[] ={"HB","HE"};


  for(unsigned int m =0; m<2; ++m){
    int count_var =0;

    for (char* &var_name: variable_array)
      {
	//std::cout << var_name << std::endl;
	count_var+=1;

	int count_hb =0;
	
	resultFile<<var_name<<"_"<<detector_part[m]<<"  = cms.vdouble( *["<<endl;

	for(unsigned int i =0; i< parameter_list.size();++i){
	  
	  if(m == 1 && parameter_list[i].eta<1.6){
	    continue;
	  }
	  else if(m == 0 && parameter_list[i].eta>1.6){
	    continue;
	  }
	  
	  if(m==0) count_hb+=1;

	  switch (count_var){
	  case 1:
	    resultFile<<parameter_list[i].mu;
	    break;
	  case 2:
	    resultFile<<parameter_list[i].sigma;
	    break;
	  case 3:
	    resultFile<<parameter_list[i].aL;
	    break;
	  case 4:
	    resultFile<<parameter_list[i].nL;
	    break;
	  case 5:
	    resultFile<<parameter_list[i].aR;
	    break;
	  case 6:
	    resultFile<<parameter_list[i].nR;
	    break;
	  }
	  resultFile<<", ";
	  //resultFile<<"  #"<<parameter_list[i].energy <<endl; 
	  if(((count_hb)%16==0 && m ==0) || (i+1)%30==0 )resultFile<<"  #"<<parameter_list[i].energy <<endl; 	
	  
	}
	resultFile<< "]), "<<endl<<endl;
      }
  }

  resultFile.close();

  return 0;
}

