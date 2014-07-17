// c++ classes
#include <iostream>

//RooFit
#include "RooRealVar.h"
#include "RooGamma.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooHist.h"
//#include "RooCmdArg.h"

//ROOT 
#include "TCanvas.h"
#include "TGaxis.h"
#include "TLegend.h"

using namespace std;
using namespace RooFit;

int main(){

  RooRealVar zero("zero","",0,0,20); 

  // --- Observable --- 
  RooRealVar depth("depth","depth (cm)",0,0,20) ; 
  
  // --- Gamma pdf --- 
  RooRealVar alpha("alpha","alpha",5.,1.,100.);
  alpha.removeMax(); // set infinite range
  RooRealVar beta("beta","beta",1.,0.,100.) ;  
  beta.removeMax(); // set infinite range
  RooRealVar mu("mu","mu",1); // no ranges means variable is fixed in fit

  RooRealVar lambda("lambda","lambda",3.,1.,100.);
  lambda.removeMax(); // set infinite range
  RooRealVar theta("theta","theta",1.5,0.,100.) ;  
  theta.removeMax(); // set infinite range

  RooGamma gamma_one("gamma_one","gamma_one pdf",depth,alpha,beta,mu) ; 
  RooGamma gamma_two("gamma_two","gamma_two pdf",depth,lambda,theta,mu) ; 

  RooRealVar fsig("fsig","signal fraction",0.5,0.,1.) ;

  RooAddPdf model("model","model",RooArgList(gamma_one,gamma_two),fsig) ;
  RooDataSet * data = model.generate(depth,100);
  fsig.setVal(0.1);
  model.fitTo(*data);

  RooPlot * depthFrame = depth.frame(mu.getVal(),depth.getMax(),20);
  zero.plotOn(depthFrame,LineStyle(kDashed),LineColor(kBlack));
  data->plotOn(depthFrame);
  model.plotOn(depthFrame,Components(gamma_one),LineStyle(kDashed),LineColor(kGreen),Name("gamma_one"));
  model.plotOn(depthFrame,Components(gamma_two),LineStyle(kDashed),LineColor(kRed),Name("gamma_two"));
  model.plotOn(depthFrame,Name("model"));
  

  cout << alpha.getVal() << " " << beta.getVal() << std::endl;
  cout << lambda.getVal() << " " << theta.getVal() << std::endl;
  cout << fsig.getVal() << " "<<depthFrame->chiSquare()<< std::endl;


  RooHist* hpull = depthFrame->pullHist();
  hpull->SetLineColor(0);
  hpull->SetMarkerColor(kRed);
  hpull->SetMarkerStyle(21);

  depthFrame->addPlotable(hpull);
  depthFrame->SetMinimum(-5);

  RooPlot * pullhistFrame = depth.frame(Title("Pull Distribution"));
  pullhistFrame->addPlotable(hpull);
  pullhistFrame->SetMaximum(5);
  pullhistFrame->SetMinimum(-5);
  

  TCanvas * can = new TCanvas("can", "can", 600, 600); 
  can->cd();
  /*
  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  TPad *pad2 = new TPad("pad2","",0,0,1,1);
  pad2->SetFillStyle(4000); //will be transparent
  */
 
  TLegend * legend = new TLegend(0.55,0.65,.9,0.9);
  legend->SetTextFont(72);
  legend->SetTextSize(0.04);
  legend->SetFillColor(kWhite);
  legend->AddEntry(depthFrame->findObject("gamma_one"),"gamma one");
  legend->AddEntry(depthFrame->findObject("gamma_two"),"gamma two");
  legend->AddEntry(depthFrame->findObject("model")," gamma");

  can->Print("test.ps[");

 
  depthFrame->Draw();
  legend->Draw();
  can->Print("test.ps");

  //pullhistFrame->Draw();

  //can->Print("test.ps");

  can->Print("test.ps]");
 

  return 0;
}


