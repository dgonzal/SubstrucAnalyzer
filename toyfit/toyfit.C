// c++ classes
#include <iostream>

//RooFit
#include "RooRealVar.h"
#include "RooGamma.h"
#include "RooDataSet.h"
#include "RooPlot.h"

//ROOT 
#include "TCanvas.h"
using namespace std;

int main(){
  // --- Observable --- 
  RooRealVar depth("depth","depth (cm)",0,0,20) ; 
  
  // --- Gamma pdf --- 
  RooRealVar alpha("alpha","alpha",5,0,100);
  alpha.removeMax(); // set infinite range
  RooRealVar beta("sigwidth","beta",1,0,100) ;  
  beta.removeMax(); // set infinite range
  RooRealVar mu("mu","mu",5); // no ranges means variable is fixed in fit
  RooGamma gamma("gamma","gamma pdf",depth,alpha,beta,mu) ; 

  RooDataSet * data = gamma.generate(depth,100);
  gamma.fitTo(*data);

  cout << alpha.getVal() << " " << beta.getVal() << std::endl;

  RooPlot * depthFrame = depth.frame(mu.getVal(),depth.getMax(),20); 
  data->plotOn(depthFrame);
  gamma.plotOn(depthFrame);

  TCanvas * can = new TCanvas("can", "can", 600, 600); 
  can->cd();
  depthFrame->Draw();
  can->Print("test.eps");
}


