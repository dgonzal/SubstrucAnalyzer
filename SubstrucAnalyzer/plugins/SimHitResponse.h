#ifndef SimHitResponse_H
#define SimHitResponse_H
//
// Original Author:  Daniel Gonzalez
//         Created:  Thu, 03 Jul 2014 11:54:05 GMT
//
//

//BOOST libs
#include "boost/tuple/tuple.hpp"

// system include files
#include <memory>
#include <string>
#include <list>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// ParticleFlow
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"

// Jets
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"

//GenParticle
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//Substructure
#include "RecoJets/JetAlgorithms/interface/CATopJetHelper.h"
#include "DataFormats/JetReco/interface/CATopJetTagInfo.h"

//Calculations
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/CandUtils/interface/pdgIdUtils.h"

//Geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

//DataFormats 
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"


// ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"



//
// class declaration
//

class SimHitResponse : public edm::EDAnalyzer {

public:
  SimHitResponse(const edm::ParameterSet& ps);
  virtual ~SimHitResponse();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  edm::EDGetTokenT<edm::PCaloHitContainer> ECAL_Hits_Label_;
  edm::EDGetTokenT<edm::PCaloHitContainer> Hits_Label_;
  std::string OutputName_; 
  std::string Module_;



  //TH1F* henergy_sum;

  TH1F* gen_energy;
  TH1F* gen_pT;

  TH1F* eta_response;
  TH1F* phi_response;
  TH1F* energy_response;

  TH1F* eb_energy;
  TH1F* ee_energy;
  TH1F* es_energy;
  TH1F* h_energy;

  TH2F* scatter_dR_depth;
  TH2F* scatter_dR_pT;
  TH2F* scatter_dR_pT_e;

  TH2F* scatter_pT_hits;
  TH2F* scatter_eta_hits;
  TH2F* scatter_pT_hitenergy;

  TH2F* eb_eta_hits;
  TH2F* ee_eta_hits;
  TH2F* es_eta_hits;
  TH2F* hcal_eta_hits;
  
  TH1F* eb_energy_response;
  TH1F* ee_energy_response;
  TH1F* es_energy_response;
  TH1F* hcal_energy_response;
  TH1F* hcal_weight_energy_response;
  TH1F* hcal_detId_energy_response;
  TH1F* ecal_detId_energy_response;

  TH2F* eb_pt_hits;
  TH2F* ee_pt_hits;
  TH2F* es_pt_hits;
  TH2F* hcal_pt_hits;
  /*
  TH2F* eb_pt_e;
  TH2F* ee_pt_e;
  TH2F* es_pt_e;
  TH2F* hcal_pt_e;
  */
  TH2F* scatter_genpT_hitenergy;
  TH2F* scatter_genenergy_hitenergy;

  TH2F* scatter_genpT_energy;   
  TH2F* scatter_genenergy_energy;

  TH2F* scatter_x_z;
  TH2F* scatter_y_z;
 
  TH3F* scatter_x_y_z; 


  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------

};
#endif
