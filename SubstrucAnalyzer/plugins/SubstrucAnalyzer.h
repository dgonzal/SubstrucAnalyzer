#ifndef SubstrucAnalyzer_H
#define SubstrucAnalyzer_H
//
// Original Author:  Daniel Gonzalez
//         Created:  Thu, 03 Jul 2014 11:54:05 GMT
//
//

// system include files
#include <memory>
#include <string>

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


// ROOT
#include "TH1F.h"
#include "TFile.h"

//
// class declaration
//

class SubstrucAnalyzer : public edm::EDAnalyzer {

public:
  SubstrucAnalyzer(const edm::ParameterSet& ps);
  virtual ~SubstrucAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------

  std::vector<edm::InputTag> jetLabels_;
  std::vector<double> jetPtMins_;
  //edm::InputTag jet_match_ ;
  double matching_radius_;
  edm::EDGetTokenT<reco::GenParticleCollection>  genParticle_;
  std::vector<const reco::GenParticle*> gentop;
  std::vector<double> gentop_eta;
  std::vector<double> gentop_phi;
  std::string OutputName_; 

  std::vector< edm::EDGetTokenT< edm::View<reco::Jet> > > jetTokens_;

  // Histograms with all jets 
  std::vector<TH1F*> jet_N;
  std::vector<TH1F*> jet_pt;
  std::vector<TH1F*> jet_eta;
  std::vector<TH1F*> jet_phi;
  std::vector<TH1F*> jet_mass;
  
  //Histograms where the jets are matched to the top quarks 
  std::vector<TH1F*> jet_matched_N;
  std::vector<TH1F*> jet_matched_pt;
  std::vector<TH1F*> jet_matched_eta;
  std::vector<TH1F*> jet_matched_phi;
  std::vector<TH1F*> jet_matched_mass;

  //Histograms with the subjets of any Jets if they exist
  std::vector<TH1F*> subjet_N;
  std::vector<TH1F*> subjet_pt;
  std::vector<TH1F*> subjet_eta;
  std::vector<TH1F*> subjet_phi;
  std::vector<TH1F*> subjet_mass;

  //-------------------------
  //Histograms for the cmsTopTag 
  //at some point these could be used for additional TopTaggers
  std::vector<TH1F*> Topjet_N;
  std::vector<TH1F*> Topjet_pt;
  std::vector<TH1F*> Topjet_eta;
  std::vector<TH1F*> Topjet_phi;
  std::vector<TH1F*> Topjet_mass;

  //Top matched in DeltaR to one of the top quarks
  std::vector<TH1F*> Topjet_matched_N;
  std::vector<TH1F*> Topjet_matched_pt;
  std::vector<TH1F*> Topjet_matched_eta;
  std::vector<TH1F*> Topjet_matched_phi;
  std::vector<TH1F*> Topjet_matched_mass;

  //Subjets of the matched Topjets
  std::vector<TH1F*> Topsubjet_matched_N;
  std::vector<TH1F*> Topsubjet_matched_pt;
  std::vector<TH1F*> Topsubjet_matched_eta;
  std::vector<TH1F*> Topsubjet_matched_phi;
  std::vector<TH1F*> Topsubjet_matched_mass;

  //Mistag rate in case a Jet is not near a top quark
  std::vector<TH1F*> Topjet_unmatched_N;
  std::vector<TH1F*> Topjet_unmatched_pt;
  std::vector<TH1F*> Topjet_unmatched_eta;
  std::vector<TH1F*> Topjet_unmatched_phi;
  std::vector<TH1F*> Topjet_unmatched_mass;

 
};
#endif
