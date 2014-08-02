
// -*- C++ -*-
//
// Package:    Substructure/SubstrucAnalyzer
// Class:      SubstrucAnalyzer
// 
/**\class SubstrucAnalyzer SubstrucAnalyzer.cc Substructure/SubstrucAnalyzer/plugins/SubstrucAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//


#include "Substructure/SubstrucAnalyzer/plugins/SubstrucAnalyzer.h"


using namespace edm;
using namespace std;
using namespace reco; 




SubstrucAnalyzer::SubstrucAnalyzer(const edm::ParameterSet& ps)
{
  //jet_match_ =ps.getParameter<edm::InputTag>("jets_to_match");


  //get the parameters from the configuration file 
  //and do what ever initialization is needed
  jetLabels_ = ps.getParameter<std::vector<edm::InputTag> >("jetLabels");
  for ( std::vector<edm::InputTag>::const_iterator jetlabel = jetLabels_.begin(),
	  jetlabelEnd = jetLabels_.end(); jetlabel != jetlabelEnd; ++jetlabel ) {
    jetTokens_.push_back( consumes<edm::View<reco::Jet>>( *jetlabel ) );
  }

  OutputName_ = ps.getParameter<string>("OutputName");

  genParticle_ = consumes<reco::GenParticleCollection>(ps.getParameter<edm::InputTag>("GenParticle"));

  jetPtMins_ = ps.getParameter<std::vector<double> > ("jetPtMins");
  matching_radius_ = ps.getParameter<double>  ("matching_radius");
  
 
}


SubstrucAnalyzer::~SubstrucAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SubstrucAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle<reco::GenParticleCollection> genParticles;
  
  iEvent.getByToken(genParticle_, genParticles );

  //find all topquarks in the GenParticleCollection
  for (reco::GenParticleCollection::const_iterator iter=genParticles->begin();iter!=genParticles->end();++iter){
    if((&*iter)->pdgId() == 6 || (&*iter)->pdgId() == -6) {
      gentop.push_back(&*iter);
      gentop_eta.push_back(iter->eta());
      gentop_phi.push_back(iter->phi());
    }
  }
 
  
  double deltaR_value = matching_radius_*matching_radius_;

  //Number of Jets that are matched/unmatched or have a cmsTopTag
  int top_counter = 0;
  int top_matched_counter = 0;
  int top_unmatched_counter = 0;

  for (unsigned int icoll = 0; icoll < jetLabels_.size(); ++icoll ) {
    edm::Handle<edm::View<reco::Jet>>  pfJetCollection;
   
    bool ValidPFJets = iEvent.getByToken(jetTokens_[icoll], pfJetCollection );
    if(!ValidPFJets) continue;
    
    edm::View<reco::Jet> const & pfjets = *pfJetCollection;
    jet_N[icoll]->Fill(pfjets.size());

    
    for ( edm::View<reco::Jet>::const_iterator jet = pfjets.begin(),jetEnd = pfjets.end(); jet != jetEnd; ++jet ) {
      if(jet->pt()<jetPtMins_[icoll]) continue;
      jet_pt[icoll]->Fill(jet->pt());
      jet_eta[icoll]->Fill(jet->eta());
      jet_phi[icoll]->Fill(jet->phi());
      jet_mass[icoll]->Fill(jet->mass());

      reco::BasicJet const * basicjet = dynamic_cast<reco::BasicJet const *>( &*jet);

      if ( jetLabels_[icoll].label() == "cmsTopTagPFJets") {
	CATopJetHelper helper(171.2, 80.4);
	reco::CATopJetProperties properties = helper(*basicjet);

	//cmsTopTag cuts performed on any recontructed jet
	if ( jet->numberOfDaughters() > 2 && properties.minMass>50 && properties.topMass<200&& properties.topMass > 140) {
	  top_counter +=1;
	  
	  Topjet_pt[icoll]->Fill(jet->pt());
	  Topjet_eta[icoll]->Fill(jet->eta());
	  Topjet_phi[icoll]->Fill(jet->phi());
	  Topjet_mass[icoll]->Fill(jet->mass());

	  unsigned int unmatch = 0;
	  
	  //matching for Topsubjets with top quarks
	  for(unsigned int it_tops=0; it_tops<gentop_eta.size(); ++it_tops){
	    if(deltaR2(jet->eta(),jet->phi(),gentop_eta[it_tops],gentop_phi[it_tops])<deltaR_value){
	      top_matched_counter +=1;
	      Topjet_matched_pt[icoll]->Fill(jet->pt());
	      Topjet_matched_eta[icoll]->Fill(jet->eta());
	      Topjet_matched_phi[icoll]->Fill(jet->phi());
	      Topjet_matched_mass[icoll]->Fill(jet->mass());

	      Top_matched_pt[icoll]->Fill( gentop[it_tops]->pt());
	      Top_matched_eta[icoll]->Fill(gentop[it_tops]->eta());
	      Top_matched_phi[icoll]->Fill(gentop[it_tops]->phi());
	      Top_matched_mass[icoll]->Fill(gentop[it_tops]->mass());

	      if ( basicjet != 0 ) 
		Topsubjet_matched_N[icoll]->Fill(jet->numberOfDaughters() );

	      for ( unsigned int ida = 0; ida < jet->numberOfDaughters(); ++ida ) {
		reco::Candidate const * subjet = jet->daughter(ida);
		Topsubjet_matched_pt  [icoll]->Fill ( subjet->pt() );
		Topsubjet_matched_eta [icoll]->Fill ( subjet->rapidity() );
		Topsubjet_matched_phi [icoll]->Fill ( subjet->phi() );
		Topsubjet_matched_mass[icoll]->Fill ( subjet->mass() );
	      }
	    }
	    else{
	      unmatch +=1;
	    }
	  }
	  //unmatched are written out after the loop over all tops
	  if(unmatch == gentop_eta.size()){
	    top_unmatched_counter=+1;
	    Topjet_unmatched_pt[icoll]->Fill(jet->pt());
	    Topjet_unmatched_eta[icoll]->Fill(jet->eta());
	    Topjet_unmatched_phi[icoll]->Fill(jet->phi());
	    Topjet_unmatched_mass[icoll]->Fill(jet->mass());
	  }
	}
      }

      Topjet_N[icoll]->Fill( top_matched_counter);
      Topjet_matched_N[icoll]->Fill( top_counter);
      Topjet_unmatched_N[icoll]->Fill(top_unmatched_counter );


      //perform matching for all jets
      int jet_top_matched_counter = 0;
      for(unsigned int it_tops=0; it_tops<gentop_eta.size(); ++it_tops){
	if(deltaR2(jet->eta(),jet->phi(),gentop_eta[it_tops],gentop_phi[it_tops])<deltaR_value){
	  jet_top_matched_counter +=1;
	  jet_matched_pt[icoll]->Fill(jet->pt());
	  jet_matched_eta[icoll]->Fill(jet->eta());
	  jet_matched_phi[icoll]->Fill(jet->phi());
	  jet_matched_mass[icoll]->Fill(jet->mass());
	}
      }
      jet_matched_N[icoll]->Fill(jet_top_matched_counter);

      if ( basicjet != 0 ) {
	subjet_N[icoll]->Fill(jet->numberOfDaughters() );
	
	for ( unsigned int ida = 0; ida < jet->numberOfDaughters(); ++ida ) {
	  reco::Candidate const * subjet = jet->daughter(ida);
	  subjet_pt  [icoll]->Fill ( subjet->pt() );
	  subjet_eta [icoll]->Fill ( subjet->rapidity() );
	  subjet_phi [icoll]->Fill ( subjet->phi() );
	  subjet_mass[icoll]->Fill ( subjet->mass() );
	}
      }
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void SubstrucAnalyzer::beginJob(){
  for (unsigned int icoll = 0; icoll < jetLabels_.size(); ++icoll  ) {
    string hname = jetLabels_[icoll].label();


    //give names "meaningfull" names to all the different histgrams I am about to fill
    //as next step this should be put into a class/function that handels that stuff for you
    string name_jet_N = hname;
    name_jet_N.append("_jet_N");
    string name_jet_pt = hname;
    name_jet_pt.append("_jet_pT");
    string name_jet_eta = hname;
    name_jet_eta.append("_jet_eta");
    string name_jet_phi = hname;
    name_jet_phi.append("_jet_phi");
    string name_jet_mass = hname;
    name_jet_mass.append("_jet_mass");

    jet_N.push_back(new TH1F(name_jet_N.c_str(),"NJets",100,0,50));
    jet_pt.push_back(new TH1F(name_jet_pt.c_str(),"Jet pT",100,0,2000));
    jet_eta.push_back(new TH1F(name_jet_eta.c_str(),"Jet eta",50,-7,7));
    jet_phi.push_back(new TH1F(name_jet_phi.c_str(),"Jet phi",50,-4,4));
    jet_mass.push_back(new TH1F(name_jet_mass.c_str(),"Jet mass",50,0,600));


    string name_Topjet_N = hname;
    name_Topjet_N.append("_Topjet_N");
    string name_Topjet_pt = hname;
    name_Topjet_pt.append("_Topjet_pT");
    string name_Topjet_eta = hname;
    name_Topjet_eta.append("_Topjet_eta");
    string name_Topjet_phi = hname;
    name_Topjet_phi.append("_Topjet_phi");
    string name_Topjet_mass = hname;
    name_Topjet_mass.append("_Topjet_mass");

    Topjet_N.push_back(new TH1F(name_Topjet_N.c_str(),"NTopjets",100,0,50));
    Topjet_pt.push_back(new TH1F(name_Topjet_pt.c_str(),"Topjet pT",100,0,2000));
    Topjet_eta.push_back(new TH1F(name_Topjet_eta.c_str(),"Topjet eta",50,-7,7));
    Topjet_phi.push_back(new TH1F(name_Topjet_phi.c_str(),"Topjet phi",50,-4,4));
    Topjet_mass.push_back(new TH1F(name_Topjet_mass.c_str(),"Topjet mass",50,0,600));


    string name_Topjet_matched_N = hname;
    name_Topjet_matched_N.append("_Topjet_matched_N");
    string name_Topjet_matched_pt = hname;
    name_Topjet_matched_pt.append("_Topjet_matched_pT");
    string name_Topjet_matched_eta = hname;
    name_Topjet_matched_eta.append("_Topjet_matched_eta");
    string name_Topjet_matched_phi = hname;
    name_Topjet_matched_phi.append("_Topjet_matched_phi");
    string name_Topjet_matched_mass = hname;
    name_Topjet_matched_mass.append("_Topjet_matched_mass");

    Topjet_matched_N.push_back(new TH1F(name_Topjet_matched_N.c_str(),"NTopjets_Matched",100,0,50));
    Topjet_matched_pt.push_back(new TH1F(name_Topjet_matched_pt.c_str(),"Topjet_Matched pT",100,0,2000));
    Topjet_matched_eta.push_back(new TH1F(name_Topjet_matched_eta.c_str(),"Topjet_Matched eta",50,-7,7));
    Topjet_matched_phi.push_back(new TH1F(name_Topjet_matched_phi.c_str(),"Topjet_Matched phi",50,-4,4));
    Topjet_matched_mass.push_back(new TH1F(name_Topjet_matched_mass.c_str(),"Topjet_Matched mass",50,0,600));
    
   
    string name_Top_matched_pt = hname;
    name_Top_matched_pt.append("_Top_matched_pT");
    string name_Top_matched_eta = hname;
    name_Top_matched_eta.append("_Top_matched_eta");
    string name_Top_matched_phi = hname;
    name_Top_matched_phi.append("_Top_matched_phi");
    string name_Top_matched_mass = hname;
    name_Top_matched_mass.append("_Top_matched_mass");

    Top_matched_pt.push_back(new TH1F(name_Top_matched_pt.c_str(),"Top_Matched pT",100,0,2000));
    Top_matched_eta.push_back(new TH1F(name_Top_matched_eta.c_str(),"Top_Matched eta",50,-7,7));
    Top_matched_phi.push_back(new TH1F(name_Top_matched_phi.c_str(),"Top_Matched phi",50,-4,4));
    Top_matched_mass.push_back(new TH1F(name_Top_matched_mass.c_str(),"Top_Matched mass",50,0,600));


    
    string name_Topjet_unmatched_N = hname;
    name_Topjet_unmatched_N.append("_Topjet_unmatched_N");
    string name_Topjet_unmatched_pt = hname;
    name_Topjet_unmatched_pt.append("_Topjet_unmatched_pT");
    string name_Topjet_unmatched_eta = hname;
    name_Topjet_unmatched_eta.append("_Topjet_unmatched_eta");
    string name_Topjet_unmatched_phi = hname;
    name_Topjet_unmatched_phi.append("_Topjet_unmatched_phi");
    string name_Topjet_unmatched_mass = hname;
    name_Topjet_unmatched_mass.append("_Topjet_unmatched_mass");

    Topjet_unmatched_N.push_back(new TH1F(name_Topjet_unmatched_N.c_str(),"NTopjets_Unmatched",100,0,50));
    Topjet_unmatched_pt.push_back(new TH1F(name_Topjet_unmatched_pt.c_str(),"Topjet_Unmatched pT",100,0,2000));
    Topjet_unmatched_eta.push_back(new TH1F(name_Topjet_unmatched_eta.c_str(),"Topjet_Unmatched eta",50,-7,7));
    Topjet_unmatched_phi.push_back(new TH1F(name_Topjet_unmatched_phi.c_str(),"Topjet_Unmatched phi",50,-4,4));
    Topjet_unmatched_mass.push_back(new TH1F(name_Topjet_unmatched_mass.c_str(),"Topjet_Unmatched mass",50,0,600));



    string name_jet_matched_N = hname;
    name_jet_matched_N.append("_jet_N");
    string name_jet_matched_pt = hname;
    name_jet_matched_pt.append("_jet_matched_pT");
    string name_jet_matched_eta = hname;
    name_jet_matched_eta.append("_jet_matched_eta");
    string name_jet_matched_phi = hname;
    name_jet_matched_phi.append("_jet_matched_phi");
    string name_jet_matched_mass = hname;
    name_jet_matched_mass.append("_jet_matched_mass");

    jet_matched_N.push_back(new TH1F(name_jet_matched_N.c_str(),"Njets_Matched",100,0,50));
    jet_matched_pt.push_back(new TH1F(name_jet_matched_pt.c_str(),"Jet_Matched pT",100,0,2000));
    jet_matched_eta.push_back(new TH1F(name_jet_matched_eta.c_str(),"Jet_Matched eta",50,-7,7));
    jet_matched_phi.push_back(new TH1F(name_jet_matched_phi.c_str(),"Jet_Matched phi",50,-4,4));
    jet_matched_mass.push_back(new TH1F(name_jet_matched_mass.c_str(),"Jet_Matched mass",50,0,600));
    
    string name_subjet_N = hname;
    name_subjet_N.append("_subjet_N");
    string name_subjet_pt = hname;
    name_subjet_pt.append("_subjet_pT");
    string name_subjet_eta = hname;
    name_subjet_eta.append("_subjet_eta");
    string name_subjet_phi = hname;
    name_subjet_phi.append("_subjet_phi");
    string name_subjet_mass = hname;
    name_subjet_mass.append("_subjet_mass");

    subjet_N.push_back(new TH1F(name_subjet_N.c_str(),"NSubjets",100,0,50));
    subjet_pt.push_back(new TH1F(name_subjet_pt.c_str(),"Subjet pT",75,0,1500));
    subjet_eta.push_back(new TH1F(name_subjet_eta.c_str(),"Subjet eta",50,-7,7));
    subjet_phi.push_back(new TH1F(name_subjet_phi.c_str(),"Subjet phi",50,-4,4));
    subjet_mass.push_back(new TH1F(name_subjet_mass.c_str(),"Subjet mass",50,0,600));

    string name_Topsubjet_matched_N = hname;
    name_Topsubjet_matched_N.append("_Topsubjet_matched_N");
    string name_Topsubjet_matched_pt = hname;
    name_Topsubjet_matched_pt.append("_Topsubjet_matched_pT");
    string name_Topsubjet_matched_eta = hname;
    name_Topsubjet_matched_eta.append("_Topsubjet_matched_eta");
    string name_Topsubjet_matched_phi = hname;
    name_Topsubjet_matched_phi.append("_Topsubjet_matched_phi");
    string name_Topsubjet_matched_mass = hname;
    name_Topsubjet_matched_mass.append("_Topsubjet_matched_mass");

    Topsubjet_matched_N.push_back(new TH1F(name_Topsubjet_matched_N.c_str(),"NTopsubjet_Matcheds",100,0,50));
    Topsubjet_matched_pt.push_back(new TH1F(name_Topsubjet_matched_pt.c_str(),"Topsubjet_Matched pT",75,0,1500));
    Topsubjet_matched_eta.push_back(new TH1F(name_Topsubjet_matched_eta.c_str(),"Topsubjet_Matched eta",50,-7,7));
    Topsubjet_matched_phi.push_back(new TH1F(name_Topsubjet_matched_phi.c_str(),"Topsubjet_Matched phi",50,-4,4));
    Topsubjet_matched_mass.push_back(new TH1F(name_Topsubjet_matched_mass.c_str(),"Topsubjet_Matched mass",50,0,600));




  }
 
}

// ------------ method called once each job just after ending the event loop  ------------
void SubstrucAnalyzer::endJob() {

  TFile *f = new TFile(OutputName_.c_str(),"recreate");

  for(unsigned int i=0; i<jet_pt.size(); ++i){
    f->mkdir(jetLabels_[i].label().c_str());
    f->cd(jetLabels_[i].label().c_str());
    jet_N[i]->Write();
    jet_pt[i]->Write();
    jet_eta[i]->Write();
    jet_phi[i]->Write();
    jet_mass[i]->Write();
    
    jet_matched_pt[i]->Write();
    jet_matched_eta[i]->Write();
    jet_matched_phi[i]->Write();
    jet_matched_mass[i]->Write();     
    
    subjet_N[i]->Write();
    subjet_pt[i]->Write();
    subjet_eta[i]->Write();
    subjet_phi[i]->Write();
    subjet_mass[i]->Write();
    
    if ( jetLabels_[i].label() == "cmsTopTagPFJets") {
      Topjet_N[i]->Write();
      Topjet_pt[i]->Write();
      Topjet_eta[i]->Write();
      Topjet_phi[i]->Write();
      Topjet_mass[i]->Write();

      Topjet_matched_N[i]->Write();
      Topjet_matched_pt[i]->Write();
      Topjet_matched_eta[i]->Write();
      Topjet_matched_phi[i]->Write();
      Topjet_matched_mass[i]->Write();

      Top_matched_pt[i]->Write();
      Top_matched_eta[i]->Write();
      Top_matched_phi[i]->Write();
      Top_matched_mass[i]->Write();
     
      Topjet_unmatched_N[i]->Write();
      Topjet_unmatched_pt[i]->Write();
      Topjet_unmatched_eta[i]->Write();
      Topjet_unmatched_phi[i]->Write();
      Topjet_unmatched_mass[i]->Write();
      
    }
  }

  delete f;
  
}

// ------------ method called when starting to processes a run  ------------
/*
void 
SubstrucAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
SubstrucAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
SubstrucAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
SubstrucAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SubstrucAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SubstrucAnalyzer);
