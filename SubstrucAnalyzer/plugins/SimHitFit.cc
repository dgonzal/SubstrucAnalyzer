// -*- C++ -*-
//
// Package:    Substructure/SubstrucAnalyzer
// Class:      SubstrucAnalyzer
// 
/**\class SimHitFit SimHitFit.cc Substructure/SimHitFit/plugins/SimHitFit.cc

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


#include "Substructure/SubstrucAnalyzer/plugins/SimHitFit.h"


using namespace edm;
using namespace std;
using namespace reco; 


struct sub_det{
  unsigned int detId;
  double energy;
  double time_energy;
  double entries;
  double EMenergy;
};



SimHitFit::SimHitFit(const edm::ParameterSet& ps)
{

  Module_ = ps.getParameter<string>("ProducerModule");
  OutputName_ = ps.getParameter<string>("OutputName");
  
  ResponseBinning_ = ps.getParameter<vector<double>>("ResponseBinning");
  

  InputTag ECALEBTag(Module_,"EcalHitsEB");
  InputTag ECALEETag(Module_,"EcalHitsEE");
  InputTag ECALESTag(Module_,"EcalHitsES");
  InputTag HCALTag  (Module_,"HcalHits"  );

  //consumes<reco::GenParticleCollection>(ps.getParameter<edm::InputTag>("GenParticle"));

  ECALEBToken = consumes<edm::PCaloHitContainer>(ECALEBTag);
  ECALEEToken = consumes<edm::PCaloHitContainer>(ECALEETag);
  ECALESToken = consumes<edm::PCaloHitContainer>(ECALESTag);
  HCALToken   = consumes<edm::PCaloHitContainer>(HCALTag);   

}


SimHitFit::~SimHitFit()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SimHitFit::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  vector<double> samplingFactors = {125.44, 125.54, 125.32, 125.13, 
				    124.46, 125.01, 125.22, 125.48, 
				    124.45, 125.9, 125.83, 127.01, 
				    126.82, 129.73, 131.83, 143.52,
				    210.55, 197.93, 186.12, 189.64, 
				    189.63, 190.28, 189.61, 189.6, 
				    190.12, 191.22, 190.9, 193.06, 
				    188.42, 188.42, 0.383, 0.368,
				    231.0, 231.0, 231.0, 231.0, 360.0,
				    360.0, 360.0, 360.0, 360.0, 360.0,
				    360.0, 360.0, 360.0, 360.0, 360.0};

  
  edm::Handle<reco::GenParticleCollection> GenParticles;
  iEvent.getByLabel("genParticles","",GenParticles);
  
  for(unsigned int i = 0; i<GenParticles->size();i++){
    reco::GenParticle gen =  GenParticles->at(i);
    gen_energy->Fill(gen.energy());
    gen_pT->Fill(gen.pt());
    gen_eta->Fill(gen.eta());
    gen_phi->Fill(gen.phi());
    gen_pdgId->Fill(gen.pdgId());
  }

  edm::Handle<edm::PCaloHitContainer> ECALEBSimHits;
  edm::Handle<edm::PCaloHitContainer> ECALEESimHits;
  edm::Handle<edm::PCaloHitContainer> ECALESSimHits;
  edm::Handle<edm::PCaloHitContainer> HCALSimHits;
  

  iEvent.getByToken(ECALEBToken,ECALEBSimHits);
  iEvent.getByToken(ECALEEToken,ECALEESimHits);
  iEvent.getByToken(ECALESToken,ECALESSimHits);
  iEvent.getByToken(HCALToken  ,HCALSimHits);
 
  vector<edm::Handle<edm::PCaloHitContainer>> CalorimeterSimHits;
  CalorimeterSimHits.push_back(ECALEBSimHits);
  CalorimeterSimHits.push_back(ECALEESimHits);
  CalorimeterSimHits.push_back(ECALESSimHits);
  CalorimeterSimHits.push_back(HCALSimHits);

  edm::ESHandle<CaloGeometry> pG;   // you might actuallly need a different type here (instead of CaloGeometry), but this will already give very reasonable results. fine tuning can happen later.
  iSetup.get<CaloGeometryRecord>().get(pG);
  const CaloGeometry* geo = pG.product();

  double hcal_energy  = 0;
  double ecal_energy = 0;
  double EMenergy = 0;
  double Hadenergy = 0;

  PCaloHit Hit;
  GlobalPoint pos;

  
  vector<sub_det> hcal_energy_depo;
  vector<sub_det> ecal_energy_depo;


  for(unsigned int p=0;p<CalorimeterSimHits.size();++p){
    edm::Handle<edm::PCaloHitContainer> SimHits = CalorimeterSimHits[p];

    for(unsigned int m =0; m<SimHits->size(); ++m){
      Hit = SimHits->at(m);
      pos = geo->getPosition(Hit.id());

      CaloDepth->Fill(Hit.depth());
      //if(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi())>0.5) continue;
      
      

      if(std::isnan(Hit.energy()) || std::isinf(Hit.energy())){
	if(std::isnan(Hit.energy())) cout<<" Hit is Nan"<<endl;
	if(std::isinf(Hit.energy())) cout<<" Hit is Inf"<<endl;
	//assert(1==0);
	//continue;
      }

      //cout<<"type: "<< p <<" Hit Energy: "<< Hit.energy()<<" Hit depth: "<<Hit.depth()<< " Hit time: "<<Hit.time()<<endl;


      double sampling_weight = 1;

      if(p==3){
	HcalDetId hcalId(Hit.id());
	if(hcalId.ietaAbs()<0) {
	  cout<<"negativ ietaAbs"<<endl;
	  //assert(0==1);
	}
	if(samplingFactors.size() < abs(hcalId.ietaAbs())){
	  cout<<"sampling Factor not covering the ietaAbs"<<endl;
	  //assert(0==1);
	}
	sampling_weight = samplingFactors[hcalId.ietaAbs()-1];
	hcal_energy_response->Fill(Hit.energy()*sampling_weight);
	hcalDepth->Fill(Hit.depth());
	HitEnergyEM->Fill(Hit.energyEM()*sampling_weight);

	hcal_energy += Hit.energy()*sampling_weight;
	EMenergy+= Hit.energyEM()*sampling_weight;
	Hadenergy+= Hit.energyHad()*sampling_weight;

	bool detId_flag =true;
	for(unsigned int m=0; m<hcal_energy_depo.size();++m){
	  if(hcal_energy_depo.at(m).detId == Hit.id()){
	    hcal_energy_depo.at(m).energy += Hit.energy();
	    hcal_energy_depo.at(m).time_energy += Hit.time()*Hit.energy();
	    hcal_energy_depo.at(m).entries +=1;
	    hcal_energy_depo.at(m).EMenergy += Hit.energyEM();
	    detId_flag = false;
	  }
	}
	if(detId_flag){
	  sub_det det_part;
	  det_part.detId = Hit.id();
	  det_part.energy = Hit.energy();
	  det_part.time_energy = Hit.time()*Hit.energy();
	  det_part.entries = 1;
	  det_part.EMenergy = Hit.energyEM();
	  hcal_energy_depo.push_back(det_part);
	}
      } 
      else{
	ecal_energy += Hit.energy()*sampling_weight;
	EMenergy+= Hit.energyEM()*sampling_weight;
	Hadenergy+= Hit.energyHad()*sampling_weight;
	ecal_energy_response->Fill(Hit.energy()*sampling_weight);
	ecalDepth->Fill(Hit.depth());
	HitEnergyEM->Fill(Hit.energyEM()*sampling_weight);
	
	bool detId_flag =true;
	for(unsigned int m=0; m<ecal_energy_depo.size();++m){
	  if(ecal_energy_depo.at(m).detId == Hit.id()){
	    ecal_energy_depo.at(m).energy += Hit.energy();
	    ecal_energy_depo.at(m).time_energy += Hit.time()*Hit.energy();
	    ecal_energy_depo.at(m).entries +=1;
	    ecal_energy_depo.at(m).EMenergy += Hit.energyEM();
	    detId_flag = false;
	  }
	}
	if(detId_flag){
	  sub_det det_part;
	  det_part.detId = Hit.id();
	  det_part.energy = Hit.energy();
	  det_part.time_energy = Hit.time()*Hit.energy();
	  det_part.entries = 1;
	  det_part.EMenergy = Hit.energyEM();
	  ecal_energy_depo.push_back(det_part);
	}
      }
    
      eta_response->Fill(pos.eta());
      phi_response->Fill(pos.phi());
      simhit_energy_response->Fill(Hit.energy()*sampling_weight);
      scatter_x_y_z->Fill(pos.x(),pos.y(),pos.z());
    }
  }
  
  for(unsigned int m=0; m<ecal_energy_depo.size();++m){
    ecal_detId_energy_response->Fill(ecal_energy_depo.at(m).energy);
    ecal_detId_genenergy->Fill(GenParticles->at(0).energy(),ecal_energy_depo.at(m).energy);
   
    if(ecal_energy_depo.at(m).energy>0.){
      //cout<<"ecal time energy: "<<ecal_energy_depo.at(m).time_energy<<" energy: "<<ecal_energy_depo.at(m).energy<<" average time: "<<ecal_energy_depo.at(m).time_energy/ecal_energy_depo.at(m).energy<<endl;
      ecal_detId_EnergyEM->Fill(ecal_energy_depo.at(m).EMenergy);

      ecalTime->Fill(ecal_energy_depo.at(m).time_energy/ecal_energy_depo.at(m).energy/ecal_energy_depo.at(m).entries);
      CaloTime->Fill(ecal_energy_depo.at(m).time_energy/ecal_energy_depo.at(m).energy/ecal_energy_depo.at(m).entries);

      ecalTime_Eweight->Fill(ecal_energy_depo.at(m).time_energy/ecal_energy_depo.at(m).energy,ecal_energy_depo.at(m).energy);
      CaloTime_Eweight->Fill(ecal_energy_depo.at(m).time_energy/ecal_energy_depo.at(m).energy,ecal_energy_depo.at(m).energy);

    }
  }

  //cout<<hcal_energy_depo.size()<<endl;

 
  for(unsigned int m=0; m<hcal_energy_depo.size();++m){
    hcal_detId_energy_response->Fill(hcal_energy_depo.at(m).energy);
    hcal_detId_genenergy->Fill(GenParticles->at(0).energy(),hcal_energy_depo.at(m).energy);
    
    if(hcal_energy_depo.at(m).energy>0.){
      //cout<<" entries "<<hcal_energy_depo.at(m).entries  <<" hcal time energy: "<<hcal_energy_depo.at(m).time_energy<<" energy: "<<hcal_energy_depo.at(m).energy<<" average time: "<<hcal_energy_depo.at(m).time_energy/hcal_energy_depo.at(m).energy<<endl;

      if(hcal_energy_depo.at(m).energy!=hcal_energy_depo.at(m).energy) assert(0==1);
      hcal_detId_EnergyEM->Fill(hcal_energy_depo.at(m).EMenergy);

      hcalTime->Fill(hcal_energy_depo.at(m).time_energy/hcal_energy_depo.at(m).energy);
      CaloTime->Fill(hcal_energy_depo.at(m).time_energy/hcal_energy_depo.at(m).energy);

      hcalTime_Eweight->Fill(hcal_energy_depo.at(m).time_energy/hcal_energy_depo.at(m).energy,hcal_energy_depo.at(m).energy);
      CaloTime_Eweight->Fill(hcal_energy_depo.at(m).time_energy/hcal_energy_depo.at(m).energy,hcal_energy_depo.at(m).energy);



    }
  }

  relativCalenergy->Fill(ecal_energy/(ecal_energy+hcal_energy));	
  ratioEcalHcal_energy->Fill(ecal_energy/hcal_energy);
  ratioEMTotal_energy->Fill(EMenergy/(ecal_energy+hcal_energy));
  genenergy_relativCalenergy->Fill(GenParticles->at(0).energy(),(ecal_energy)/(ecal_energy+hcal_energy));
  
  energy_response->Fill(hcal_energy+ecal_energy);

  TotalEnergyEM->Fill(EMenergy);
  TotalEnergyHad->Fill(Hadenergy);

  eta_gen=GenParticles->at(0).eta();  phi_gen = GenParticles->at(0).phi(); e_ecal=ecal_energy ; e_hcal=hcal_energy ; e_calo=hcal_energy +ecal_energy ;
  eEM_calo = EMenergy; eHad_calo = Hadenergy;

  if(e_calo>200){
    cout<<"too high energy"<<endl;
    //assert(0==1);
  }
  tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void SimHitFit::beginJob(){

  gen_energy = new TH1F("gen_energy","Energy of Genparticles",3001,-0.5,3000.5);
  gen_pT = new TH1F("gen_pT","pT of Genparticles",3001,-0.5,3000.5);
  gen_eta = new TH1F("gen_eta","#Eta of Genparticles",100,0,3);
  gen_phi = new TH1F("gen_phi","#phi of Genparticles",100,-3.14159265359,3.14159265359);
  gen_pdgId = new TH1F("gen_pdgId","PdgId of Genparticles",1001,-0.5,1000);

  energy_response = new TH1F("energy_response","Energy Response",ResponseBinning_[0],ResponseBinning_[1],ResponseBinning_[2]);

  relativCalenergy = new TH1F("relativCalenergy","relativ Energy in ECAL",300,0,1);
  genenergy_relativCalenergy = new TH2F("genenergy_relativCalenergy","",150,0,20,300,0,1);

  eta_response = new TH1F("eta_response","#eta",50,-4,4);
  phi_response = new TH1F("phi_response","#phi",50,-4,4);
  simhit_energy_response = new TH1F("simhit_energy_response","E",150,0,70);

  ecal_energy_response = new TH1F("eb_energy_response","E",150,0,40);
  hcal_energy_response = new TH1F("hcal_energy_response","E",150,0,50);
 
  hcal_detId_energy_response = new TH1F("hcal_detId_energy_response","E",250,0,10);
  ecal_detId_energy_response = new TH1F("ecal_detId_energy_response","E",250,0,10);

  hcal_detId_genenergy = new TH2F("hcal_detId_genenergy","E_{DetId}, E_{gen}",150,0,20,250,0,10);
  ecal_detId_genenergy = new TH2F("ecal_detId_genenergy","E_{DetId}, E_{gen}",150,0,20,250,0,10);

 
  scatter_x_y_z = new TH3F("scatter_x_y_z","",150,-500,500,150,-500,500,150,-800,800);

  //----
  HitEnergyEM = new TH1F("hit_energy_em","EM Energy of Hits",150,0,8);
  TotalEnergyEM = new TH1F("energy_em","EM Energy",ResponseBinning_[0],ResponseBinning_[1],ResponseBinning_[2]);
  TotalEnergyHad = new TH1F("energy_had","Hadronic Energy",ResponseBinning_[0],ResponseBinning_[1],ResponseBinning_[2]);
  hcal_detId_EnergyEM = new TH1F("hcal_detId_EnergyEM","hcal EnergyEM per detId",150,0,6);
  ecal_detId_EnergyEM = new TH1F("ecal_detId_EnergyEM","ecal EnergyEM per detId",150,0,6);
 
  CaloTime = new TH1F("calo_time","Time in the Calorimeter",250,0,150);
  ecalTime = new TH1F("ecal_time","Time in ECAL",250,0,150);
  hcalTime = new TH1F("hcal_time","Time in HCAL",250,0,150);

  CaloTime_Eweight = new TH1F("calo_time_eweight","Time Energy weight in the Calorimeter",250,0,150);
  ecalTime_Eweight = new TH1F("ecal_time_eweight","Time Energy weight in ECAL",250,0,150);
  hcalTime_Eweight = new TH1F("hcal_time_eweight","Time Energy weight in HCAL",250,0,150);

  CaloDepth = new TH1F("calo_depth","Calorimeter Depth",100,-.05,99.5);
  ecalDepth = new TH1F("ecal_depth","ECAL Depth",100,-.05,99.5);
  hcalDepth = new TH1F("hcal_depth","HCAL Depth",100,-.05,99.5);
  
  ratioEcalHcal_energy = new TH1F("ratioEcalHcal_energy","ratio Ecal/Hcal Energy",100,0,10);
  ratioEMTotal_energy = new TH1F("ratioEMTotal_energy","ratio EM Energy",100,0,2);
  //----

  f = new TFile(OutputName_.c_str(),"RECREATE");
  tree = new TTree("Total", "Energy Calorimeter info");
  tree->Branch("eta",&eta_gen,"eta_gen/D");
  tree->Branch("phi",&phi_gen,"phi_gen/D");
  tree->Branch("e_ecal",&e_ecal,"e_ecal/D");
  tree->Branch("e_hcal",&e_hcal,"e_hcal/D");
  tree->Branch("e_calo",&e_calo,"e_calo/D");
  tree->Branch("eEM_calo",&eEM_calo,"eEM_calo/D");
  tree->Branch("eHad_calo",&eHad_calo,"eHad_calo/D");
}

// ------------ method called once each job just after ending the event loop  ------------
void SimHitFit::endJob() {

  f->cd();

  tree->Write();


  gen_energy->Write();
  gen_pT->Write();
  gen_eta->Write();
  gen_phi->Write();
  gen_pdgId->Write();

  eta_response->Write();
  phi_response->Write();
  simhit_energy_response->Write();

  ecal_energy_response->Write();
  hcal_energy_response->Write();
 
  hcal_detId_energy_response->Write();
  ecal_detId_energy_response->Write();
  ecal_detId_genenergy->Write();
  hcal_detId_genenergy->Write();

  scatter_x_y_z->Write();

  relativCalenergy->Write(); 
  genenergy_relativCalenergy->Write();

  energy_response->Write();

  //----
  HitEnergyEM->Write(); 
  TotalEnergyEM->Write(); 
  TotalEnergyHad->Write(); 
  hcal_detId_EnergyEM->Write();
  ecal_detId_EnergyEM->Write(); 
 
  CaloTime->Write();
  ecalTime->Write(); 
  hcalTime->Write();

  CaloTime_Eweight->Write();
  ecalTime_Eweight->Write(); 
  hcalTime_Eweight->Write();

  CaloDepth->Write(); 
  ecalDepth->Write();
  hcalDepth->Write(); 
  
  ratioEcalHcal_energy->Write(); 
  ratioEMTotal_energy->Write();

  //---
  delete f;

}

// ------------ method called when starting to processes a run  ------------
/*
void 
SimHitFit::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
SimHitFit::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
SimHitFit::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
SimHitFit::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SimHitFit::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimHitFit);
