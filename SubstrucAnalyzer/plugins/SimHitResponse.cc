
// -*- C++ -*-
//
// Package:    Substructure/SubstrucAnalyzer
// Class:      SubstrucAnalyzer
// 
/**\class SimHitResponse SimHitResponse.cc Substructure/SimHitResponse/plugins/SimHitResponse.cc

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


#include "Substructure/SubstrucAnalyzer/plugins/SimHitResponse.h"


using namespace edm;
using namespace std;
using namespace reco; 


struct sub_det{
  unsigned int detId;
  double energy;
};



SimHitResponse::SimHitResponse(const edm::ParameterSet& ps)
{
  //jet_match_ =ps.getParameter<edm::InputTag>("jets_to_match");

  //vector<PCaloHit> "g4SimHits" "HcalHits" "SIM"     

  //SimHits_ = consumes< edm::View<std::vector<PCaloHit> > ("HcalHits");

  Module_ = ps.getParameter<string>("ProducerModule");
  OutputName_ = ps.getParameter<string>("OutputName");
  //Hits_Label_ = consumes<edm::PCaloHitContainer>(ps.getParameter<edm::InputTag>("HitLabel"));
  //ECAL_Hits_Label_ = consumes<edm::PCaloHitContainer>(ps.getParameter<edm::InputTag>("ECALHitLabel"));

}


SimHitResponse::~SimHitResponse()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SimHitResponse::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edm::PCaloHitContainer> ECALEBSimHits;
  edm::Handle<edm::PCaloHitContainer> ECALEESimHits;
  edm::Handle<edm::PCaloHitContainer> ECALESSimHits;
  edm::Handle<edm::PCaloHitContainer> HCALSimHits;

  edm::Handle<reco::GenParticleCollection> GenParticles;

  iEvent.getByLabel(Module_,"HcalHits",HCALSimHits); 
  iEvent.getByLabel(Module_,"EcalHitsEB",ECALEBSimHits); 
  iEvent.getByLabel(Module_,"EcalHitsEE",ECALEESimHits); 
  iEvent.getByLabel(Module_,"EcalHitsES",ECALESSimHits); 

  iEvent.getByLabel("genParticles","",GenParticles);

  
  for(unsigned int i = 0; i<GenParticles->size();i++){
    reco::GenParticle gen =  GenParticles->at(i);
    std::cout<<" particleID "<<gen.pdgId() <<" gen pT "<< gen.pt()<<" phi "<<gen.phi()<< " eta "<<gen.eta()<< " e "<< gen.energy()<<std::endl;
    }
  


  gen_energy->Fill(GenParticles->at(0).energy());
  gen_pT->Fill(GenParticles->at(0).pt());

  cout<<"NSimHits ECALEB "<<ECALEBSimHits->size()<< "  ECALES "<<ECALESSimHits->size()<< "  ECALEE "<<ECALEESimHits->size()<< " HCAL " <<HCALSimHits->size() <<endl;

  //if(!found_simhit && ) return;

  edm::ESHandle<CaloGeometry> pG;   // you might actuallly need a different type here (instead of CaloGeometry), but this will already give very reasonable results. fine tuning can happen later.
  iSetup.get<CaloGeometryRecord>().get(pG);
  const CaloGeometry* geo = pG.product();

  double hcal_energy  = 0;
  double hcal_energy_weight  = 0;
  double ecaleb_energy =0;
  double ecalee_energy =0;
  double ecales_energy =0;

  PCaloHit Hit;
  GlobalPoint pos;

  double depth=0;

  double max_posX=0;
  double max_posY=0;


  vector<sub_det> ecal_energy_depo;


  for(unsigned int i = 0; i<ECALEBSimHits->size(); ++i){
    Hit = ECALEBSimHits->at(i);
    pos = geo->getPosition(Hit.id());
    depth = sqrt(pos.x()*pos.x()+pos.y()*pos.y());
    
    bool detId_flag =true;
    for(unsigned int m=0; m<ecal_energy_depo.size();++m){
      if(ecal_energy_depo.at(m).detId == Hit.id()){
	ecal_energy_depo.at(m).energy+= Hit.energy();
	detId_flag = false;
      }
    }
    if(detId_flag){
      sub_det det_part;
      det_part.detId = Hit.id();
      det_part.energy = Hit.energy();
      ecal_energy_depo.push_back(det_part);
    }

    if(max_posX<pos.x()) max_posX=pos.x();
    if(max_posY<pos.y()) max_posY=pos.y();

    scatter_dR_depth->Fill(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()),depth);
    scatter_dR_pT->Fill(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()),GenParticles->at(0).pt());
    scatter_dR_pT_e->Fill(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()),GenParticles->at(0).pt(),Hit.energy());

    eta_response->Fill( pos.eta());
    phi_response->Fill( pos.phi());

    energy_response->Fill(Hit.energy());
    eb_energy_response->Fill(Hit.energy());

    ecaleb_energy+= Hit.energy();

    scatter_genpT_hitenergy->Fill(GenParticles->at(0).pt(),Hit.energy() );
    scatter_genenergy_hitenergy->Fill(GenParticles->at(0).energy(), Hit.energy());

  }

 
  for(unsigned int i = 0; i<ECALEESimHits->size(); ++i){
    Hit = ECALEESimHits->at(i);
    pos = geo->getPosition(Hit.id());
    depth = sqrt(pos.x()*pos.x()+pos.y()*pos.y());
    
    bool detId_flag =true;
    for(unsigned int m=0; m<ecal_energy_depo.size();++m){
      if(ecal_energy_depo.at(m).detId == Hit.id()){
	ecal_energy_depo.at(m).energy+=Hit.energy();
	detId_flag = false;
      }
    }
    if(detId_flag){
      sub_det det_part;
      det_part.detId = Hit.id();
      det_part.energy = Hit.energy();
      ecal_energy_depo.push_back(det_part);
    }

    if(max_posX<pos.x()) max_posX=pos.x();
    if(max_posY<pos.y()) max_posY=pos.y();

    scatter_dR_depth->Fill(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()),depth);
    scatter_dR_pT->Fill(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()),GenParticles->at(0).pt());
    scatter_dR_pT_e->Fill(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()),GenParticles->at(0).pt(),Hit.energy());

    eta_response->Fill( pos.eta());
    phi_response->Fill( pos.phi());

    ecalee_energy+= Hit.energy();

    energy_response->Fill(Hit.energy());
    ee_energy_response->Fill(Hit.energy());

    scatter_genpT_hitenergy->Fill(GenParticles->at(0).pt(),Hit.energy());
    scatter_genenergy_hitenergy->Fill(GenParticles->at(0).energy(), Hit.energy());


  }

  for(unsigned int i = 0; i<ECALESSimHits->size(); ++i){
    Hit = ECALESSimHits->at(i);
    pos = geo->getPosition(Hit.id());
    depth = sqrt(pos.x()*pos.x()+pos.y()*pos.y());
        
    bool detId_flag =true;
    for(unsigned int m=0; m<ecal_energy_depo.size();++m){
      if(ecal_energy_depo.at(m).detId == Hit.id()){
	ecal_energy_depo.at(m).energy+=Hit.energy();
	detId_flag = false;
      }
    }
    if(detId_flag){
      sub_det det_part;
      det_part.detId = Hit.id();
      det_part.energy = Hit.energy();
      ecal_energy_depo.push_back(det_part);
    }

    if(max_posX<pos.x()) max_posX=pos.x();
    if(max_posY<pos.y()) max_posY=pos.y();

    scatter_dR_depth->Fill(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()),depth);
    scatter_dR_pT->Fill(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()),GenParticles->at(0).pt());
    scatter_dR_pT_e->Fill(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()),GenParticles->at(0).pt(),Hit.energy());

    eta_response->Fill( pos.eta());
    phi_response->Fill( pos.phi());

    ecales_energy+= Hit.energy();

    energy_response->Fill(Hit.energy());
    es_energy_response->Fill(Hit.energy());



    scatter_genpT_hitenergy->Fill(GenParticles->at(0).pt(),Hit.energy());
    scatter_genenergy_hitenergy->Fill(GenParticles->at(0).energy(), Hit.energy());

  }


  vector<sub_det> det_energy_depo;
  double samplingFactors[] = {125.44, 125.54, 125.32,125.13, 124.46,

                    125.01, 125.22, 125.48, 124.45, 125.9,

                    125.83, 127.01, 126.82, 129.73, 131.83,

                    143.52, 210.55, 197.93, 186.12,189.64, 189.63,

                    190.28, 189.61, 189.6, 190.12, 191.22,

                    190.9, 193.06, 188.42, 188.42};



  for (unsigned int i = 0; i< HCALSimHits->size(); ++i){
    bool detId_flag =true;

    Hit = HCALSimHits->at(i);
    HcalDetId hcalId(Hit.id());
    //cout<<hcalId.ietaAbs()<<endl;
    
    double weighted_energy;


    
    weighted_energy= Hit.energy()*samplingFactors[hcalId.ietaAbs()-1];
    //if(hcalId.ietaAbs() >  16 ) weighted_energy= Hit.energy()*180;


    for(unsigned int m=0; m<det_energy_depo.size();++m){
      if(det_energy_depo.at(m).detId == Hit.id()){
	det_energy_depo.at(m).energy+= weighted_energy;
	detId_flag = false;
      }
    }
    if(detId_flag){
      sub_det det_part;
      det_part.detId = Hit.id();
      det_part.energy = weighted_energy;
      det_energy_depo.push_back(det_part);
    }


    hcal_energy+=Hit.energy();
    hcal_energy_weight+=weighted_energy;

    pos = geo->getPosition(Hit.id());
    depth = sqrt(pos.x()*pos.x()+pos.y()*pos.y());
    
    if(max_posX<pos.x()) max_posX=pos.x();
    if(max_posY<pos.y()) max_posY=pos.y();

    scatter_dR_depth->Fill(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()),depth);
    scatter_dR_pT->Fill(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()),GenParticles->at(0).pt());
    scatter_dR_pT_e->Fill(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()),GenParticles->at(0).pt(),Hit.energy());

    scatter_genpT_hitenergy->Fill(GenParticles->at(0).pt(), weighted_energy);
    scatter_genenergy_hitenergy->Fill(GenParticles->at(0).energy(), weighted_energy);

    eta_response->Fill( pos.eta());
    phi_response->Fill( pos.phi());
    energy_response->Fill(weighted_energy);
    hcal_energy_response->Fill(Hit.energy());

    hcal_weight_energy_response->Fill(weighted_energy);

    //cout<<"energy: "<<hit_energy<<" phi: "<<hit_phi<<" eta: "<<hit_eta<<endl;

    //cout<<pos.x()<<" "<<pos.y()<<" "<<pos.z()<<endl;
    
    //if(HCALSimHits->size()>1000){

    scatter_x_z->Fill(pos.x(),pos.z());
    scatter_y_z->Fill(pos.y(),pos.z());
    scatter_x_y_z->Fill(pos.x(),pos.y(),pos.z());
      //}

  }

  for(unsigned int m=0; m<ecal_energy_depo.size();++m){
    ecal_detId_energy_response->Fill(ecal_energy_depo.at(m).energy);
    ecal_detId_genenergy->Fill(GenParticles->at(0).energy(),ecal_energy_depo.at(m).energy);
  }

  for(unsigned int m=0; m<det_energy_depo.size();++m){
    hcal_detId_energy_response->Fill(det_energy_depo.at(m).energy);
    hcal_detId_genenergy->Fill(GenParticles->at(0).energy(),det_energy_depo.at(m).energy);
  }

  eb_eta_hits->Fill(GenParticles->at(0).eta(),ECALEBSimHits->size(),ecaleb_energy);
  ee_eta_hits->Fill(GenParticles->at(0).eta(),ECALEESimHits->size(),ecales_energy);
  es_eta_hits->Fill(GenParticles->at(0).eta(),ECALESSimHits->size(),ecalee_energy);
  hcal_eta_hits->Fill(GenParticles->at(0).eta(),HCALSimHits->size(),hcal_energy_weight);

  eb_pt_hits->Fill(GenParticles->at(0).pt(),ECALEBSimHits->size(),ecaleb_energy);
  ee_pt_hits->Fill(GenParticles->at(0).pt(),ECALEESimHits->size(),ecales_energy);
  es_pt_hits->Fill(GenParticles->at(0).pt(),ECALESSimHits->size(),ecalee_energy);
  hcal_pt_hits->Fill(GenParticles->at(0).pt(),HCALSimHits->size(),hcal_energy_weight);

  eb_pt_hits->Fill(GenParticles->at(0).pt(),ECALEBSimHits->size(),ecaleb_energy);
  ee_pt_hits->Fill(GenParticles->at(0).pt(),ECALEESimHits->size(),ecales_energy);
  es_pt_hits->Fill(GenParticles->at(0).pt(),ECALESSimHits->size(),ecalee_energy);
  hcal_pt_hits->Fill(GenParticles->at(0).pt(),HCALSimHits->size(),hcal_energy_weight);
  

  scatter_genpT_energy->Fill(GenParticles->at(0).pt(),ecaleb_energy+ecalee_energy+ecales_energy+hcal_energy_weight);    
  scatter_genenergy_energy->Fill(GenParticles->at(0).energy(),(ecaleb_energy+ecalee_energy+ecales_energy+hcal_energy_weight-GenParticles->at(0).energy())/GenParticles->at(0).energy());    
  
  relativCalenergy->Fill((ecaleb_energy+ecalee_energy+ecales_energy)/(ecaleb_energy+ecalee_energy+ecales_energy+hcal_energy_weight));	
  genenergy_relativCalenergy->Fill(GenParticles->at(0).energy(),(ecaleb_energy+ecalee_energy+ecales_energy)/(ecaleb_energy+ecalee_energy+ecales_energy+hcal_energy_weight));

  scatter_eta_hits->Fill(GenParticles->at(0).eta(),ECALEBSimHits->size()+ECALEESimHits->size()+ECALESSimHits->size()+HCALSimHits->size(),hcal_energy_weight+ecaleb_energy+ecales_energy+ecalee_energy);
  scatter_pT_hits->Fill(GenParticles->at(0).pt(),ECALEBSimHits->size()+ECALEESimHits->size()+ECALESSimHits->size()+HCALSimHits->size(),hcal_energy_weight+ecaleb_energy+ecales_energy+ecalee_energy);
  scatter_pT_hitenergy->Fill(GenParticles->at(0).pt(),hcal_energy_weight+ecaleb_energy+ecales_energy+ecalee_energy);
  

  std::cout<<"hcal energy "<<hcal_energy<<" corrected "<<hcal_energy_weight<< " ecal EB energy "<< ecaleb_energy<< " ecal ES energy "<< ecales_energy << " ecal EE energy "<< ecalee_energy << std::endl;

  std::cout<<" max (x,y) "<<max_posX<<" "<<max_posY<<std::endl;

  if(GenParticles->at(0).energy()>8 &&GenParticles->at(0).energy()<10)  barrel_9gev_energy->Fill(hcal_energy_weight+ecaleb_energy+ecales_energy+ecalee_energy);

  //henergy_sum->Fill(hcal_energy);

  if(ECALEBSimHits->size()>0)eb_energy->Fill(ecaleb_energy);
  if(ECALEESimHits->size()>0)ee_energy->Fill(ecalee_energy);
  if(ECALESSimHits->size()>0)es_energy->Fill(ecales_energy);
  if(HCALSimHits->size()>0)h_energy->Fill(hcal_energy);


}


// ------------ method called once each job just before starting event loop  ------------
void SimHitResponse::beginJob(){


  barrel_9gev_energy = new TH1F("barrel_9gev_energy","Energy 9GeV K+",50,0,20);

  gen_energy = new TH1F("gen_energy","Energy of Genparticles",50,0,25);
  gen_pT = new TH1F("gen_pT","pT of Genparticles",50,0,20);

  relativCalenergy = new TH1F("relativCalenergy","relativ Energy in ECAL",300,0,1);
  genenergy_relativCalenergy = new TH2F("genenergy_relativCalenergy","",150,0,20,300,0,1);

  eta_response = new TH1F("eta_response","#eta",50,-4,4);
  phi_response = new TH1F("phi_response","#phi",50,-4,4);
  energy_response = new TH1F("energy_response","E",150,0,0.002);

  eb_energy_response = new TH1F("eb_energy_response","E",150,0,2);//0.002
  ee_energy_response = new TH1F("ee_energy_response","E",150,0,2);
  es_energy_response = new TH1F("es_energy_response","E",150,0,2);
  hcal_energy_response = new TH1F("hcal_energy_response","E",150,0,2);
  hcal_weight_energy_response = new TH1F("hcal_weight_energy_response","E",150,0,2);
  hcal_detId_energy_response = new TH1F("hcal_detId_energy_response","E",250,0,10);
  ecal_detId_energy_response = new TH1F("ecal_detId_energy_response","E",250,0,10);

  hcal_detId_genenergy = new TH2F("hcal_detId_genenergy","E_{DetId}, E_{gen}",150,0,20,250,0,10);
  ecal_detId_genenergy = new TH2F("ecal_detId_genenergy","E_{DetId}, E_{gen}",150,0,20,250,0,10);

  eb_energy = new TH1F("eb_energy","E",100,0,40);
  ee_energy = new TH1F("ee_energy","E",100,0,1);
  es_energy = new TH1F("es_energy","E",100,0,0.1);
  h_energy= new TH1F("hcal_energy","E",100,0,0.5);

  scatter_dR_depth = new TH2F("scatter_dR_depth","",50,0,6,250,0,500);
  scatter_dR_pT = new TH2F("scatter_dR_pT","",50,0,6,50,0,20);
  scatter_dR_pT_e = new TH2F("scatter_dR_pT_e","",50,0,6,50,0,20);

  scatter_pT_hits = new TH2F("scatter_pT_hits","",50,0,10,100,0,650);
  scatter_eta_hits = new TH2F("scatter_eta_hits","",50,-2.5,2.5,100,0,650);
  scatter_pT_hitenergy = new TH2F("scatter_pT_hitenergy","",150,0,10,150,0,10);

  eb_eta_hits = new TH2F("eb_eta_hits","",50,-2.5,2.5,100,0,650);
  ee_eta_hits = new TH2F("ee_eta_hits","",50,-2.5,2.5,100,0,650);
  es_eta_hits = new TH2F("es_eta_hits","",50,-2.5,2.5,100,0,650);
  hcal_eta_hits = new TH2F("hcal_eta_hits","",50,-2.5,2.5,100,0,650);

  eb_pt_hits = new TH2F("eb_pt_hits","",150,0,10,100,0,650);
  ee_pt_hits = new TH2F("ee_pt_hits","",150,0,10,100,0,650);
  es_pt_hits = new TH2F("es_pt_hits","",150,0,10,100,0,650);
  hcal_pt_hits = new TH2F("hcal_pt_hits","",150,0,10,100,0,650);

  scatter_genpT_hitenergy     = new TH2F("genpt_hitenergy","",150,0,20,150,0,21);
  scatter_genenergy_hitenergy = new TH2F("genenergy_hitenergy","",150,0,21,150,0,21);

  scatter_genpT_energy     = new TH2F("genpt_energy","",150,0,21,150,0,45);
  scatter_genenergy_energy = new TH2F("genenergy_energy","",150,0,21,150,-1.1,5);


  scatter_x_z = new TH2F("scatter_x_z","",150,-500,500,150,-800,800);
  scatter_y_z = new TH2F("scatter_y_z","",150,-500,500,150,-800,800);
  scatter_x_y_z = new TH3F("scatter_x_y_z","",150,-500,500,150,-500,500,150,-800,800);
}

// ------------ method called once each job just after ending the event loop  ------------
void SimHitResponse::endJob() {

  TFile *f = new TFile(OutputName_.c_str(),"recreate");
  gen_energy->Write();
  gen_pT->Write();
  
  eta_response->Write();
  phi_response->Write();
  energy_response->Write();

  eb_energy_response->Write();
  ee_energy_response->Write();
  es_energy_response->Write();
  hcal_energy_response->Write();
  hcal_weight_energy_response->Write();
  hcal_detId_energy_response->Write();
  ecal_detId_energy_response->Write();

  scatter_x_z->Write();
  scatter_y_z->Write();
  scatter_x_y_z->Write();

  scatter_dR_depth->Write();
  scatter_dR_pT->Write();
  scatter_dR_pT_e->Write();

  scatter_pT_hits->Write();
  scatter_eta_hits->Write();
  scatter_pT_hitenergy->Write();

  eb_eta_hits->Write();
  ee_eta_hits->Write(); 
  es_eta_hits->Write(); 
  hcal_eta_hits->Write();

  eb_pt_hits->Write();
  ee_pt_hits->Write(); 
  es_pt_hits->Write(); 
  hcal_pt_hits->Write();

  eb_energy->Write(); 
  ee_energy->Write();  
  es_energy->Write(); 
  h_energy->Write();

  scatter_genpT_hitenergy->Write();   
  scatter_genenergy_hitenergy->Write(); 

  scatter_genpT_energy->Write();       
  scatter_genenergy_energy->Write();   


  relativCalenergy->Write(); 
  genenergy_relativCalenergy->Write();

  barrel_9gev_energy->Write();


  ecal_detId_genenergy->Write();
  hcal_detId_genenergy->Write();

  delete f;

}

// ------------ method called when starting to processes a run  ------------
/*
void 
SimHitResponse::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
SimHitResponse::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
SimHitResponse::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
SimHitResponse::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SimHitResponse::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimHitResponse);
