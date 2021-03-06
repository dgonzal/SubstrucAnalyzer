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
  
  ResponseBinning_ = ps.getParameter<vector<double>>("ResponseBinning");
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
  double samplingFactors[] = {125.44, 125.54, 125.32, 125.13, 
			      124.46, 125.01, 125.22, 125.48, 
			      124.45, 125.9, 125.83, 127.01, 
			      126.82, 129.73, 131.83, 143.52,
			      210.55, 197.93, 186.12, 189.64, 
			      189.63, 190.28, 189.61, 189.6, 
			      190.12, 191.22, 190.9, 193.06, 
			      188.42, 188.42};

  count_events+=1;


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
    //std::cout<<"Event "<<count_events <<" particleID "<<gen.pdgId() <<" gen pT "<< gen.pt()<<" phi "<<gen.phi()<< " eta "<<gen.eta()<< " e "<< gen.energy()<<std::endl;
    gen_energy->Fill(gen.energy());
    gen_pT->Fill(gen.pt());
    gen_eta->Fill(gen.eta());
    gen_phi->Fill(gen.phi());
    gen_pdgId->Fill(gen.pdgId());
    

  }


  //cout<<"NSimHits ECALEB "<<ECALEBSimHits->size()<< "  ECALES "<<ECALESSimHits->size()<< "  ECALEE "<<ECALEESimHits->size()<< " HCAL " <<HCALSimHits->size() <<endl;

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


  TH2F* phi_eta_vec = new TH2F("vec_phi_eta","#phi #eta",100,-3.5,3.5,100,-4,4);
  TH2F* delPhi_delEta_vec = new TH2F("vec_delphi_deleta","#Delta #phi, #Delta #eta",100,0,5,100,0,5);
  TH2F* energy_delR_vec = new TH2F("vec_energy_delR","Energy #Delta R" ,100,0,2,100,0,2);
  TH3F* x_y_z_vec = new TH3F("vec_x_y_z","x,y,z",100,-200,200,100,-200,200,100,-40,40);

  
  for(unsigned int p=0;p<4 &&count_events==-1 ;++p){
    edm::Handle<edm::PCaloHitContainer> SimHits;
    
    switch(p){
    case 0 : 
      SimHits = ECALEBSimHits;break;
    case 1 :
      SimHits = ECALEESimHits;break;
    case 2 :
      SimHits = ECALESSimHits;break;
    case 3 :
      SimHits = HCALSimHits;break;
      
    }
    
    for(unsigned int m =0; m<SimHits->size(); ++m){
      Hit = SimHits->at(m);
      pos = geo->getPosition(Hit.id());
      

      phi_eta_vec->Fill(pos.eta(),pos.phi());
      delPhi_delEta_vec->Fill(deltaR(0,pos.phi(),0,GenParticles->at(0).phi()),deltaR(pos.eta(),0,GenParticles->at(0).eta(),0));
      x_y_z_vec->Fill(pos.x(),pos.y(),pos.z());


      if(p==3){
	HcalDetId hcalId(Hit.id());

	cout<<m<<" Energy " <<Hit.energy()*samplingFactors[hcalId.ietaAbs()-1]<<" DeltaR "<<deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi())<<endl;

	energy_delR_vec->Fill( Hit.energy()*samplingFactors[hcalId.ietaAbs()-1],deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()));
      }
      else 
	energy_delR_vec->Fill( Hit.energy(),deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()));	
    }
  }
  if(count_events==61){
    vec_phi_eta.push_back(phi_eta_vec);
    vec_delPhi_delEta.push_back(delPhi_delEta_vec);
    vec_x_y_z.push_back( x_y_z_vec);
    vec_energy_delR.push_back(energy_delR_vec);
  }

  vector<sub_det> ecal_energy_depo;


  for(unsigned int i = 0; i<ECALEBSimHits->size(); ++i){
    Hit = ECALEBSimHits->at(i);
    pos = geo->getPosition(Hit.id());
    depth = sqrt(pos.x()*pos.x()+pos.y()*pos.y());
    
    //if(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi())>0.5) continue;

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

    simhit_energy_response->Fill(Hit.energy());
    eb_energy_response->Fill(Hit.energy());

    ecaleb_energy+= Hit.energy();

    scatter_genpT_hitenergy->Fill(GenParticles->at(0).pt(),Hit.energy() );
    scatter_genenergy_hitenergy->Fill(GenParticles->at(0).energy(), Hit.energy());

  }

 
  for(unsigned int i = 0; i<ECALEESimHits->size(); ++i){
    Hit = ECALEESimHits->at(i);
    pos = geo->getPosition(Hit.id());
    depth = sqrt(pos.x()*pos.x()+pos.y()*pos.y());
    
    //if(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi())>0.5) continue;


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

    simhit_energy_response->Fill(Hit.energy());
    ee_energy_response->Fill(Hit.energy());

    scatter_genpT_hitenergy->Fill(GenParticles->at(0).pt(),Hit.energy());
    scatter_genenergy_hitenergy->Fill(GenParticles->at(0).energy(), Hit.energy());


  }

  for(unsigned int i = 0; i<ECALESSimHits->size(); ++i){
    Hit = ECALESSimHits->at(i);
    pos = geo->getPosition(Hit.id());
    depth = sqrt(pos.x()*pos.x()+pos.y()*pos.y());
        
    //if(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi())>0.5) continue;


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

    simhit_energy_response->Fill(Hit.energy());
    es_energy_response->Fill(Hit.energy());



    scatter_genpT_hitenergy->Fill(GenParticles->at(0).pt(),Hit.energy());
    scatter_genenergy_hitenergy->Fill(GenParticles->at(0).energy(), Hit.energy());

  }


  vector<sub_det> det_energy_depo;

  
  for (unsigned int i = 0; i< HCALSimHits->size(); ++i){
    Hit = HCALSimHits->at(i);
    pos = geo->getPosition(Hit.id());
    depth = sqrt(pos.x()*pos.x()+pos.y()*pos.y());
    /*
    if(count_events==10)
      cout<<deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi())<< " ? "<<bool(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi())>0.5)<<endl;
    */
    //if(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi())>0.5) continue;


    bool detId_flag =true;

    
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


    if(max_posX<pos.x()) max_posX=pos.x();
    if(max_posY<pos.y()) max_posY=pos.y();

    scatter_dR_depth->Fill(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()),depth);
    scatter_dR_pT->Fill(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()),GenParticles->at(0).pt());
    scatter_dR_pT_e->Fill(deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()),GenParticles->at(0).pt(),Hit.energy());

    scatter_genpT_hitenergy->Fill(GenParticles->at(0).pt(), weighted_energy);
    scatter_genenergy_hitenergy->Fill(GenParticles->at(0).energy(), weighted_energy);

    eta_response->Fill( pos.eta());
    phi_response->Fill( pos.phi());
    simhit_energy_response->Fill(weighted_energy);
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
    pos = geo->getPosition(ecal_energy_depo.at(m).detId);
    ecal_detId_energy_response->Fill(ecal_energy_depo.at(m).energy);
    ecal_detId_genenergy->Fill(GenParticles->at(0).energy(),ecal_energy_depo.at(m).energy);
    ecal_detId_deltaR->Fill(GenParticles->at(0).energy(),deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()));
    ecal_detId_energy_deltaR->Fill(ecal_energy_depo.at(m).energy,deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()));
    ecal_detId_phi_eta->Fill(pos.phi(),pos.eta());
  }

  ecal_detId_energy_Num->Fill(GenParticles->at(0).energy(),ecal_energy_depo.size());


  for(unsigned int m=0; m<det_energy_depo.size();++m){
    pos = geo->getPosition(det_energy_depo.at(m).detId);
    hcal_detId_energy_response->Fill(det_energy_depo.at(m).energy);
    hcal_detId_genenergy->Fill(GenParticles->at(0).energy(),det_energy_depo.at(m).energy);
    hcal_detId_deltaR->Fill(GenParticles->at(0).energy(),deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()));
    hcal_detId_energy_deltaR->Fill(det_energy_depo.at(m).energy,deltaR(pos.eta(),pos.phi(),GenParticles->at(0).eta(),GenParticles->at(0).phi()));
    hcal_detId_phi_eta->Fill(pos.phi(),pos.eta());
  }

  hcal_detId_energy_Num->Fill(GenParticles->at(0).energy(),det_energy_depo.size());

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
  
 /*
  std::cout<<"hcal energy "<<hcal_energy<<" corrected "<<hcal_energy_weight<< " ecal EB energy "<< ecaleb_energy<< " ecal ES energy "<< ecales_energy << " ecal EE energy "<< ecalee_energy << std::endl;

  std::cout<<" max (x,y) "<<max_posX<<" "<<max_posY<<std::endl;

 
  if(hcal_energy_weight+ecaleb_energy+ecales_energy+ecalee_energy<1){ 
    std::cout<<"------------------------"<<std::endl;
    std::cout<<"Event "<<count_events<<" particleID "<< GenParticles->at(0).pdgId() <<" gen pT "<< GenParticles->at(0).pt()<<" phi "<<GenParticles->at(0).phi()<< " eta "<<GenParticles->at(0).eta()<< " e "<< GenParticles->at(0).energy()<<std::endl;
    cout<<"NSimHits ECALEB "<<ECALEBSimHits->size()<< "  ECALES "<<ECALESSimHits->size()<< "  ECALEE "<<ECALEESimHits->size()<< " HCAL " <<HCALSimHits->size()<<" energy "<<hcal_energy_weight+ecaleb_energy+ecales_energy+ecalee_energy <<endl;
    std::cout<<"------------------------"<<std::endl;
  }
  */
  
  energy_response->Fill(hcal_energy_weight+ecaleb_energy+ecales_energy+ecalee_energy);

 
  if(ECALEBSimHits->size()>0)eb_energy->Fill(ecaleb_energy);
  if(ECALEESimHits->size()>0)ee_energy->Fill(ecalee_energy);
  if(ECALESSimHits->size()>0)es_energy->Fill(ecales_energy);
  if(HCALSimHits->size()>0)h_energy->Fill(hcal_energy);

  delete phi_eta_vec;
  delete delPhi_delEta_vec; 
  delete energy_delR_vec;
  delete x_y_z_vec; 





}


// ------------ method called once each job just before starting event loop  ------------
void SimHitResponse::beginJob(){
  count_events=0;

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
  simhit_energy_response = new TH1F("simhit_energy_response","E",150,0,0.002);

  eb_energy_response = new TH1F("eb_energy_response","E",150,0,2);//0.002
  ee_energy_response = new TH1F("ee_energy_response","E",150,0,2);
  es_energy_response = new TH1F("es_energy_response","E",150,0,2);
  hcal_energy_response = new TH1F("hcal_energy_response","E",150,0,2);
  hcal_weight_energy_response = new TH1F("hcal_weight_energy_response","E",150,0,2);
  hcal_detId_energy_response = new TH1F("hcal_detId_energy_response","E",250,0,10);
  ecal_detId_energy_response = new TH1F("ecal_detId_energy_response","E",250,0,10);

  hcal_detId_genenergy = new TH2F("hcal_detId_genenergy","E_{DetId}, E_{gen}",150,0,20,250,0,10);
  ecal_detId_genenergy = new TH2F("ecal_detId_genenergy","E_{DetId}, E_{gen}",150,0,20,250,0,10);

  hcal_detId_energy_deltaR = new TH2F("hcal_detId_energy_deltaR","E_{DetId}, #Delta R",150,0,10,250,0,10);
  ecal_detId_energy_deltaR = new TH2F("ecal_detId_energy_deltaR","E_{DetId}, #Delta R",150,0,10,250,0,10);

  ecal_detId_deltaR = new TH2F("ecal_detId_deltaR","E_{Gen}, #Delta R",150,0,20,250,0,6);
  hcal_detId_deltaR = new TH2F("hcal_detId_deltaR","E_{Gen}, #Delta R",150,0,20,250,0,6);

  ecal_detId_phi_eta = new TH2F("ecal_detId_phi_eta","#phi, #eta",150,-4,4,250,-4,4);
  hcal_detId_phi_eta = new TH2F("hcal_detId_phi_eta","#phi, #eta",150,-4,4,250,-4,4);

  ecal_detId_energy_Num = new TH2F("ecal_detId_energy_Num","Eg{Gen}, Number of DetId",150,0,20,250,0,100);
  hcal_detId_energy_Num = new TH2F("hcal_detId_energy_Num","E_{Gen}, Number of DetId",150,0,20,250,0,100);

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
  gen_eta->Write();
  gen_phi->Write();
  gen_pdgId->Write();

  eta_response->Write();
  phi_response->Write();
  simhit_energy_response->Write();

  eb_energy_response->Write();
  ee_energy_response->Write();
  es_energy_response->Write();
  hcal_energy_response->Write();
  hcal_weight_energy_response->Write();
  hcal_detId_energy_response->Write();
  ecal_detId_energy_response->Write();

  ecal_detId_deltaR->Write();
  hcal_detId_deltaR->Write();

  hcal_detId_energy_deltaR->Write();
  ecal_detId_energy_deltaR->Write();

  ecal_detId_phi_eta->Write();
  hcal_detId_phi_eta->Write(); 

  ecal_detId_energy_Num->Write();
  hcal_detId_energy_Num->Write();

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

  energy_response->Write();


  ecal_detId_genenergy->Write();
  hcal_detId_genenergy->Write();



  //std::cout<< vec_phi_eta.size()<<std::endl;
  /*
  //single shower
  for(unsigned int i = 0; i< vec_phi_eta.size() &&vec_phi_eta.size()>0 ;i++){
    stringstream dir_name;
    dir_name<< "Single_Shower_"<<i+1;
    f->mkdir(dir_name.str().c_str());
    f->cd(dir_name.str().c_str());

    vec_phi_eta[i]->Write();
    vec_delPhi_delEta[i]->Write();
    vec_x_y_z[i]->Write();
    vec_energy_delR[i]->Write();
  }
  */
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
