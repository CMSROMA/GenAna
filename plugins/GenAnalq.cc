// -*- C++ -*-
//
// Package:    MyGeneration/GenAnalq
// Class:      GenAnalq
// 
/**\class GenAnalq GenAnalq.cc MyGeneration/GenAnalq/plugins/GenAnalq.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Francesco Santanastasio
//         Created:  Tue, 31 Oct 2017 17:49:46 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"

//
// class declaration
//

#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <cassert>
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TF1.h"


// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

using namespace std;
using namespace edm;
using namespace reco;

class GenAnalq : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit GenAnalq(const edm::ParameterSet&);
      ~GenAnalq();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      int DefineBranches();

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      void Initialize();

      // ----------member data ---------------------------

  edm::Service<TFileService> fs_;
  
  edm::EDGetTokenT<reco::GenParticleCollection> genPartInputToken_;
  edm::EDGetTokenT<reco::GenJetCollection> genJetsAK4InputToken_;
  edm::EDGetTokenT<reco::GenJetCollection> genJetsAK8InputToken_;

  TH1F *h1_LQ_mass;

  TTree *outTree_;
  float LQ_mass_;

  int id_p = 2212;
  int id_LQ = 9911561;
  int id_photon = 22;
  int id_e = 11;
  int id_mu = 13;
  int id_tau = 15;
  int id_d = 1;
  int id_u = 2;
  int id_s = 3;
  int id_c = 4;
  int id_b = 5;
  int id_t = 6;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenAnalq::GenAnalq(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  usesResource("TFileService");
  genPartInputToken_ = (consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticlesInputTag")));
  genJetsAK4InputToken_ = (consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJetsAK4InputTag")));
  genJetsAK8InputToken_ = (consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJetsAK8InputTag")));
}


GenAnalq::~GenAnalq()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void
GenAnalq::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  Initialize();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif

   edm::Handle<reco::GenParticleCollection> genParticles;
   if (!iEvent.isRealData())
     iEvent.getByToken(genPartInputToken_, genParticles);

   edm::Handle<reco::GenJetCollection> genJetsAK4;
   if (!iEvent.isRealData())
     iEvent.getByToken(genJetsAK4InputToken_,genJetsAK4);

   edm::Handle<reco::GenJetCollection> genJetsAK8;
   if (!iEvent.isRealData())
     iEvent.getByToken(genJetsAK8InputToken_,genJetsAK8);

   TLorentzVector Gamma, LQ, LepIni, QuarkIni, LepFin, QuarkFin, Lep;
   // l q --> LQ --> l q 
   //Gamma = Photon emitted from one proton (missing in some events)
   //LepFromGamma = Final state electron from initial state Gamma 
   //LQ = Leptoquark
   //LepIni = Initial state lepton
   //QuarkIni = Initial state quark
   //LepFin = Final state lepton (daughter of LQ)
   //QuarkFin = Final state quark (daughter of LQ)

   if( genParticles.isValid() ) {

     for( reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it ) 
       {

	 //LQ
	 if (it->pdgId()==id_LQ)
	   {
	     LQ.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());	 
	   }

	 //Gamma
	 if (it->pdgId()==id_photon && it->mother()->pdgId()==id_p)
	   {
	     Gamma.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());	 
	   }

	 //Final state lepton
	 if ( (fabs(it->pdgId())==id_e || fabs(it->pdgId())==id_mu || fabs(it->pdgId())==id_tau)  && it->mother()->pdgId()==id_LQ)
	   {
	     LepFin.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());	 
	   }

	 //Final state quark
	 if ( ( fabs(it->pdgId())==id_u || fabs(it->pdgId())==id_c || fabs(it->pdgId())==id_t 
		|| fabs(it->pdgId())==id_d || fabs(it->pdgId())==id_b || fabs(it->pdgId())==id_s )
	      && it->mother()->pdgId()==id_LQ ) 
	   {
	     QuarkFin.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());	 
	   }
	 
       }

     //Fill histograms     
     h1_LQ_mass->Fill(LQ.M());

     //Fill tree variables
     LQ_mass_ = LQ.M();
    
   }

   //---- Fill Tree ---
   //add if statement here to filter events
   outTree_->Fill();     
   //------------------

}//end analyze for each event

// ------------ method called once each job just before starting event loop  ------------
void 
GenAnalq::beginJob()
{  
  // Book the histograms
  h1_LQ_mass = fs_->make<TH1F>("h1_LQ_mass","h1_LQ_mass",1000,0,10000);

  // Book the tree and define branches
  outTree_ = fs_->make<TTree>("events","events");
  DefineBranches();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenAnalq::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenAnalq::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

int GenAnalq::DefineBranches()
{
  //outTree_->Branch("run"                  ,&run_               ,"run_/I");
  //outTree_->Branch("met"                  ,&met_               ,"met_/F");
  //outTree_->Branch("gen_eta","vector<float>" , &gen_eta      );

  outTree_->Branch("LQ_mass"                  ,&LQ_mass_               ,"LQ_mass_/F");

  return 0;
}

void GenAnalq::Initialize()
{
  LQ_mass_            = -999;
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenAnalq);
