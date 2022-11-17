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
#include "TH2F.h"
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

  TH1F *h1_LQ_mass, *h1_lep_pt, *h1_lep_eta, *h1_q_pt, *h1_q_eta, *h1_gamma_pt, *h1_gamma_pt__pos, *h1_lepFromGamma_pt, *h1_lepFromGamma_eta;
  TH1F *h1_xip, *h1_xiLQ, *h1_xip_over_xiLQ_minus_one, *h1_xip_over_xiLQ_minus_one__inPPS, *h1_xip_over_xiLQ_minus_one__plusLep, *h1_xip_over_xiLQ_minus_one__inPPS__plusLep; 
  TH2F *h2_lqIn__vs__lqOut_mass, *h2_xip__vs__xiLQ;
  TH1F *h1_LQreco_mass, *h1_jet_pt, *h1_jet_eta;

  TTree *outTree_;
  int event_;
  float LQ_E_,LQ_mass_, LQ_pt_, LQ_eta_, LQ_phi_, LQ_id_, LQ_px_, LQ_py_, LQ_pz_;
  float gamma_E_,gamma_mass_, gamma_pt_, gamma_eta_, gamma_phi_, gamma_id_, gamma_from_p_, gamma_px_, gamma_py_, gamma_pz_;
  float jet1_mass_, jet1_pt_, jet1_eta_, jet1_phi_;
  float lepFromGamma_E_, lepFromGamma_mass_, lepFromGamma_pt_, lepFromGamma_eta_, lepFromGamma_phi_, lepFromGamma_id_;
  float pIn_E_, pIn_mass_, pIn_px_, pIn_py_, pIn_pz_, pIn_id_;
  float pOut_E_, pOut_mass_, pOut_px_, pOut_py_, pOut_pz_, pOut_id_;
  float xi_p_, xi_LQ_, xi_LQ_wrong_, xi_LQ_plusLep_, xi_LQ_wrong_plusLep_;
  float lepIn_E_,lepIn_mass_, lepIn_pt_, lepIn_eta_, lepIn_phi_, lepIn_id, lepIn_id_, lepIn_px_, lepIn_py_, lepIn_pz_;
  float qIn_E_,qIn_mass_, qIn_pt_, qIn_eta_, qIn_phi_, qIn_id_, qIn_px_, qIn_py_, qIn_pz_;
  float lepOut_E_,lepOut_mass_, lepOut_pt_, lepOut_eta_, lepOut_phi_, lepOut_id_;
  float qOut_E_,qOut_mass_, qOut_pt_, qOut_eta_, qOut_phi_, qOut_id_;
  float lqIn_mass_, lqOut_mass_;

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
  int status_intermediate = 11;

  float max_xi_PPS = 0.2;

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

   TLorentzVector LQ, LQplusLep, pIn, pOut, Gamma, LepFromGamma, LepIn, QuarkIn, LepOut, QuarkOut;
   // p --> p+Gamma , Gamma -->ll ,  l q --> LQ --> l q 
   //Gamma = Photon emitted from one proton (missing in some events?)
   //LepFromGamma = Final state lepton from initial state Gamma 
   //LQ = Leptoquark
   //LepIn = Incoming lepton
   //QuarkIn = Incoming quark
   //LepOut = Outgoing lepton (daughter of LQ)
   //QuarkOut = Outgoing quark (daughter of LQ)

   event_ = iEvent.id().event();
   //cout << "----------- Event: " << event_ << endl;

   if( genParticles.isValid() ) {

     int nLQ = 0;
     int nLepGamma = 0;
     int nGamma = 0;
     int nLepIn = 0;
     int nQuarkIn = 0;
     int nLepOut = 0;
     int nQuarkOut = 0;
     int gammaFromProton = -9;

     for( reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it ) 
       {
	 
	 //LQ
	 if (it->pdgId()==id_LQ)
	   {
	     //LQ.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());	 	     
	     LQ.SetPxPyPzE(it->px(),it->py(),it->pz(),it->energy());	 	     
	     nLQ++;
	     
	     LQ_E_ = LQ.E(); 
	     LQ_mass_ = LQ.M(); 
	     LQ_pt_ = LQ.Pt(); 
	     LQ_eta_ = LQ.Eta(); 
	     LQ_phi_ = LQ.Phi();  
	     LQ_id_ = it->pdgId();  
	     LQ_px_ = it->px();  
	     LQ_py_ = it->py();  
	     LQ_pz_ = it->pz();  
	   }
	 
	 //Gamma and Proton
	 if (it->pdgId()==id_photon 
	     && it->status()==status_intermediate 
	     && ( it->mother()->pdgId()==id_p || 
		  ( fabs(it->mother()->pdgId())==id_u || fabs(it->mother()->pdgId())==id_c || fabs(it->mother()->pdgId())==id_t 
		    || fabs(it->mother()->pdgId())==id_d || fabs(it->mother()->pdgId())==id_s || fabs(it->mother()->pdgId())==id_b ) 
		  )
	     )
	   {
	     //Gamma.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());	 
	     Gamma.SetPxPyPzE(it->px(),it->py(),it->pz(),it->energy());	 	     

	     //Find the proton that emitted the gamma
	     for( reco::GenParticleCollection::const_iterator myit = genParticles->begin(); myit != genParticles->end(); ++myit ) 
	       {
		 if(myit->pdgId()==id_p && it->pz()*myit->pz()>0.)
		   {
		     pIn.SetPxPyPzE(myit->px(),myit->py(),myit->pz(),myit->energy()); 	 	     
		     break;
		   }
	       }
	     //if(it->pz()>0)
	     //  pIn.SetPxPyPzE(0.,0.,+6500,sqrt(6500*6500+0.938*0.938)); //proton +z LHC 13 TeV	 	     
	     //if(it->pz()<0)
	     //  pIn.SetPxPyPzE(0.,0.,-6500,sqrt(6500*6500+0.938*0.938)); //proton -z LHC 13 TeV	     
	     nGamma++;
	     
	     gamma_E_ = Gamma.E(); 
	     gamma_mass_ = Gamma.M(); 
	     gamma_pt_ = Gamma.Pt(); 
	     gamma_eta_ = Gamma.Eta(); 
	     gamma_phi_ = Gamma.Phi();  
	     gamma_id_ = it->pdgId();  
	     gamma_px_ = it->px();
	     gamma_py_ = it->py();
	     gamma_pz_ = it->pz();   
	     
	     if(it->mother()->pdgId()==id_p) //gamma from proton
	       {
		 gamma_from_p_ = 1;
		 gammaFromProton = 1;
	       }
	     if( ( fabs(it->mother()->pdgId())==id_u || fabs(it->mother()->pdgId())==id_c || fabs(it->mother()->pdgId())==id_t 
		   || fabs(it->mother()->pdgId())==id_d || fabs(it->mother()->pdgId())==id_s || fabs(it->mother()->pdgId())==id_b )
		 ) //gamma from quark
	       {
		 gamma_from_p_ = 0;
		 gammaFromProton = 0;
	       }
	     
	     pIn_E_ = pIn.E(); 
	     pIn_mass_ = pIn.M();  
	     pIn_px_ = pIn.Px();  
	     pIn_py_ = pIn.Py();  
	     pIn_pz_ = pIn.Pz();   
	     pIn_id_ = id_p;  
	     
	     pOut = pIn - Gamma;
	     pOut_E_ = pOut.E(); 
	     pOut_mass_ = pOut.M(); 
	     pOut_px_ = pOut.Px(); 
	     pOut_py_ = pOut.Py(); 
	     pOut_pz_ = pOut.Pz(); 
	     pOut_id_ = pIn_id_; 
	   }
	 
	 //Final state lepton from Gamma	 
	 if ( (fabs(it->pdgId())==id_e || fabs(it->pdgId())==id_mu || fabs(it->pdgId())==id_tau) and it->status()==1 )
	   {

	     if(it->numberOfMothers()>0)
	       {
		 const Candidate * thismother = it->mother();
		 
		 while(true)
		   {
		     if(thismother->numberOfMothers()==0)
		       break;

		     if(thismother->pdgId()==id_LQ)
		       break;
		     
		     if (thismother->pdgId()==id_photon 
			 && thismother->status()==status_intermediate 
			 && ( thismother->mother()->pdgId()==id_p || 
			      ( fabs(thismother->mother()->pdgId())==id_u || fabs(thismother->mother()->pdgId())==id_c || fabs(thismother->mother()->pdgId())==id_t 
				|| fabs(thismother->mother()->pdgId())==id_d || fabs(thismother->mother()->pdgId())==id_s || fabs(thismother->mother()->pdgId())==id_b ) )  
			 )
		       {
			 //lepton from gamma found
			 nLepGamma++;

			 LepFromGamma.SetPxPyPzE(it->px(),it->py(),it->pz(),it->energy());	 	     

			 lepFromGamma_E_ = LepFromGamma.E(); 
			 lepFromGamma_mass_ = LepFromGamma.M(); 
			 lepFromGamma_pt_ = LepFromGamma.Pt(); 
			 lepFromGamma_eta_ = LepFromGamma.Eta(); 
			 lepFromGamma_phi_ = LepFromGamma.Phi();  
			 lepFromGamma_id_ = it->pdgId();  

			 break;
		       }
		     else
		       {
			 thismother = thismother->mother();
		       }
		   }

	       }

	   }
	 
	 //Incoming lepton
	 if ( (fabs(it->pdgId())==id_e || fabs(it->pdgId())==id_mu || fabs(it->pdgId())==id_tau) )
	   {

	     int dauId = -999;
	     if( it->numberOfDaughters()>0 )
	       {
		 const Candidate * d = it->daughter( 0 );
		 dauId = d->pdgId();
	       }

	     if(dauId==id_LQ)
	       {
		 //LepIn.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());	 
		 LepIn.SetPxPyPzE(it->px(),it->py(),it->pz(),it->energy());	 	     
		 nLepIn++;
		 
		 lepIn_E_ = it->energy(); 
		 lepIn_mass_ = it->mass(); 
		 lepIn_pt_ = LepIn.Pt(); 
		 lepIn_eta_ = LepIn.Eta(); 
		 lepIn_phi_ = LepIn.Phi();  
		 lepIn_id_ = it->pdgId();  
		 lepIn_px_ = it->px(); 
		 lepIn_py_ = it->py(); 
		 lepIn_pz_ = it->pz(); 
	       }
	   }

	 //Incoming quark
	 if ( ( fabs(it->pdgId())==id_u || fabs(it->pdgId())==id_c || fabs(it->pdgId())==id_t 
		|| fabs(it->pdgId())==id_d || fabs(it->pdgId())==id_s || fabs(it->pdgId())==id_b ) ) 
	   {
	     int dauId = -999;
	     if( it->numberOfDaughters()>0 )
	       {
		 const Candidate * d = it->daughter( 0 );
		 dauId = d->pdgId();
	       }

	     if(dauId==id_LQ)
	       {
		 //QuarkIn.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());	 
		 QuarkIn.SetPxPyPzE(it->px(),it->py(),it->pz(),it->energy());	 	     
		 nQuarkIn++;
		 
		 qIn_E_ = it->energy(); 
		 qIn_mass_ = it->mass(); 
		 qIn_pt_ = QuarkIn.Pt(); 
		 qIn_eta_ = QuarkIn.Eta(); 
		 qIn_phi_ = QuarkIn.Phi();  
		 qIn_id_ = it->pdgId();  
		 qIn_px_ = it->px();  
		 qIn_py_ = it->py();  
		 qIn_pz_ = it->pz();  
	       }
	   }

	 //Outgoing lepton
	 if ( (fabs(it->pdgId())==id_e || fabs(it->pdgId())==id_mu || fabs(it->pdgId())==id_tau)  && it->mother()->pdgId()==id_LQ)
	   {
	     //LepOut.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());	 
	     LepOut.SetPxPyPzE(it->px(),it->py(),it->pz(),it->energy());	 	     
	     nLepOut++;

	     lepOut_E_ = LepOut.E(); 
	     lepOut_mass_ = LepOut.M(); 
	     lepOut_pt_ = LepOut.Pt(); 
	     lepOut_eta_ = LepOut.Eta(); 
	     lepOut_phi_ = LepOut.Phi();  
	     lepOut_id_ = it->pdgId();  
	   }

	 //Outgoing quark
	 if ( ( fabs(it->pdgId())==id_u || fabs(it->pdgId())==id_c || fabs(it->pdgId())==id_t 
		|| fabs(it->pdgId())==id_d || fabs(it->pdgId())==id_b || fabs(it->pdgId())==id_s )
	      && it->mother()->pdgId()==id_LQ ) 
	   {
	     //QuarkOut.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());	 
	     QuarkOut.SetPxPyPzE(it->px(),it->py(),it->pz(),it->energy());	 	     
	     nQuarkOut++;

	     qOut_E_ = QuarkOut.E(); 
	     qOut_mass_ = QuarkOut.M(); 
	     qOut_pt_ = QuarkOut.Pt(); 
	     qOut_eta_ = QuarkOut.Eta(); 
	     qOut_phi_ = QuarkOut.Phi();  
	     qOut_id_ = it->pdgId();   
	   }

	 if(nGamma>0 && nLepGamma>0)
	   {

	     LQplusLep = LQ + LepFromGamma;

	     xi_p_ = 1 - ( pOut.Pz() / pIn.Pz() ) ; 
	     if( pOut.Pz() > 0 )
	       {
		 xi_LQ_ = ( LQ.E() + LQ.Pz() ) / fabs( 2*pIn.Pz() ) ; 
		 xi_LQ_wrong_ = ( LQ.E() - LQ.Pz() ) / fabs( 2*pIn.Pz() ) ; 

		 xi_LQ_plusLep_ = ( LQplusLep.E() + LQplusLep.Pz() ) / fabs( 2*pIn.Pz() ) ; 
		 xi_LQ_wrong_plusLep_ = ( LQplusLep.E() - LQplusLep.Pz() ) / fabs( 2*pIn.Pz() ) ; 				 
	       } 
	     else
	       {
		 xi_LQ_ = ( LQ.E() - LQ.Pz() ) / fabs( 2*pIn.Pz() ) ; 
		 xi_LQ_wrong_ = ( LQ.E() + LQ.Pz() ) / fabs( 2*pIn.Pz() ) ; 

		 xi_LQ_plusLep_ = ( LQplusLep.E() - LQplusLep.Pz() ) / fabs( 2*pIn.Pz() ) ; 
		 xi_LQ_wrong_plusLep_ = ( LQplusLep.E() + LQplusLep.Pz() ) / fabs( 2*pIn.Pz() ) ; 				 
	       }
	   }

	 //Exit from loop when particles are found
	 if(nLQ>0 && nLepGamma>0 && nGamma>0 && nLepIn>0 && nQuarkIn>0 && nLepOut>0 && nQuarkOut>0)
	   break;	   
	 
       }// end loop over gen particles


     //Fill tree variables
     lqIn_mass_ = (LepIn + QuarkIn).M();
     lqOut_mass_ = (LepOut + QuarkOut).M();

     //Fill histograms     
     h1_LQ_mass->Fill(LQ.M());
     h1_lep_pt->Fill(LepOut.Pt());
     h1_lep_eta->Fill(LepOut.Eta());
     h1_q_pt->Fill(QuarkOut.Pt());
     h1_q_eta->Fill(QuarkOut.Eta());
     h1_gamma_pt->Fill(gamma_pt_);
     
     h2_lqIn__vs__lqOut_mass->Fill(lqOut_mass_,lqIn_mass_);

     //Fill only if photon is coming from proton
     if(nGamma>0 && nLepGamma>0 && gammaFromProton==1)
       {
	 h1_gamma_pt__pos->Fill(Gamma.Pt());
	 h1_lepFromGamma_pt->Fill(LepFromGamma.Pt());
	 h1_lepFromGamma_eta->Fill(LepFromGamma.Eta());
	 h1_xip->Fill(xi_p_);
	 h1_xiLQ->Fill(xi_LQ_);

	 h1_xip_over_xiLQ_minus_one->Fill( (xi_p_/xi_LQ_) - 1 );
	 if(xi_p_<max_xi_PPS)
	   h1_xip_over_xiLQ_minus_one__inPPS->Fill( (xi_p_/xi_LQ_) - 1 );

	 h1_xip_over_xiLQ_minus_one__plusLep->Fill( (xi_p_/xi_LQ_plusLep_) - 1 );
	 if(xi_p_<max_xi_PPS)
	   h1_xip_over_xiLQ_minus_one__inPPS__plusLep->Fill( (xi_p_/xi_LQ_plusLep_) - 1 );

	 h2_xip__vs__xiLQ->Fill(xi_LQ_,xi_p_);	 
       }


     if( genJetsAK8.isValid() && genJetsAK8->size()>=2 ) 
       {

	 TLorentzVector AK8Jet1, AK8Jet2, LQreco;
	 AK8Jet1.SetPtEtaPhiM(genJetsAK8->at(0).pt(),genJetsAK8->at(0).eta(),genJetsAK8->at(0).phi(),genJetsAK8->at(0).mass());
	 AK8Jet2.SetPtEtaPhiM(genJetsAK8->at(1).pt(),genJetsAK8->at(1).eta(),genJetsAK8->at(1).phi(),genJetsAK8->at(1).mass());

	 if( AK8Jet1.DeltaR(LepOut) > AK8Jet2.DeltaR(LepOut) )
	   {
	     jet1_mass_ = AK8Jet1.M();
	     jet1_pt_ = AK8Jet1.Pt();
	     jet1_eta_ = AK8Jet1.Eta();
	     jet1_phi_ = AK8Jet1.Phi();

	     LQreco = AK8Jet1 + LepOut;
	     
	     h1_LQreco_mass->Fill(LQreco.M());
	     h1_jet_pt->Fill(AK8Jet1.Pt());
	     h1_jet_eta->Fill(AK8Jet1.Eta());
	   }
	 else
	   {
	     jet1_mass_ = AK8Jet2.M();
	     jet1_pt_ = AK8Jet2.Pt();
	     jet1_eta_ = AK8Jet2.Eta();
	     jet1_phi_ = AK8Jet2.Phi();

	     LQreco = AK8Jet2 + LepOut;
	     
	     h1_LQreco_mass->Fill(LQreco.M());
	     h1_jet_pt->Fill(AK8Jet2.Pt());
	     h1_jet_eta->Fill(AK8Jet2.Eta());
	   }

       }

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
  h1_lep_pt = fs_->make<TH1F>("h1_lep_pt","h1_lep_pt",1000,0,10000);
  h1_lep_eta = fs_->make<TH1F>("h1_lep_eta","h1_lep_eta",100,-5,5);
  h1_q_pt = fs_->make<TH1F>("h1_q_pt","h1_q_pt",1000,0,10000);
  h1_q_eta = fs_->make<TH1F>("h1_q_eta","h1_q_eta",100,-5,5);
  h1_gamma_pt = fs_->make<TH1F>("h1_gamma_pt","h1_gamma_pt",300,-2000,1000);
  h1_gamma_pt__pos = fs_->make<TH1F>("h1_gamma_pt__pos","h1_gamma_pt__pos",50,0,10);
  h1_lepFromGamma_pt = fs_->make<TH1F>("h1_lepFromGamma_pt","h1_lepFromGamma_pt",200,0,2000);
  h1_lepFromGamma_eta = fs_->make<TH1F>("h1_lepFromGamma_eta","h1_lepFromGamma_eta",100,-10,10);
  h1_xip = fs_->make<TH1F>("h1_xip","h1_xip",100,0,1);
  h1_xiLQ = fs_->make<TH1F>("h1_xiLQ","h1_xiLQ",100,0,1);
  h1_xip_over_xiLQ_minus_one = fs_->make<TH1F>("h1_xip_over_xiLQ_minus_one","h1_xip_over_xiLQ_minus_one",2000,-10,10);
  h1_xip_over_xiLQ_minus_one__inPPS = fs_->make<TH1F>("h1_xip_over_xiLQ_minus_one__inPPS","h1_xip_over_xiLQ_minus_one__inPPS",2000,-10,10);
  h1_xip_over_xiLQ_minus_one__plusLep = fs_->make<TH1F>("h1_xip_over_xiLQ_minus_one__plusLep","h1_xip_over_xiLQ_minus_one__plusLep",2000,-10,10);
  h1_xip_over_xiLQ_minus_one__inPPS__plusLep = fs_->make<TH1F>("h1_xip_over_xiLQ_minus_one__inPPS__plusLep","h1_xip_over_xiLQ_minus_one__inPPS__plusLep",2000,-10,10);

  h2_lqIn__vs__lqOut_mass = fs_->make<TH2F>("h2_lqIn__vs__lqOut_mass","h2_lqIn__vs__lqOut_mass",500,-10000,10000,500,-30000,30000);
  h2_xip__vs__xiLQ = fs_->make<TH2F>("h2_xip__vs__xiLQ","h2_xip__vs__xiLQ",50,0,1,50,0,1);

  h1_LQreco_mass = fs_->make<TH1F>("h1_LQreco_mass","h1_LQreco_mass",1000,0,10000);
  h1_jet_pt = fs_->make<TH1F>("h1_jet_pt","h1_jet_pt",1000,0,10000);
  h1_jet_eta = fs_->make<TH1F>("h1_jet_eta","h1_jet_eta",100,-5,5);

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
  outTree_->Branch("event"                  ,&event_               ,"event_/I");
  //outTree_->Branch("gen_eta","vector<float>" , &gen_eta      );

  outTree_->Branch("LQ_E"                     ,&LQ_E_                  ,"LQ_E_/F");
  outTree_->Branch("LQ_mass"                  ,&LQ_mass_               ,"LQ_mass_/F");
  outTree_->Branch("LQ_pt"                    ,&LQ_pt_                 ,"LQ_pt_/F");
  outTree_->Branch("LQ_eta"                   ,&LQ_eta_                ,"LQ_eta_/F");
  outTree_->Branch("LQ_phi"                   ,&LQ_phi_                ,"LQ_phi_/F");
  outTree_->Branch("LQ_id"                    ,&LQ_id_                 ,"LQ_id_/F");
  outTree_->Branch("LQ_px"                    ,&LQ_px_                 ,"LQ_px_/F");
  outTree_->Branch("LQ_py"                    ,&LQ_py_                 ,"LQ_py_/F");
  outTree_->Branch("LQ_pz"                    ,&LQ_pz_                 ,"LQ_pz_/F");

  outTree_->Branch("gamma_E"                     ,&gamma_E_                  ,"gamma_E_/F");
  outTree_->Branch("gamma_mass"                  ,&gamma_mass_               ,"gamma_mass_/F");
  outTree_->Branch("gamma_pt"                    ,&gamma_pt_                 ,"gamma_pt_/F");
  outTree_->Branch("gamma_eta"                   ,&gamma_eta_                ,"gamma_eta_/F");
  outTree_->Branch("gamma_phi"                   ,&gamma_phi_                ,"gamma_phi_/F");
  outTree_->Branch("gamma_id"                    ,&gamma_id_                 ,"gamma_id_/F");
  outTree_->Branch("gamma_from_p_"               ,&gamma_from_p_             ,"gamma_from_p_/F");
  outTree_->Branch("gamma_px"                    ,&gamma_px_                 ,"gamma_px_/F");
  outTree_->Branch("gamma_py"                    ,&gamma_py_                 ,"gamma_py_/F");
  outTree_->Branch("gamma_pz"                    ,&gamma_pz_                 ,"gamma_pz_/F");

  outTree_->Branch("lepFromGamma_E"                     ,&lepFromGamma_E_                  ,"lepFromGamma_E_/F");
  outTree_->Branch("lepFromGamma_mass"                  ,&lepFromGamma_mass_               ,"lepFromGamma_mass_/F");
  outTree_->Branch("lepFromGamma_pt"                    ,&lepFromGamma_pt_                 ,"lepFromGamma_pt_/F");
  outTree_->Branch("lepFromGamma_eta"                   ,&lepFromGamma_eta_                ,"lepFromGamma_eta_/F");
  outTree_->Branch("lepFromGamma_phi"                   ,&lepFromGamma_phi_                ,"lepFromGamma_phi_/F");
  outTree_->Branch("lepFromGamma_id"                    ,&lepFromGamma_id_                 ,"lepFromGamma_id_/F");

  outTree_->Branch("pIn_E"                     ,&pIn_E_                  ,"pIn_E_/F");
  outTree_->Branch("pIn_mass"                  ,&pIn_mass_               ,"pIn_mass_/F");
  outTree_->Branch("pIn_px"                    ,&pIn_px_                 ,"pIn_px_/F");
  outTree_->Branch("pIn_py"                    ,&pIn_py_                 ,"pIn_py_/F");
  outTree_->Branch("pIn_pz"                    ,&pIn_pz_                 ,"pIn_pz_/F");
  outTree_->Branch("pIn_id"                    ,&pIn_id_                 ,"pIn_id_/F");

  outTree_->Branch("pOut_E"                     ,&pOut_E_                  ,"pOut_E_/F");
  outTree_->Branch("pOut_mass"                  ,&pOut_mass_               ,"pOut_mass_/F");
  outTree_->Branch("pOut_px"                    ,&pOut_px_                 ,"pOut_px_/F");
  outTree_->Branch("pOut_py"                    ,&pOut_py_                 ,"pOut_py_/F");
  outTree_->Branch("pOut_pz"                    ,&pOut_pz_                 ,"pOut_pz_/F");
  outTree_->Branch("pOut_id"                    ,&pOut_id_                 ,"pOut_id_/F");

  outTree_->Branch("jet1_mass"                  ,&jet1_mass_               ,"jet1_mass_/F");
  outTree_->Branch("jet1_pt"                    ,&jet1_pt_                 ,"jet1_pt_/F");
  outTree_->Branch("jet1_eta"                   ,&jet1_eta_                ,"jet1_eta_/F");
  outTree_->Branch("jet1_phi"                   ,&jet1_phi_                ,"jet1_phi_/F");

  outTree_->Branch("xi_p"                       ,&xi_p_                    ,"xi_p_/F");
  outTree_->Branch("xi_LQ"                      ,&xi_LQ_                   ,"xi_LQ_/F");
  outTree_->Branch("xi_LQ_wrong"                ,&xi_LQ_wrong_             ,"xi_LQ_wrong_/F");
  outTree_->Branch("xi_LQ_plusLep"                      ,&xi_LQ_plusLep_                   ,"xi_LQ_plusLep_/F");
  outTree_->Branch("xi_LQ_wrong_plusLep"                ,&xi_LQ_wrong_plusLep_             ,"xi_LQ_wrong_plusLep_/F");

  outTree_->Branch("lepIn_E"                     ,&lepIn_E_                  ,"lepIn_E_/F");
  outTree_->Branch("lepIn_mass"                  ,&lepIn_mass_               ,"lepIn_mass_/F");
  outTree_->Branch("lepIn_pt"                    ,&lepIn_pt_                 ,"lepIn_pt_/F");
  outTree_->Branch("lepIn_eta"                   ,&lepIn_eta_                ,"lepIn_eta_/F");
  outTree_->Branch("lepIn_phi"                   ,&lepIn_phi_                ,"lepIn_phi_/F");
  outTree_->Branch("lepIn_id"                    ,&lepIn_id_                 ,"lepIn_id_/F");
  outTree_->Branch("lepIn_px"                    ,&lepIn_px_                 ,"lepIn_px_/F");
  outTree_->Branch("lepIn_py"                    ,&lepIn_py_                 ,"lepIn_py_/F");
  outTree_->Branch("lepIn_pz"                    ,&lepIn_pz_                 ,"lepIn_pz_/F");

  outTree_->Branch("qIn_E"                     ,&qIn_E_                  ,"qIn_E_/F");
  outTree_->Branch("qIn_mass"                  ,&qIn_mass_               ,"qIn_mass_/F");
  outTree_->Branch("qIn_pt"                    ,&qIn_pt_                 ,"qIn_pt_/F");
  outTree_->Branch("qIn_eta"                   ,&qIn_eta_                ,"qIn_eta_/F");
  outTree_->Branch("qIn_phi"                   ,&qIn_phi_                ,"qIn_phi_/F");
  outTree_->Branch("qIn_id"                    ,&qIn_id_                 ,"qIn_id_/F");
  outTree_->Branch("qIn_px"                    ,&qIn_px_                 ,"qIn_px_/F");
  outTree_->Branch("qIn_py"                    ,&qIn_py_                 ,"qIn_py_/F");
  outTree_->Branch("qIn_pz"                    ,&qIn_pz_                 ,"qIn_pz_/F");

  outTree_->Branch("lepOut_E"                     ,&lepOut_E_                  ,"lepOut_E_/F");
  outTree_->Branch("lepOut_mass"                  ,&lepOut_mass_               ,"lepOut_mass_/F");
  outTree_->Branch("lepOut_pt"                    ,&lepOut_pt_                 ,"lepOut_pt_/F");
  outTree_->Branch("lepOut_eta"                   ,&lepOut_eta_                ,"lepOut_eta_/F");
  outTree_->Branch("lepOut_phi"                   ,&lepOut_phi_                ,"lepOut_phi_/F");
  outTree_->Branch("lepOut_id"                    ,&lepOut_id_                 ,"lepOut_id_/F");

  outTree_->Branch("qOut_E"                     ,&qOut_E_                  ,"qOut_E_/F");
  outTree_->Branch("qOut_mass"                  ,&qOut_mass_               ,"qOut_mass_/F");
  outTree_->Branch("qOut_pt"                    ,&qOut_pt_                 ,"qOut_pt_/F");
  outTree_->Branch("qOut_eta"                   ,&qOut_eta_                ,"qOut_eta_/F");
  outTree_->Branch("qOut_phi"                   ,&qOut_phi_                ,"qOut_phi_/F");
  outTree_->Branch("qOut_id"                    ,&qOut_id_                 ,"qOut_id_/F");

  outTree_->Branch("lqIn_mass"                    ,&lqIn_mass_                 ,"lqIn_mass_/F");
  outTree_->Branch("lqOut_mass"                   ,&lqOut_mass_                ,"lqOut_mass_/F");

  return 0;
}

void GenAnalq::Initialize()
{
  event_             = -999;

  LQ_E_               = -999;
  LQ_mass_            = -999;
  LQ_pt_              = -999; 
  LQ_eta_             = -999; 
  LQ_phi_             = -999; 
  LQ_id_              = -999; 
  LQ_px_              = -999; 
  LQ_py_              = -999;
  LQ_pz_              = -999; 

  gamma_E_               = -999;
  gamma_mass_            = -999;
  gamma_pt_              = -999; 
  gamma_eta_             = -999; 
  gamma_phi_             = -999; 
  gamma_id_              = -999; 
  gamma_from_p_          = -999; 
  gamma_px_              = -999; 
  gamma_py_              = -999; 
  gamma_pz_              = -999; 

  lepFromGamma_E_               = -999;
  lepFromGamma_mass_            = -999;
  lepFromGamma_pt_              = -999; 
  lepFromGamma_eta_             = -999; 
  lepFromGamma_phi_             = -999; 
  lepFromGamma_id_              = -999; 

  pIn_E_               = -999;
  pIn_mass_            = -999;
  pIn_px_              = -999; 
  pIn_py_             = -999; 
  pIn_pz_             = -999; 
  pIn_id_              = -999; 

  pOut_E_               = -999;
  pOut_mass_            = -999;
  pOut_px_              = -999; 
  pOut_py_             = -999; 
  pOut_pz_             = -999; 
  pOut_id_              = -999; 

  jet1_mass_            = -999;
  jet1_pt_              = -999; 
  jet1_eta_             = -999; 
  jet1_phi_             = -999; 

  xi_p_                 = -999; 
  xi_LQ_                = -999;  
  xi_LQ_wrong_          = -999;         
  xi_LQ_plusLep_                = -999;  
  xi_LQ_wrong_plusLep_          = -999;         

  lepIn_E_               = -999;
  lepIn_mass_            = -999;
  lepIn_pt_              = -999; 
  lepIn_eta_             = -999; 
  lepIn_phi_             = -999; 
  lepIn_id_              = -999; 
  lepIn_px_              = -999; 
  lepIn_py_              = -999; 
  lepIn_pz_              = -999; 

  qIn_E_               = -999;
  qIn_mass_            = -999;
  qIn_pt_              = -999; 
  qIn_eta_             = -999; 
  qIn_phi_             = -999; 
  qIn_id_              = -999; 
  qIn_px_              = -999; 
  qIn_py_              = -999; 
  qIn_pz_              = -999; 
  qIn_pz_              = -999; 

  lepOut_E_               = -999;
  lepOut_mass_            = -999;
  lepOut_pt_              = -999; 
  lepOut_eta_             = -999; 
  lepOut_phi_             = -999; 
  lepOut_id_              = -999; 

  qOut_E_               = -999;
  qOut_mass_            = -999;
  qOut_pt_              = -999; 
  qOut_eta_             = -999; 
  qOut_phi_             = -999; 
  qOut_id_              = -999; 

  lqIn_mass_              = -999; 
  lqOut_mass_              = -999; 

}

//define this as a plug-in
DEFINE_FWK_MODULE(GenAnalq);
