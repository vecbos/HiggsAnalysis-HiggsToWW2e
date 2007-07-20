// -*- C++ -*-
/*

   Description: <one line class summary>

   Implementation:
   This module just pick up the MET object from the event and put 
   it again in the event
 
*/
//
// Original Author:  Emanuele Di Marco
//         Created:  Mon Apr  2 17:28:56 CEST 2007
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"

#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWMetProducer.h"

using namespace std;
using namespace edm;
using namespace reco;

// --- Constructors ---
HWWMetProducer::HWWMetProducer(const edm::ParameterSet& iConfig)  {
  produceGenMET_   = iConfig.getUntrackedParameter<bool>("produceGenMET", false);
  genMetsTag_      = iConfig.getParameter<edm::InputTag>("genMetsTag");
  recoMetsTag_     = iConfig.getParameter<edm::InputTag>("recoMetsTag");
  genMetsOut_      = iConfig.getParameter<std::string>("genMetsOut");
  recoMetsOut_     = iConfig.getParameter<std::string>("recoMetsOut");
  genMetCandsOut_  = iConfig.getParameter<std::string>("genMetCandsOut");
  recoMetCandsOut_ = iConfig.getParameter<std::string>("recoMetCandsOut");
  min_Pt_          = iConfig.getUntrackedParameter<double>("MinPt", 0.0);
  outputFile_      = iConfig.getUntrackedParameter<string>("outputFile", "");

  produces<CaloMETCollection>(recoMetsOut_).setBranchAlias(recoMetsOut_);
  produces<CandidateCollection>(recoMetCandsOut_).setBranchAlias(recoMetCandsOut_);
  if(produceGenMET_) {
    produces<GenMETCollection>(genMetsOut_).setBranchAlias(genMetsOut_);
    produces<CandidateCollection>(genMetCandsOut_).setBranchAlias(genMetCandsOut_);
  }
}

HWWMetProducer::HWWMetProducer()   {
  produces<CaloMETCollection>();
}

// --- default destructor ---
HWWMetProducer::~HWWMetProducer() {}

// ------------ method called to produce the data  ------------
// --- produce the original Calo(Gen)MET collection plus a correspondent CandidateCollection
// --- for now only reco candidates are put into event
// ------------------------------------------------------------
void HWWMetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)  {
  // Reco MET
  Handle<CaloMETCollection> recoMetsHandle;
  iEvent.getByLabel(recoMetsTag_, recoMetsHandle);
  const CaloMETCollection *recoMetColl = recoMetsHandle.product();
  const CaloMET recoMet = recoMetColl->front();
  auto_ptr<CaloMETCollection> recoMetOut( new CaloMETCollection(*recoMetColl) );
  iEvent.put(recoMetOut,recoMetsOut_);
    
  // now produce the candidate collection
  auto_ptr<CandidateCollection> recoMetCandColl(new CandidateCollection);
  recoMetCandColl->reserve(recoMetColl->size());
  for(int i=0; i<(int)recoMetColl->size(); i++) {
    const CaloMET *itCaloMET = (CaloMET*) &recoMetColl->at(i);
    LeafCandidate *metCandColl = new LeafCandidate(0, Particle::LorentzVector(itCaloMET->px(), itCaloMET->py(), itCaloMET->pz(), itCaloMET->energy())); 
    recoMetCandColl->push_back(metCandColl);
  }
  iEvent.put(recoMetCandColl,recoMetCandsOut_);

  // Fill the monitoring histograms
  if(outputFile_.size() != 0) {
    hRecoMET_phi->Fill( recoMet.phi() );
    hRecoMET_MET->Fill( recoMet.energy() * sin(recoMet.theta()) );
    hRecoMET_sig->Fill( recoMet.mEtSig() );
    LogDebug("HWWMetProducer") << "Reco MET phi = " << recoMet.phi();
    LogDebug("HWWMetProducer") << "Reco MET et = " << recoMet.energy() * sin(recoMet.theta());
  }
    

  // Generated MET
  if(produceGenMET_) {
    Handle<GenMETCollection> genMetsHandle;
    iEvent.getByLabel(genMetsTag_, genMetsHandle);
    const GenMETCollection *genMetColl = genMetsHandle.product();
    const GenMET genMet = genMetColl->front();
    auto_ptr<GenMETCollection> genMetOut( new GenMETCollection(*genMetColl) );
    iEvent.put(genMetOut,genMetsOut_);

    // now produce the candidate collection
    auto_ptr<CandidateCollection> genMetCandColl(new CandidateCollection);
    genMetCandColl->reserve(genMetColl->size());
    for(int i=0; i<(int)genMetColl->size(); i++) {
      const GenMET *itGenMET = (GenMET*) &genMetColl->at(i);
      LeafCandidate *metCandColl = new LeafCandidate(0, Particle::LorentzVector(itGenMET->px(), itGenMET->py(), itGenMET->pz(), itGenMET->energy())); 
      genMetCandColl->push_back(metCandColl);
    }
    iEvent.put(genMetCandColl,genMetCandsOut_);


    // Fill the monitoring histograms
    if(outputFile_.size() != 0) {
      hGenMET_phi->Fill( genMet.phi() );
      hGenMET_MET->Fill( genMet.energy() * sin(genMet.theta()) );
      hGenMET_sig->Fill( genMet.mEtSig() );
      LogDebug("HWWMetProducer") << "Gen MET phi = " << genMet.phi();
      LogDebug("HWWMetProducer") << "Gen MET et = " << genMet.energy() * sin(genMet.theta());
    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
HWWMetProducer::beginJob(const edm::EventSetup&)
{
  if(outputFile_.size() != 0) {
    // Make the output files
    LogInfo("HWWMetProducer") << "MET producer monitoring histograms will be saved to '"
				<< outputFile_.c_str() << "'";
    dataFile_ = new TFile(outputFile_.c_str(), "RECREATE");
    // Book the histograms
    BookHistos();
  }
  else LogInfo("HWWMetProducer") << "MET producer monitoring histograms will NOT be saved";
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HWWMetProducer::endJob() {
  if(outputFile_.size() != 0) {
    dataFile_->Write();
    dataFile_->Close();
  }
}

void 
HWWMetProducer::BookHistos() {
  hRecoMET_phi = new TH1F("hRecoMET_phi", "Reconstructed MET phi", 50,-M_PI,M_PI);
  hRecoMET_MET = new TH1F("hRecoMET_MET", "Reconstructed MET met", 50,0,100.);
  hRecoMET_sig = new TH1F("hRecoMET_sig", "Reconstructed MET significance", 50,0.,1.);
  hGenMET_phi = new TH1F("hGenMET_phi", "Generated MET phi", 50,-M_PI,M_PI); 
  hGenMET_MET = new TH1F("hGenMET_MET", "Generated MET met", 50,0,100.);
  hGenMET_sig = new TH1F("hGenMET_sig", "Generated MET significance", 50,0.,1.);
}

