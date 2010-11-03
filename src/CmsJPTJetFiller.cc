//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsJPTJetFiller
//
// Original Author:  Emanuele Di Marco
//         Created:  Mon Sep  29 11:01:00 CEST 2008
//-----------------------------------------------------------------------

// system include files
#include <memory>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JPTJet.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsJPTJetFiller.h"

#include <string>

using namespace edm;
using namespace reco;


//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------


CmsJPTJetFiller::CmsJPTJetFiller(CmsTree *cmsTree, 
                                 int maxTracks, int maxMCTracks,
                                 bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsJPTJetFillerData)
{
  cmstree=cmsTree;

  saveJetBTag_ = false;

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();
}


//--------------
// Destructor --
//--------------

CmsJPTJetFiller::~CmsJPTJetFiller() {

  // delete here the vector ptr's
  delete privateData_->chargedHadronEnergy;
  delete privateData_->neutralHadronEnergy;
  delete privateData_->chargedEmEnergy;
  delete privateData_->neutralEmEnergy;
  delete privateData_->chargedMultiplicity;
  delete privateData_->muonMultiplicity;
  delete privateData_->elecMultiplicity;
  delete privateData_->ZSPCor;
  delete privateData_->Eta2momtr;
  delete privateData_->Phi2momtr;
  delete privateData_->combinedSecondaryVertexBJetTags;
  delete privateData_->combinedSecondaryVertexMVABJetTags;
  delete privateData_->jetBProbabilityBJetTags;
  delete privateData_->jetProbabilityBJetTags;
  delete privateData_->simpleSecondaryVertexHighEffBJetTags;
  delete privateData_->simpleSecondaryVertexHighPurBJetTags;
  delete privateData_->softMuonBJetTags;
  delete privateData_->softMuonByIP3dBJetTags;
  delete privateData_->softMuonByPtBJetTags;
  delete privateData_->softElectronBJetTags;
  delete privateData_->softElectronByIP3dBJetTags;
  delete privateData_->softElectronByPtBJetTags;
  delete privateData_->trackCountingHighPurBJetTags;
  delete privateData_->trackCountingHighEffBJetTags;
  delete privateData_->uncorrEnergy;

}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out

void CmsJPTJetFiller::saveJetBTag(bool what) { saveJetBTag_=what; }

// Set boolean control options for quantities that are written out

void CmsJPTJetFiller::writeCollectionToTree(edm::InputTag collectionTag,
                                            const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                            const std::string &columnPrefix, const std::string &columnSuffix,
                                            bool dumpData,
                                            edm::InputTag uncorrectedCollectionTag) {

  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsJPTJetFiller") << "Can't get candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();

  const edm::View<reco::Candidate> *uncorrectedCollection = 0;
  if(uncorrectedCollectionTag.label() != std::string("")) {
    edm::Handle< edm::View<reco::Candidate> > uncorrectedCollectionHandle;
    try { iEvent.getByLabel(uncorrectedCollectionTag, uncorrectedCollectionHandle); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsJPTJetFiller") << "Can't get candidate collection: " << collectionTag; }
    uncorrectedCollection = uncorrectedCollectionHandle.product();
  }

  privateData_->clearTrkVectors();

  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      edm::LogInfo("CmsJPTJetFiller") << "Track length " << collection->size() 
                                      << " is too long for declared max length for tree "
                                      << maxTracks_ << " and no output flag is set."
                                      << " No tracks written to tuple for this event ";
      return;
    }
    
    if((int)collection->size() > maxTracks_){
      edm::LogInfo("CmsJPTJetFiller") << "Track length " << collection->size() 
                                      << " is too long for declared max length for tree "
                                      << maxTracks_ 
                                      << ". Collection will be truncated ";
    }
  
    *(privateData_->ncand) = collection->size();

    edm::Handle<reco::JetTagCollection> combinedSecondaryVertexBJetTags,
      combinedSecondaryVertexMVABJetTags,
      jetBProbabilityBJetTags,
      jetProbabilityBJetTags,
      simpleSecondaryVertexHighEffBJetTags,
      simpleSecondaryVertexHighPurBJetTags,
      softMuonBJetTags,
      softMuonByIP3dBJetTags,
      softMuonByPtBJetTags,
      softElectronBJetTags,
      softElectronByIP3dBJetTags,
      softElectronByPtBJetTags,
      trackCountingHighPurBJetTags,
      trackCountingHighEffBJetTags;

    if(saveJetBTag_) {
      iEvent.getByLabel("newCombinedSecondaryVertexBJPTJetTags", combinedSecondaryVertexBJetTags);
      iEvent.getByLabel("newCombinedSecondaryVertexMVABJPTJetTags", combinedSecondaryVertexMVABJetTags);
      iEvent.getByLabel("newJetBProbabilityBJPTJetTags", jetBProbabilityBJetTags);
      iEvent.getByLabel("newJetProbabilityBJPTJetTags", jetProbabilityBJetTags);
      iEvent.getByLabel("newSimpleSecondaryVertexHighEffBJPTJetTags", simpleSecondaryVertexHighEffBJetTags);
      iEvent.getByLabel("newSimpleSecondaryVertexHighPurBJPTJetTags", simpleSecondaryVertexHighPurBJetTags);
      iEvent.getByLabel("newSoftMuonBJPTJetTags", softMuonBJetTags);
      iEvent.getByLabel("newSoftMuonByIP3dBJPTJetTags", softMuonByIP3dBJetTags);
      iEvent.getByLabel("newSoftMuonByPtBJPTJetTags", softMuonByPtBJetTags);
      iEvent.getByLabel("newSoftElectronBJPTJetTags", softElectronBJetTags);
      iEvent.getByLabel("newSoftElectronByIP3dBJPTJetTags", softElectronByIP3dBJetTags);
      iEvent.getByLabel("newSoftElectronByPtBJPTJetTags", softElectronByPtBJetTags);
      iEvent.getByLabel("newTrackCountingHighPurBJPTJetTags", trackCountingHighPurBJetTags);
      iEvent.getByLabel("newTrackCountingHighEffBJPTJetTags", trackCountingHighEffBJetTags);
    }

    int index = 0;
    edm::View<reco::Candidate>::const_iterator cand;
    for(cand=collection->begin(); cand!=collection->end(); cand++) {
      // fill basic kinematics
      if(saveCand_) writeCandInfo(&(*cand),iEvent,iSetup);
      // fill jet extra informations

      // em, had fractions
      const JPTJet *thisJPTJet = dynamic_cast< const JPTJet * > ( &(*cand) );
      if( thisJPTJet != 0 ) { 
        privateData_->chargedHadronEnergy->push_back( thisJPTJet->chargedHadronEnergy() );
        privateData_->neutralHadronEnergy->push_back( thisJPTJet->neutralHadronEnergy() );
        privateData_->chargedEmEnergy->push_back( thisJPTJet->chargedEmEnergy() );
        privateData_->neutralEmEnergy->push_back( thisJPTJet->neutralEmEnergy() );
        privateData_->chargedMultiplicity->push_back( thisJPTJet->chargedMultiplicity() );
        privateData_->muonMultiplicity->push_back( thisJPTJet->muonMultiplicity() );
        privateData_->elecMultiplicity->push_back( thisJPTJet->elecMultiplicity() );
        privateData_->ZSPCor->push_back( thisJPTJet->getZSPCor() );

        JPTJet::Specific spec = thisJPTJet->getSpecific();
        privateData_->Eta2momtr->push_back( spec.Eta2momtr );
        privateData_->Phi2momtr->push_back( spec.Phi2momtr );
      }
      else {
        privateData_->chargedHadronEnergy->push_back( -1. );
        privateData_->neutralHadronEnergy->push_back( -1. );
        privateData_->chargedEmEnergy->push_back( -1. );
        privateData_->neutralEmEnergy->push_back( -1. );
        privateData_->chargedMultiplicity->push_back( -1. );
        privateData_->muonMultiplicity->push_back( -1. );
        privateData_->elecMultiplicity->push_back( -1. );
        privateData_->ZSPCor->push_back( -1. );
        privateData_->Eta2momtr->push_back( -1. );
        privateData_->Phi2momtr->push_back( -1. );
      }

      // fill the btag algorithms output
      if(saveJetBTag_) {
        privateData_->combinedSecondaryVertexBJetTags->push_back( (*combinedSecondaryVertexBJetTags)[index].second );
        privateData_->combinedSecondaryVertexMVABJetTags->push_back( (*combinedSecondaryVertexMVABJetTags)[index].second );
        privateData_->jetBProbabilityBJetTags->push_back( (*jetBProbabilityBJetTags)[index].second );
        privateData_->jetProbabilityBJetTags->push_back( (*jetProbabilityBJetTags)[index].second );
        privateData_->simpleSecondaryVertexHighEffBJetTags->push_back( (*simpleSecondaryVertexHighEffBJetTags)[index].second );
        privateData_->simpleSecondaryVertexHighPurBJetTags->push_back( (*simpleSecondaryVertexHighPurBJetTags)[index].second );
        privateData_->softMuonBJetTags->push_back( (*softMuonBJetTags)[index].second );
        privateData_->softMuonByIP3dBJetTags->push_back( (*softMuonByIP3dBJetTags)[index].second );
        privateData_->softMuonByPtBJetTags->push_back( (*softMuonByPtBJetTags)[index].second );
        privateData_->softElectronBJetTags->push_back( (*softElectronBJetTags)[index].second );
        privateData_->softElectronByIP3dBJetTags->push_back( (*softElectronByIP3dBJetTags)[index].second );
        privateData_->softElectronByPtBJetTags->push_back( (*softElectronByPtBJetTags)[index].second );
        privateData_->trackCountingHighPurBJetTags->push_back( (*trackCountingHighPurBJetTags)[index].second );
        privateData_->trackCountingHighEffBJetTags->push_back( (*trackCountingHighEffBJetTags)[index].second );
      } else {
        privateData_->combinedSecondaryVertexBJetTags->push_back( -1. );
        privateData_->combinedSecondaryVertexMVABJetTags->push_back( -1. );
        privateData_->jetBProbabilityBJetTags->push_back( -1. );
        privateData_->jetProbabilityBJetTags->push_back( -1. );
        privateData_->simpleSecondaryVertexHighEffBJetTags->push_back( -1. );
        privateData_->simpleSecondaryVertexHighPurBJetTags->push_back( -1. );
        privateData_->softMuonBJetTags->push_back( 1. );
        privateData_->softMuonByIP3dBJetTags->push_back( -1. );
        privateData_->softMuonByPtBJetTags->push_back( -1. );
        privateData_->softElectronBJetTags->push_back( -1. );
        privateData_->softElectronByIP3dBJetTags->push_back( -1. );
        privateData_->softElectronByPtBJetTags->push_back( -1. );
        privateData_->trackCountingHighPurBJetTags->push_back( -1. );
        privateData_->trackCountingHighEffBJetTags->push_back( -1. );
      }

      // if an uncorrected jet collection is provided, save also the uncorrected energy
      if(uncorrectedCollection) {
        dumpUncorrEnergy_ = true;
        float rawEnergy = -1.;
        edm::View<reco::Candidate>::const_iterator cand2;
        for(cand2=uncorrectedCollection->begin(); cand2!=uncorrectedCollection->end(); cand2++) {
          const JPTJet *uncorrectedJPTJet = dynamic_cast< const JPTJet * > ( &(*cand2) );
          // corrected and uncorrected jets differ only for jet PT 
          if( thisJPTJet->eta() == uncorrectedJPTJet->eta() && thisJPTJet->phi() == uncorrectedJPTJet->phi() ) {
            rawEnergy = uncorrectedJPTJet->energy();
            break;
          }
        }
        privateData_->uncorrEnergy->push_back(rawEnergy);
      } else dumpUncorrEnergy_ = false;
      index++; 
    }
  }
  else {
    *(privateData_->ncand) = 0;
  }
  
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the 
  // tree 
  int blockSize = (collection) ? collection->size() : 0;
  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");
  
  if(saveCand_) treeCandInfo(columnPrefix,columnSuffix);
  treeJetInfo(columnPrefix,columnSuffix);

  if(dumpData) cmstree->dumpData();
	
}









void CmsJPTJetFiller::treeJetInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"chargedHadronEnergy"+colSuffix).c_str(), *privateData_->chargedHadronEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"neutralHadronEnergy"+colSuffix).c_str(), *privateData_->neutralHadronEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chargedEmEnergy"+colSuffix).c_str(), *privateData_->chargedEmEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"neutralEmEnergy"+colSuffix).c_str(), *privateData_->neutralEmEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chargedMultiplicity"+colSuffix).c_str(), *privateData_->chargedMultiplicity, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"muonMultiplicity"+colSuffix).c_str(), *privateData_->muonMultiplicity, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"elecMultiplicity"+colSuffix).c_str(), *privateData_->elecMultiplicity, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"ZSPCor"+colSuffix).c_str(), *privateData_->ZSPCor, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"Eta2momtr"+colSuffix).c_str(), *privateData_->Eta2momtr, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"Phi2momtr"+colSuffix).c_str(), *privateData_->Phi2momtr, nCandString.c_str(), 0, "Reco");
  if(saveJetBTag_) {
    cmstree->column((colPrefix+"combinedSecondaryVertexBJetTags"+colSuffix).c_str(), *privateData_->combinedSecondaryVertexBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"combinedSecondaryVertexMVABJetTags"+colSuffix).c_str(), *privateData_->combinedSecondaryVertexMVABJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"jetBProbabilityBJetTags"+colSuffix).c_str(), *privateData_->jetBProbabilityBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"jetProbabilityBJetTags"+colSuffix).c_str(), *privateData_->jetProbabilityBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"simpleSecondaryVertexHighEffBJetTags"+colSuffix).c_str(), *privateData_->simpleSecondaryVertexHighEffBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"simpleSecondaryVertexHighPurBJetTags"+colSuffix).c_str(), *privateData_->simpleSecondaryVertexHighPurBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"softMuonBJetTags"+colSuffix).c_str(), *privateData_->softMuonBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"softMuonByIP3dBJetTags"+colSuffix).c_str(), *privateData_->softMuonByIP3dBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"softMuonByPtBJetTags"+colSuffix).c_str(), *privateData_->softMuonByPtBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"softElectronBJetTags"+colSuffix).c_str(), *privateData_->softElectronBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"softElectronByIP3dBJetTags"+colSuffix).c_str(), *privateData_->softElectronByIP3dBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"softElectronByPtBJetTags"+colSuffix).c_str(), *privateData_->softElectronByPtBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackCountingHighPurBJetTags"+colSuffix).c_str(), *privateData_->trackCountingHighPurBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackCountingHighEffBJetTags"+colSuffix).c_str(), *privateData_->trackCountingHighEffBJetTags, nCandString.c_str(), 0, "Reco");
  }
  if(dumpUncorrEnergy_) {
    cmstree->column((colPrefix+"uncorrEnergy"+colSuffix).c_str(), *privateData_->uncorrEnergy, nCandString.c_str(), 0, "Reco");
  }

}







void CmsJPTJetFillerData::initialise() {

  initialiseCandidate();
  chargedHadronEnergy = new vector<float>;
  neutralHadronEnergy = new vector<float>;
  chargedEmEnergy = new vector<float>;
  neutralEmEnergy = new vector<float>;
  chargedMultiplicity = new vector<float>;
  muonMultiplicity = new vector<float>;
  elecMultiplicity = new vector<float>;
  ZSPCor = new vector<float>;
  Eta2momtr = new vector<float>;
  Phi2momtr = new vector<float>;
  combinedSecondaryVertexBJetTags = new vector<float>;
  combinedSecondaryVertexMVABJetTags = new vector<float>;
  jetBProbabilityBJetTags = new vector<float>;
  jetProbabilityBJetTags = new vector<float>;
  simpleSecondaryVertexHighEffBJetTags = new vector<float>;
  simpleSecondaryVertexHighPurBJetTags = new vector<float>;
  softMuonBJetTags = new vector<float>;
  softMuonByIP3dBJetTags = new vector<float>;
  softMuonByPtBJetTags = new vector<float>;
  softElectronBJetTags = new vector<float>;
  softElectronByIP3dBJetTags = new vector<float>;
  softElectronByPtBJetTags = new vector<float>;
  trackCountingHighPurBJetTags = new vector<float>;
  trackCountingHighEffBJetTags = new vector<float>;
  uncorrEnergy = new vector<float>;

}

void CmsJPTJetFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();
  chargedHadronEnergy->clear();
  neutralHadronEnergy->clear();
  chargedEmEnergy->clear();
  neutralEmEnergy->clear();
  chargedMultiplicity->clear();
  muonMultiplicity->clear();
  elecMultiplicity->clear();
  ZSPCor->clear();
  Eta2momtr->clear();
  Phi2momtr->clear();
  combinedSecondaryVertexBJetTags->clear();
  combinedSecondaryVertexMVABJetTags->clear();
  jetBProbabilityBJetTags->clear();
  jetProbabilityBJetTags->clear();
  simpleSecondaryVertexHighEffBJetTags->clear();
  simpleSecondaryVertexHighPurBJetTags->clear();
  softMuonBJetTags->clear();
  softMuonByIP3dBJetTags->clear();
  softMuonByPtBJetTags->clear();
  softElectronBJetTags->clear();
  softElectronByIP3dBJetTags->clear();
  softElectronByPtBJetTags->clear();
  trackCountingHighPurBJetTags->clear();
  trackCountingHighEffBJetTags->clear();
  uncorrEnergy->clear();

}
