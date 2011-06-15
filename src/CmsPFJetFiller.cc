//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsPFJetFiller
//
// Original Author:  Emanuele Di Marco
//         Created:  Mon Sep  29 11:01:00 CEST 2008
//-----------------------------------------------------------------------

// system include files
#include <memory>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFJetFiller.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include <string>

using namespace edm;
using namespace reco;



QGLikelihoodVars computeQGLikelihoodVars( const PFJet* pfjet, float R=0., float ptratio=0.);


//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------


CmsPFJetFiller::CmsPFJetFiller(CmsTree *cmsTree, 
			   int maxTracks, int maxMCTracks,
			   bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsPFJetFillerData)
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

CmsPFJetFiller::~CmsPFJetFiller() {

  // delete here the vector ptr's
  delete privateData_->chargedHadronEnergy;
  delete privateData_->neutralHadronEnergy;
  delete privateData_->photonEnergy;
  delete privateData_->electronEnergy;
  delete privateData_->muonEnergy;
  delete privateData_->HFHadronEnergy;
  delete privateData_->HFEMEnergy;
  delete privateData_->chargedHadronMultiplicity;
  delete privateData_->neutralHadronMultiplicity;
  delete privateData_->photonMultiplicity;
  delete privateData_->electronMultiplicity;
  delete privateData_->muonMultiplicity;
  delete privateData_->HFHadronMultiplicity;
  delete privateData_->HFEMMultiplicity;
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
  delete privateData_->L2L3CorrEnergy;
  delete privateData_->area;

  // for backward compatibility with existing trees
  delete privateData_->chargedEmEnergy;
  delete privateData_->neutralEmEnergy;

  // additional variables for Marini's likelihood:
  delete privateData_->ptD;
  delete privateData_->rmsCand;
}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out

void CmsPFJetFiller::saveJetBTag(bool what) { saveJetBTag_=what; }

void CmsPFJetFiller::setBTags(edm::ParameterSet btagcollections) { BTagCollections_ = btagcollections; }

// Set boolean control options for quantities that are written out

void CmsPFJetFiller::writeCollectionToTree(edm::InputTag collectionTag,
                                           const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                           const std::string &columnPrefix, const std::string &columnSuffix,
                                           bool dumpData,
                                           edm::InputTag uncorrectedCollectionTag, edm::InputTag L2L3correctedCollectionTag) {

  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFJetFiller") << "Can't get candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();

  const edm::View<reco::Candidate> *uncorrectedCollection = 0;
  if(uncorrectedCollectionTag.label() != std::string("")) {
    edm::Handle< edm::View<reco::Candidate> > uncorrectedCollectionHandle;
    try { iEvent.getByLabel(uncorrectedCollectionTag, uncorrectedCollectionHandle); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFJetFiller") << "Can't get candidate collection: " << uncorrectedCollectionTag; }
    uncorrectedCollection = uncorrectedCollectionHandle.product();
  }

  const edm::View<reco::Candidate> *L2L3correctedCollection = 0;
  if(L2L3correctedCollectionTag.label() != std::string("")) {
    edm::Handle< edm::View<reco::Candidate> > L2L3correctedCollectionHandle;
    try { iEvent.getByLabel(L2L3correctedCollectionTag, L2L3correctedCollectionHandle); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFJetFiller") << "Can't get candidate collection: " << L2L3correctedCollectionTag; }
    L2L3correctedCollection = L2L3correctedCollectionHandle.product();
  }

  privateData_->clearTrkVectors();

  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      edm::LogInfo("CmsPFJetFiller") << "Track length " << collection->size() 
				     << " is too long for declared max length for tree "
				     << maxTracks_ << " and no output flag is set."
				     << " No tracks written to tuple for this event ";
      return;
    }
    
    if((int)collection->size() > maxTracks_){
      edm::LogInfo("CmsPFJetFiller") << "Track length " << collection->size() 
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
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("combinedSecondaryVertexBJetTags"), combinedSecondaryVertexBJetTags);
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("combinedSecondaryVertexMVABJetTags"), combinedSecondaryVertexMVABJetTags);
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("jetBProbabilityBJetTags"), jetBProbabilityBJetTags);
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("jetProbabilityBJetTags"), jetProbabilityBJetTags);
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("simpleSecondaryVertexHighEffBJetTags"), simpleSecondaryVertexHighEffBJetTags);
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("simpleSecondaryVertexHighPurBJetTags"), simpleSecondaryVertexHighPurBJetTags);
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("softMuonBJetTags"), softMuonBJetTags);
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("softMuonByIP3dBJetTags"), softMuonByIP3dBJetTags);
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("softMuonByPtBJetTags"), softMuonByPtBJetTags);
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("softElectronBJetTags"), softElectronBJetTags);
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("softElectronByIP3dBJetTags"), softElectronByIP3dBJetTags);
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("softElectronByPtBJetTags"), softElectronByPtBJetTags);
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("trackCountingHighPurBJetTags"), trackCountingHighPurBJetTags);
      iEvent.getByLabel(BTagCollections_.getParameter<edm::InputTag>("trackCountingHighEffBJetTags"), trackCountingHighEffBJetTags);
    }

    int index = 0;
    edm::View<reco::Candidate>::const_iterator cand;
    for(cand=collection->begin(); cand!=collection->end(); cand++) {
      // fill basic kinematics
      if(saveCand_) writeCandInfo(&(*cand),iEvent,iSetup);
      // fill jet extra informations

      // em, had fractions
      const PFJet *thisPFJet = dynamic_cast< const PFJet * > ( &(*cand) );
      if( thisPFJet != 0 ) { 
        privateData_->chargedHadronEnergy->push_back( thisPFJet->chargedHadronEnergy() );
        privateData_->neutralHadronEnergy->push_back( thisPFJet->neutralHadronEnergy() );
        privateData_->photonEnergy->push_back( thisPFJet->photonEnergy() );
        privateData_->electronEnergy->push_back( thisPFJet->electronEnergy() );
        privateData_->muonEnergy->push_back( thisPFJet->muonEnergy() );
        privateData_->HFHadronEnergy->push_back( thisPFJet->HFHadronEnergy() );
        privateData_->HFEMEnergy->push_back( thisPFJet->HFEMEnergy() );
        privateData_->chargedHadronMultiplicity->push_back( thisPFJet->chargedHadronMultiplicity() );
        privateData_->neutralHadronMultiplicity->push_back( thisPFJet->neutralHadronMultiplicity() );
        privateData_->photonMultiplicity->push_back( thisPFJet->photonMultiplicity() );
        privateData_->electronMultiplicity->push_back( thisPFJet->electronMultiplicity() );
        privateData_->muonMultiplicity->push_back( thisPFJet->muonMultiplicity() );
        privateData_->HFHadronMultiplicity->push_back( thisPFJet->HFHadronMultiplicity() );
        privateData_->HFEMMultiplicity->push_back( thisPFJet->HFEMMultiplicity() );
        
        // for backward compatibility with existing trees
        privateData_->chargedEmEnergy->push_back( thisPFJet->chargedEmEnergy() );
        privateData_->neutralEmEnergy->push_back( thisPFJet->neutralEmEnergy() );
        privateData_->area->push_back( thisPFJet->jetArea() );

        // compute marini's variables
        QGLikelihoodVars qgvars = computeQGLikelihoodVars(thisPFJet);
        privateData_->ptD->push_back( qgvars.ptD );
        privateData_->rmsCand->push_back( qgvars.rmsCand );

      }
      else {
        privateData_->chargedHadronEnergy->push_back( -1. );
        privateData_->neutralHadronEnergy->push_back( -1. );
        privateData_->photonEnergy->push_back( -1. );
        privateData_->electronEnergy->push_back( -1. );
        privateData_->muonEnergy->push_back( -1. );
        privateData_->chargedHadronMultiplicity->push_back( -1. );
        privateData_->neutralHadronMultiplicity->push_back( -1. );
        privateData_->photonMultiplicity->push_back( -1. );
        privateData_->electronMultiplicity->push_back( -1. );
        privateData_->muonMultiplicity->push_back( -1. );

        // for backward compatibility with existing trees
        privateData_->chargedEmEnergy->push_back( -1. );
        privateData_->neutralEmEnergy->push_back( -1. );
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
          const PFJet *uncorrectedPFJet = dynamic_cast< const PFJet * > ( &(*cand2) );
          // corrected and uncorrected jets differ only for jet PT 
          if(  thisPFJet->jetArea() == uncorrectedPFJet->jetArea() ) {
            rawEnergy = uncorrectedPFJet->energy();
            break;
          }
        }
        privateData_->uncorrEnergy->push_back(rawEnergy);
      } else dumpUncorrEnergy_ = false;

      // if an L2L3corrected jet collection is provided, save also the L2L3corrected energy
      if(L2L3correctedCollection) {
        dumpL2L3CorrEnergy_ = true;
        float L2L3Energy = -999.;
        edm::View<reco::Candidate>::const_iterator cand2;
        for(cand2=L2L3correctedCollection->begin(); cand2!=L2L3correctedCollection->end(); cand2++) {
          const PFJet *L2L3correctedPFJet = dynamic_cast< const PFJet * > ( &(*cand2) );
          // corrected and uncorrected jets differ only for jet PT 
          if(  thisPFJet->jetArea() == L2L3correctedPFJet->jetArea() ) {
            L2L3Energy = L2L3correctedPFJet->energy();
            break;
          }
        }
        privateData_->L2L3CorrEnergy->push_back(L2L3Energy);
      } else dumpL2L3CorrEnergy_ = false;
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









void CmsPFJetFiller::treeJetInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"chargedHadronEnergy"+colSuffix).c_str(), *privateData_->chargedHadronEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"neutralHadronEnergy"+colSuffix).c_str(), *privateData_->neutralHadronEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"photonEnergy"+colSuffix).c_str(), *privateData_->photonEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"electronEnergy"+colSuffix).c_str(), *privateData_->electronEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"muonEnergy"+colSuffix).c_str(), *privateData_->muonEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"HFHadronEnergy"+colSuffix).c_str(), *privateData_->HFHadronEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"HFEMEnergy"+colSuffix).c_str(), *privateData_->HFEMEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chargedHadronMultiplicity"+colSuffix).c_str(), *privateData_->chargedHadronMultiplicity, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"neutralHadronMultiplicity"+colSuffix).c_str(), *privateData_->neutralHadronMultiplicity, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"photonMultiplicity"+colSuffix).c_str(), *privateData_->photonMultiplicity, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"electronMultiplicity"+colSuffix).c_str(), *privateData_->electronMultiplicity, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"muonMultiplicity"+colSuffix).c_str(), *privateData_->muonMultiplicity, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"HFHadronMultiplicity"+colSuffix).c_str(), *privateData_->HFHadronMultiplicity, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"HFEMMultiplicity"+colSuffix).c_str(), *privateData_->HFEMMultiplicity, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"area"+colSuffix).c_str(), *privateData_->area, nCandString.c_str(), 0, "Reco");

  // for backward compatibility with existing trees 
  cmstree->column((colPrefix+"chargedEmEnergy"+colSuffix).c_str(), *privateData_->chargedEmEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"neutralEmEnergy"+colSuffix).c_str(), *privateData_->neutralEmEnergy, nCandString.c_str(), 0, "Reco");
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
  if(dumpL2L3CorrEnergy_) {
    cmstree->column((colPrefix+"L2L3CorrEnergy"+colSuffix).c_str(), *privateData_->L2L3CorrEnergy, nCandString.c_str(), 0, "Reco");
  }
  cmstree->column((colPrefix+"rmsCand"+colSuffix).c_str(), *privateData_->rmsCand, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"ptD"+colSuffix).c_str(), *privateData_->ptD, nCandString.c_str(), 0, "Reco");

}







void CmsPFJetFillerData::initialise() {

  initialiseCandidate();
  chargedHadronEnergy = new vector<float>;
  neutralHadronEnergy = new vector<float>;
  chargedEmEnergy = new vector<float>;
  neutralEmEnergy = new vector<float>;
  photonEnergy = new vector<float>;
  electronEnergy = new vector<float>;
  muonEnergy = new vector<float>;
  HFHadronEnergy = new vector<float>;
  HFEMEnergy = new vector<float>;
  chargedHadronMultiplicity = new vector<int>;
  neutralHadronMultiplicity = new vector<int>;
  photonMultiplicity = new vector<int>;
  electronMultiplicity = new vector<int>;
  muonMultiplicity = new vector<int>;
  HFHadronMultiplicity = new vector<int>;
  HFEMMultiplicity = new vector<int>;
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
  L2L3CorrEnergy = new vector<float>;
  area = new vector<float>;
  ptD = new vector<float>;
  rmsCand = new vector<float>;

}

void CmsPFJetFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();
  chargedHadronEnergy->clear();
  neutralHadronEnergy->clear();
  chargedEmEnergy->clear();
  neutralEmEnergy->clear();
  photonEnergy->clear();
  electronEnergy->clear();
  muonEnergy->clear();
  HFHadronEnergy->clear();
  HFEMEnergy->clear();
  chargedHadronMultiplicity->clear();
  neutralHadronMultiplicity->clear();
  photonMultiplicity->clear();
  electronMultiplicity->clear();
  muonMultiplicity->clear();
  HFHadronMultiplicity->clear();
  HFEMMultiplicity->clear();
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
  L2L3CorrEnergy->clear();
  area->clear();
  ptD->clear();
  rmsCand->clear();

}



QGLikelihoodVars computeQGLikelihoodVars( const PFJet* pfjet, float R, float ptratio ) {

  std::vector<fastjet::PseudoJet> *input_particles = new std::vector<fastjet::PseudoJet>;

  for(long int i=0;i<pfjet->nConstituents();++i)
     {
     input_particles->push_back(fastjet::PseudoJet( pfjet->getJetConstituentsQuick()[i]->px(),
                          pfjet->getJetConstituentsQuick()[i]->py(),
                          pfjet->getJetConstituentsQuick()[i]->pz(),
                          pfjet->getJetConstituentsQuick()[i]->energy()
                      ) );
     }

  fastjet::JetDefinition ak_def(fastjet::antikt_algorithm, R);
  fastjet::JetDefinition ak_def1(fastjet::antikt_algorithm, 0.04);
  fastjet::JetDefinition ak_def2(fastjet::antikt_algorithm, 0.25);
  
  fastjet::ClusterSequence seq(*input_particles, ak_def);
  fastjet::ClusterSequence seq1(*input_particles,ak_def1);
  fastjet::ClusterSequence seq2(*input_particles,ak_def2);
  
  float jtpt = pfjet->pt();

  //recompute ptmin
  float ptmin = ptratio * jtpt;
  
  //now I want to know how many subjets there are with pt>ptmin
  vector<fastjet::PseudoJet> inclusive_jets;
  
  inclusive_jets = sorted_by_pt(seq1.inclusive_jets(jtpt*0.30));
  //int nsubjets1=inclusive_jets.size();
  
  inclusive_jets = sorted_by_pt(seq2.inclusive_jets(jtpt*0.05));
  //  int nsubjets2=inclusive_jets.size();
  
  inclusive_jets = sorted_by_pt(seq.inclusive_jets(ptmin));
  //  int nsubjets =  inclusive_jets.size(); 
  
  //cerco il raggio del jet
  //the next definition of WEIGHT is used in the radius: Pt^2
  Double_t jtpt_s=0;//sum to get the jtpt that I see after clustering
  if(inclusive_jets.size()>0)
    {
    jtpt_s=0;
    for(unsigned int i = 0; i<inclusive_jets.size(); i++ )
    		jtpt_s+=inclusive_jets.at(i).perp();
  }
  #define WEIGHT (inclusive_jets.at(i).perp()*inclusive_jets.at(i).perp())
  float r=0.;
  if(inclusive_jets.size()>0)
  {
  Double_t eta0=0,phi0=0,Sum=0.0;
  for(unsigned int i=0;i<inclusive_jets.size();i++)
    {
    Sum+= WEIGHT ;
    eta0+= WEIGHT * inclusive_jets.at(i).eta();
    //problema con phi
    if(1.5<inclusive_jets.at(0).phi()&&inclusive_jets.at(0).phi()<4.6)
    	phi0+= WEIGHT * inclusive_jets.at(i).phi();
    else
    	phi0+= WEIGHT * inclusive_jets.at(i).phi_std();
    }
  eta0/=Sum;
  phi0/=Sum;
  r=0.;
  for(unsigned int i=0;i<inclusive_jets.size();i++)
    {
    if(1.5<inclusive_jets.at(0).phi()&&inclusive_jets.at(0).phi()<4.6)
    r+= ((inclusive_jets.at(i).eta()-eta0)*(inclusive_jets.at(i).eta()-eta0) 
    		+ (inclusive_jets.at(i).phi()-phi0)*(inclusive_jets.at(i).phi()-phi0)) * WEIGHT ;
    else
    r+=((inclusive_jets.at(i).eta()-eta0)*(inclusive_jets.at(i).eta()-eta0) 
    		+ (inclusive_jets.at(i).phi_std()-phi0)*(inclusive_jets.at(i).phi_std()-phi0))*WEIGHT;
    }
    r/=Sum;
  }else r=0;
  //I want to compute the PtD =\sqrt( \sum_j (Pt_j/Pt_jet) ^ 2)
  float PtD=0.0;
  for(unsigned int i=0 ; i < inclusive_jets.size(); i++ )
    {
    PtD+=inclusive_jets.at(i).perp2()/(jtpt_s*jtpt_s);
    }
  PtD=sqrt(PtD);
  
//  //I compute sub1ptratio
//if (inclusive_jets.size()>0)
//	sub1ptratio=inclusive_jets.at(0).perp()/jtpt_s;
//else 
//	sub1ptratio=0;

  QGLikelihoodVars vars;
  vars.ptD = PtD;
  vars.rmsCand = r;

  return vars;

}
