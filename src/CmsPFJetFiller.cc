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

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFJetFiller.h"

#include <string>

using namespace edm;
using namespace reco;


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
  delete privateData_->chargedEmEnergy;
  delete privateData_->neutralEmEnergy;
  delete privateData_->neutralMultiplicity;
  delete privateData_->chargedMultiplicity;
  delete privateData_->muonMultiplicity;
  delete privateData_->uncorrEnergy;

}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out

void CmsPFJetFiller::writeCollectionToTree(edm::InputTag collectionTag,
                                           const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                           const std::string &columnPrefix, const std::string &columnSuffix,
                                           bool dumpData,
                                           edm::InputTag uncorrectedCollectionTag) {

  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFJetFiller") << "Can't get candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();

  const edm::View<reco::Candidate> *uncorrectedCollection = 0;
  if(uncorrectedCollectionTag.label() != std::string("")) {
    edm::Handle< edm::View<reco::Candidate> > uncorrectedCollectionHandle;
    try { iEvent.getByLabel(uncorrectedCollectionTag, uncorrectedCollectionHandle); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFJetFiller") << "Can't get candidate collection: " << collectionTag; }
    uncorrectedCollection = uncorrectedCollectionHandle.product();
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
        privateData_->chargedEmEnergy->push_back( thisPFJet->chargedEmEnergy() );
        privateData_->neutralEmEnergy->push_back( thisPFJet->neutralEmEnergy() );
        privateData_->neutralMultiplicity->push_back( thisPFJet->neutralMultiplicity() );
        privateData_->chargedMultiplicity->push_back( thisPFJet->chargedMultiplicity() );
        privateData_->muonMultiplicity->push_back( thisPFJet->muonMultiplicity() );
      }
      else {
        privateData_->chargedHadronEnergy->push_back( -1. );
        privateData_->neutralHadronEnergy->push_back( -1. );
        privateData_->chargedEmEnergy->push_back( -1. );
        privateData_->neutralEmEnergy->push_back( -1. );
        privateData_->neutralMultiplicity->push_back( -1. );
        privateData_->chargedMultiplicity->push_back( -1. );
        privateData_->muonMultiplicity->push_back( -1. );
      }

      // if an uncorrected jet collection is provided, save also the uncorrected energy
      if(uncorrectedCollection) {
        dumpUncorrEnergy_ = true;
        float rawEnergy = -1.;
        edm::View<reco::Candidate>::const_iterator cand2;
        for(cand2=uncorrectedCollection->begin(); cand2!=uncorrectedCollection->end(); cand2++) {
          const PFJet *uncorrectedPFJet = dynamic_cast< const PFJet * > ( &(*cand2) );
          // corrected and uncorrected jets differ only for jet PT 
          if( thisPFJet->neutralEmEnergy() == uncorrectedPFJet->neutralEmEnergy() ) {
            rawEnergy = uncorrectedPFJet->energy();
            break;
          }
        }
        privateData_->uncorrEnergy->push_back(rawEnergy);
      } else dumpUncorrEnergy_ = false;

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
  cmstree->column((colPrefix+"chargedEmEnergy"+colSuffix).c_str(), *privateData_->chargedEmEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"neutralEmEnergy"+colSuffix).c_str(), *privateData_->neutralEmEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"neutralMultiplicity"+colSuffix).c_str(), *privateData_->neutralMultiplicity, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chargedMultiplicity"+colSuffix).c_str(), *privateData_->chargedMultiplicity, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"muonMultiplicity"+colSuffix).c_str(), *privateData_->muonMultiplicity, nCandString.c_str(), 0, "Reco");
  if(dumpUncorrEnergy_) {
    cmstree->column((colPrefix+"uncorrEnergy"+colSuffix).c_str(), *privateData_->uncorrEnergy, nCandString.c_str(), 0, "Reco");
  }

}







void CmsPFJetFillerData::initialise() {

  initialiseCandidate();
  chargedHadronEnergy = new vector<float>;
  neutralHadronEnergy = new vector<float>;
  chargedEmEnergy = new vector<float>;
  neutralEmEnergy = new vector<float>;
  neutralMultiplicity = new vector<float>;
  chargedMultiplicity = new vector<float>;
  muonMultiplicity = new vector<float>;
  uncorrEnergy = new vector<float>;

}

void CmsPFJetFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();
  chargedHadronEnergy->clear();
  neutralHadronEnergy->clear();
  chargedEmEnergy->clear();
  neutralEmEnergy->clear();
  neutralMultiplicity->clear();
  chargedMultiplicity->clear();
  muonMultiplicity->clear();
  uncorrEnergy->clear();

}
