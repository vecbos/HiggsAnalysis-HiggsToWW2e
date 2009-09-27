//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsJetFiller
//
// Original Author:  Emanuele Di Marco
//         Created:  Mon Apr  9 11:01:00 CEST 2007
//-----------------------------------------------------------------------

// system include files
#include <memory>

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

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsJetFiller.h"

#include <TTree.h>

#include <string>

using namespace edm;
using namespace reco;


//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------


CmsJetFiller::CmsJetFiller(CmsTree *cmsTree, 
			   edm::InputTag jetVertexAlphaCollection,
			   int maxTracks, int maxMCTracks,
			   bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsJetFillerData)
{
  cmstree=cmsTree;
  jetVertexAlphaCollection_=jetVertexAlphaCollection;

  saveJetExtras_=true;

  saveJetBTag_ = false;
  
  trkIndexName_ = new std::string("n");
  
  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();
}

CmsJetFiller::CmsJetFiller(CmsTree *cmsTree, 
			   edm::InputTag jetVertexAlphaCollection,
			   bool fatTree, 
			   int maxTracks, int maxMCTracks,
			   bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,fatTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsJetFillerData)
{
  cmstree=cmsTree;
  jetVertexAlphaCollection_=jetVertexAlphaCollection;

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();
}


//--------------
// Destructor --
//--------------

CmsJetFiller::~CmsJetFiller() {

  // delete here the vector ptr's
  delete privateData_->emFrac;
  delete privateData_->hadFrac;
  delete privateData_->alpha;
  delete privateData_->combinedSecondaryVertexBJetTags;
  delete privateData_->combinedSecondaryVertexMVABJetTags;
  delete privateData_->jetBProbabilityBJetTags;
  delete privateData_->jetProbabilityBJetTags;
  delete privateData_->simpleSecondaryVertexBJetTags;
  delete privateData_->softElectronBJetTags;
  delete privateData_->softMuonBJetTags;
  delete privateData_->softMuonNoIPBJetTags;
  delete privateData_->trackCountingHighPurBJetTags;
  delete privateData_->trackCountingHighEffBJetTags;

}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out

void CmsJetFiller::saveJetExtras(bool what) { saveJetExtras_=what; }

void CmsJetFiller::saveJetBTag(bool what) { saveJetBTag_=what; }

void CmsJetFiller::writeCollectionToTree(edm::InputTag collectionTag,
					 const edm::Event& iEvent, const edm::EventSetup& iSetup,
					 const std::string &columnPrefix, const std::string &columnSuffix,
					 bool dumpData) {

  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsJetFiller") << "Can't get candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();
  
  privateData_->clearTrkVectors();

  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      edm::LogInfo("CmsJetFiller") << "Track length " << collection->size() 
			       << " is too long for declared max length for tree "
			       << maxTracks_ << " and no output flag is set."
			       << " No tracks written to tuple for this event ";
      return;
    }

    if((int)collection->size() > maxTracks_){
      edm::LogInfo("CmsJetFiller") << "Track length " << collection->size() 
			       << " is too long for declared max length for tree "
			       << maxTracks_ 
			       << ". Collection will be truncated ";
    }
  
    *(privateData_->ncand) = collection->size();

    typedef std::vector<double> jetVtxQualCollection;
    Handle<jetVtxQualCollection> JV_alpha;
    jetVtxQualCollection::const_iterator jetVtxAlphaItr;
    if(saveJetExtras_) {
      try { iEvent.getByLabel(jetVertexAlphaCollection_,JV_alpha); }
      catch ( cms::Exception& ex ) { 
	edm::LogWarning("CmsJetFiller") << "Can't get jet vertex alpha collection " 
				   << jetVertexAlphaCollection_;   
      }     
      jetVtxAlphaItr = JV_alpha->begin();
    }

    //    typedef AssociationVector<CandidateBaseRefProd, vector<float> > btagCollection;
    //    Handle<btagCollection> combinedSecondaryVertexBJetTags, 
    edm::Handle<reco::JetTagCollection> combinedSecondaryVertexBJetTags,
      combinedSecondaryVertexMVABJetTags,
      jetBProbabilityBJetTags,
      jetProbabilityBJetTags,
      simpleSecondaryVertexBJetTags,
      softElectronBJetTags,
      softMuonBJetTags,
      softMuonNoIPBJetTags,
      trackCountingHighPurBJetTags,
      trackCountingHighEffBJetTags;

    if(saveJetBTag_) {
      iEvent.getByLabel("newCombinedSecondaryVertexBJetTags", combinedSecondaryVertexBJetTags);
      iEvent.getByLabel("newCombinedSecondaryVertexMVABJetTags", combinedSecondaryVertexMVABJetTags);
      iEvent.getByLabel("newJetBProbabilityBJetTags", jetBProbabilityBJetTags);
      iEvent.getByLabel("newJetProbabilityBJetTags", jetProbabilityBJetTags);
      iEvent.getByLabel("newSimpleSecondaryVertexBJetTags", simpleSecondaryVertexBJetTags);
      iEvent.getByLabel("newSoftElectronBJetTags", softElectronBJetTags);
      iEvent.getByLabel("newSoftMuonBJetTags", softMuonBJetTags);
      iEvent.getByLabel("newSoftMuonNoIPBJetTags", softMuonNoIPBJetTags);
      iEvent.getByLabel("newTrackCountingHighPurBJetTags", trackCountingHighPurBJetTags);
      iEvent.getByLabel("newTrackCountingHighEffBJetTags", trackCountingHighEffBJetTags);
    }

    int index = 0;
    edm::View<reco::Candidate>::const_iterator cand;
    for(cand=collection->begin(); cand!=collection->end(); cand++) {
      const CaloJet *thisRecoJet = dynamic_cast< const CaloJet * > ( &(*cand) );
      // fill basic kinematics
      if(saveCand_) writeCandInfo(&(*cand),iEvent,iSetup);
      // fill jet extra informations
      if(saveJetExtras_) { 
	privateData_->alpha->push_back(*jetVtxAlphaItr);
	jetVtxAlphaItr++;

	// em, had fractions
	if( thisRecoJet != 0 ) { 
	  privateData_->emFrac->push_back( thisRecoJet->emEnergyFraction() );
	  privateData_->hadFrac->push_back( thisRecoJet->energyFractionHadronic() );
	}
	else {
	  privateData_->emFrac->push_back( -1.);
	  privateData_->hadFrac->push_back( -1.);
	}
      }
      else {
	privateData_->alpha->push_back( -1. );
	privateData_->emFrac->push_back( -1. );
	privateData_->hadFrac->push_back( -1. );
      }
      // fill the btag algorithms output
      if(saveJetBTag_) {
        privateData_->combinedSecondaryVertexBJetTags->push_back( (*combinedSecondaryVertexBJetTags)[index].second );
        privateData_->combinedSecondaryVertexMVABJetTags->push_back( (*combinedSecondaryVertexMVABJetTags)[index].second );
        privateData_->jetBProbabilityBJetTags->push_back( (*jetBProbabilityBJetTags)[index].second );
        privateData_->jetProbabilityBJetTags->push_back( (*jetProbabilityBJetTags)[index].second );
        privateData_->simpleSecondaryVertexBJetTags->push_back( (*simpleSecondaryVertexBJetTags)[index].second );
        privateData_->softElectronBJetTags->push_back( (*softElectronBJetTags)[index].second );
        privateData_->softMuonBJetTags->push_back( (*softMuonBJetTags)[index].second );
        privateData_->softMuonNoIPBJetTags->push_back( (*softMuonNoIPBJetTags)[index].second );
        privateData_->trackCountingHighPurBJetTags->push_back( (*trackCountingHighPurBJetTags)[index].second );
        privateData_->trackCountingHighEffBJetTags->push_back( (*trackCountingHighEffBJetTags)[index].second );
      } else {
        privateData_->combinedSecondaryVertexBJetTags->push_back( -1. );
        privateData_->combinedSecondaryVertexMVABJetTags->push_back( -1. );
        privateData_->jetBProbabilityBJetTags->push_back( -1. );
        privateData_->jetProbabilityBJetTags->push_back( -1. );
        privateData_->simpleSecondaryVertexBJetTags->push_back( -1. );
        privateData_->softElectronBJetTags->push_back( -1. );
        privateData_->softMuonBJetTags->push_back( 1. );
        privateData_->softMuonNoIPBJetTags->push_back( -1. );
        privateData_->trackCountingHighPurBJetTags->push_back( -1. );
        privateData_->trackCountingHighEffBJetTags->push_back( -1. );
      }
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
  if(saveJetExtras_) treeJetInfo(columnPrefix,columnSuffix);

  if(dumpData) cmstree->dumpData();
	
}









void CmsJetFiller::treeJetInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"alpha"+colSuffix).c_str(), *privateData_->alpha, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"emFrac"+colSuffix).c_str(), *privateData_->emFrac, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hadFrac"+colSuffix).c_str(), *privateData_->hadFrac, nCandString.c_str(), 0, "Reco");
  if(saveJetBTag_) {
    cmstree->column((colPrefix+"combinedSecondaryVertexBJetTags"+colSuffix).c_str(), *privateData_->combinedSecondaryVertexBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"combinedSecondaryVertexMVABJetTags"+colSuffix).c_str(), *privateData_->combinedSecondaryVertexMVABJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"jetBProbabilityBJetTags"+colSuffix).c_str(), *privateData_->jetBProbabilityBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"jetProbabilityBJetTags"+colSuffix).c_str(), *privateData_->jetProbabilityBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"simpleSecondaryVertexBJetTags"+colSuffix).c_str(), *privateData_->simpleSecondaryVertexBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"softElectronBJetTags"+colSuffix).c_str(), *privateData_->softElectronBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"softMuonBJetTags"+colSuffix).c_str(), *privateData_->softMuonBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"softMuonNoIPBJetTags"+colSuffix).c_str(), *privateData_->softMuonNoIPBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackCountingHighPurBJetTags"+colSuffix).c_str(), *privateData_->trackCountingHighPurBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackCountingHighEffBJetTags"+colSuffix).c_str(), *privateData_->trackCountingHighEffBJetTags, nCandString.c_str(), 0, "Reco");
  }
}







void CmsJetFillerData::initialise() {

  initialiseCandidate();
  alpha = new vector<float>;
  emFrac = new vector<float>;
  hadFrac = new vector<float>;
  combinedSecondaryVertexBJetTags = new vector<float>;
  combinedSecondaryVertexMVABJetTags = new vector<float>;
  jetBProbabilityBJetTags = new vector<float>;
  jetProbabilityBJetTags = new vector<float>;
  simpleSecondaryVertexBJetTags = new vector<float>;
  softElectronBJetTags = new vector<float>;
  softMuonBJetTags = new vector<float>;
  softMuonNoIPBJetTags = new vector<float>;
  trackCountingHighPurBJetTags = new vector<float>;
  trackCountingHighEffBJetTags = new vector<float>;

}




void CmsJetFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();
  alpha->clear();
  emFrac->clear();
  hadFrac->clear();
  combinedSecondaryVertexBJetTags->clear();
  combinedSecondaryVertexMVABJetTags->clear();
  jetBProbabilityBJetTags->clear();
  jetProbabilityBJetTags->clear();
  simpleSecondaryVertexBJetTags->clear();
  softElectronBJetTags->clear();
  softMuonBJetTags->clear();
  softMuonNoIPBJetTags->clear();
  trackCountingHighPurBJetTags->clear();
  trackCountingHighEffBJetTags->clear();

}
