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

#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
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
			   int maxTracks, int maxMCTracks,
			   bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsJetFillerData)
{
  cmstree=cmsTree;

  saveJetExtras_=true;

  saveJetBTag_ = false;
  
  trkIndexName_ = new std::string("n");
  
  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();
}

CmsJetFiller::CmsJetFiller(CmsTree *cmsTree, 
			   bool fatTree, 
			   int maxTracks, int maxMCTracks,
			   bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,fatTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsJetFillerData)
{
  cmstree=cmsTree;

  trkIndexName_ = new std::string("n");

  saveJetExtras_=true;

  saveJetBTag_ = false;

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
  delete privateData_->Id;
  delete privateData_->nHit;
  delete privateData_->nHit90;
  delete privateData_->fHPD;
  delete privateData_->covEtaEta;
  delete privateData_->covPhiPhi;
  delete privateData_->fLS;
  delete privateData_->fOOT;
  delete privateData_->combinedSecondaryVertexBJetTags;
  delete privateData_->combinedSecondaryVertexMVABJetTags;
  delete privateData_->jetBProbabilityBJetTags;
  delete privateData_->jetProbabilityBJetTags;
  delete privateData_->simpleSecondaryVertexHighEffBJetTags;
  delete privateData_->simpleSecondaryVertexHighPurBJetTags;
  delete privateData_->softMuonBJetTags;
  delete privateData_->softMuonByIP3dBJetTags;
  delete privateData_->softMuonByPtBJetTags;
  //  delete privateData_->softElectronBJetTags;
  delete privateData_->softElectronByIP3dBJetTags;
  delete privateData_->softElectronByPtBJetTags;
  delete privateData_->trackCountingHighPurBJetTags;
  delete privateData_->trackCountingHighEffBJetTags;
  delete privateData_->trackCountingVeryHighEffBJetTags;
  delete privateData_->uncorrEnergy;
  delete privateData_->L2L3CorrEnergy;
  delete privateData_->area;

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
					 bool dumpData,
                                         edm::InputTag uncorrectedCollectionTag, edm::InputTag L2L3correctedCollectionTag) {

  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsJetFiller") << "Can't get candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();

  const edm::View<reco::CaloJet> *uncorrectedCollection = 0;
  if(uncorrectedCollectionTag.label() != std::string("")) {
    edm::Handle< edm::View<reco::CaloJet> > uncorrectedCollectionHandle;
    try { iEvent.getByLabel(uncorrectedCollectionTag, uncorrectedCollectionHandle); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsJetFiller") << "Can't get candidate collection: " << uncorrectedCollectionTag; }
    uncorrectedCollection = uncorrectedCollectionHandle.product();
  }

  const edm::View<reco::CaloJet> *L2L3correctedCollection = 0;
  if(L2L3correctedCollectionTag.label() != std::string("")) {
    edm::Handle< edm::View<reco::CaloJet> > L2L3correctedCollectionHandle;
    try { iEvent.getByLabel(L2L3correctedCollectionTag, L2L3correctedCollectionHandle); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsJetFiller") << "Can't get candidate collection: " << L2L3correctedCollectionTag; }
    L2L3correctedCollection = L2L3correctedCollectionHandle.product();
  }
  
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

    edm::Handle<reco::JetTagCollection> combinedSecondaryVertexBJetTags,
      combinedSecondaryVertexMVABJetTags,
      jetBProbabilityBJetTags,
      jetProbabilityBJetTags,
      simpleSecondaryVertexHighEffBJetTags,
      simpleSecondaryVertexHighPurBJetTags,
      softMuonBJetTags,
      softMuonByIP3dBJetTags,
      softMuonByPtBJetTags,
      //      softElectronBJetTags,
      softElectronByIP3dBJetTags,
      softElectronByPtBJetTags,
      trackCountingHighPurBJetTags,
      trackCountingHighEffBJetTags,
      trackCountingVeryHighEffBJetTags;

    if(saveJetBTag_) {
      iEvent.getByLabel("newCombinedSecondaryVertexBJetTags", combinedSecondaryVertexBJetTags);
      iEvent.getByLabel("newCombinedSecondaryVertexMVABJetTags", combinedSecondaryVertexMVABJetTags);
      iEvent.getByLabel("newJetBProbabilityBJetTags", jetBProbabilityBJetTags);
      iEvent.getByLabel("newJetProbabilityBJetTags", jetProbabilityBJetTags);
      iEvent.getByLabel("newSimpleSecondaryVertexHighEffBJetTags", simpleSecondaryVertexHighEffBJetTags);
      iEvent.getByLabel("newSimpleSecondaryVertexHighPurBJetTags", simpleSecondaryVertexHighPurBJetTags);
      iEvent.getByLabel("newSoftMuonBJetTags", softMuonBJetTags);
      iEvent.getByLabel("newSoftMuonByIP3dBJetTags", softMuonByIP3dBJetTags);
      iEvent.getByLabel("newSoftMuonByPtBJetTags", softMuonByPtBJetTags);
      // iEvent.getByLabel("newSoftElectronBJetTags", softElectronBJetTags);
      iEvent.getByLabel("newSoftElectronByIP3dBJetTags", softElectronByIP3dBJetTags);
      iEvent.getByLabel("newSoftElectronByPtBJetTags", softElectronByPtBJetTags);
      iEvent.getByLabel("newTrackCountingHighPurBJetTags", trackCountingHighPurBJetTags);
      iEvent.getByLabel("newTrackCountingHighEffBJetTags", trackCountingHighEffBJetTags);
      iEvent.getByLabel("newTrackCountingVeryHighEffBJetTags", trackCountingVeryHighEffBJetTags);
    }

    edm::Handle<reco::JetIDValueMap> hJetIDMap;
    if(saveJetExtras_) {
      if ( iEvent.getByLabel( "ak5JetID", hJetIDMap ) ) {}
      else edm::LogWarning("CmsJetFiller") << "Can't get ak5JetID collection." << endl;
    }

    JetIDSelectionFunctor jetIDFunctorPURE1( JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::MINIMAL );
    JetIDSelectionFunctor jetIDFunctorPURE2( JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::LOOSE );
    JetIDSelectionFunctor jetIDFunctorPURE3( JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::TIGHT );
    JetIDSelectionFunctor jetIDFunctorDQM1( JetIDSelectionFunctor::DQM09, JetIDSelectionFunctor::MINIMAL );
    JetIDSelectionFunctor jetIDFunctorDQM2( JetIDSelectionFunctor::DQM09, JetIDSelectionFunctor::LOOSE );
    JetIDSelectionFunctor jetIDFunctorDQM3( JetIDSelectionFunctor::DQM09, JetIDSelectionFunctor::TIGHT );

    int index = 0;
    edm::View<reco::Candidate>::const_iterator cand;
    for(cand=collection->begin(); cand!=collection->end(); cand++) {
      const CaloJet *thisRecoJet = dynamic_cast< const CaloJet * > ( &(*cand) );
      // fill basic kinematics
      if(saveCand_) writeCandInfo(&(*cand),iEvent,iSetup);
      // fill jet extra informations
      if(saveJetExtras_) { 

	// em, had fractions and Jet ID
	if( thisRecoJet != 0 ) { 
	  privateData_->emFrac->push_back( thisRecoJet->emEnergyFraction() );
	  privateData_->hadFrac->push_back( thisRecoJet->energyFractionHadronic() );
          privateData_->area->push_back( thisRecoJet->jetArea() );
	}
	else {
	  privateData_->emFrac->push_back( -1.);
	  privateData_->hadFrac->push_back( -1.);
          privateData_->area->push_back( -1. );
	}
      }
      else {
	privateData_->emFrac->push_back( -1. );
	privateData_->hadFrac->push_back( -1. );
        privateData_->area->push_back( -1. );
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
        //        privateData_->softElectronBJetTags->push_back( (*softElectronBJetTags)[index].second );
        privateData_->softElectronByIP3dBJetTags->push_back( (*softElectronByIP3dBJetTags)[index].second );
        privateData_->softElectronByPtBJetTags->push_back( (*softElectronByPtBJetTags)[index].second );
        privateData_->trackCountingHighPurBJetTags->push_back( (*trackCountingHighPurBJetTags)[index].second );
        privateData_->trackCountingHighEffBJetTags->push_back( (*trackCountingHighEffBJetTags)[index].second );
        privateData_->trackCountingVeryHighEffBJetTags->push_back( (*trackCountingVeryHighEffBJetTags)[index].second );
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
        //        privateData_->softElectronBJetTags->push_back( -1. );
        privateData_->softElectronByIP3dBJetTags->push_back( -1. );
        privateData_->softElectronByPtBJetTags->push_back( -1. );
        privateData_->trackCountingHighPurBJetTags->push_back( -1. );
        privateData_->trackCountingHighEffBJetTags->push_back( -1. );
        privateData_->trackCountingVeryHighEffBJetTags->push_back( -1. );
      }

      // if an uncorrected jet collection is provided, save also the uncorrected energy
      if(uncorrectedCollection) {
        dumpUncorrEnergy_ = true;
        float rawEnergy = -999.;
        float fHPD, covEtaEta, covPhiPhi, fLS, fOOT;
        fHPD = covEtaEta = covPhiPhi = fLS = fOOT = -999.;
        int Id, nHit, nHit90;
        Id = nHit = nHit90 = -999;
        edm::View<reco::CaloJet>::const_iterator cand2;
        for(cand2=uncorrectedCollection->begin(); cand2!=uncorrectedCollection->end(); cand2++) {
          const CaloJet *uncorrectedRecoJet = dynamic_cast< const CaloJet * > ( &(*cand2) );
          // corrected and uncorrected jets differ only for jet PT 
          if( thisRecoJet->jetArea() == uncorrectedRecoJet->jetArea() ) {

            rawEnergy = uncorrectedRecoJet->energy();

            if(saveJetExtras_) { // jet ID has to be done on uncorrected jets
              unsigned int idx = cand2 - uncorrectedCollection->begin();
              edm::RefToBase<reco::CaloJet> jetRef = uncorrectedCollection->refAt(idx);
              reco::JetID const & jetId = (*hJetIDMap)[ jetRef ];

              fHPD = jetId.fHPD;
              fLS = jetId.fLS;
              fOOT = jetId.fHFOOT;
              nHit = jetIDFunctorPURE1.count_hits( cand2->getCaloConstituents() );
              nHit90 = jetId.n90Hits;
              covEtaEta = cand2->etaetaMoment();
              covPhiPhi = cand2->phiphiMoment();

              int passedPURE1 = ( jetIDFunctorPURE1( *cand2, jetId ) ) ? 1 : 0;
              int passedPURE2 = ( jetIDFunctorPURE2( *cand2, jetId ) ) ? 1 : 0;
              int passedPURE3 = ( jetIDFunctorPURE3( *cand2, jetId ) ) ? 1 : 0;
              int passedDQM1 = ( jetIDFunctorDQM1( *cand2, jetId ) ) ? 1 : 0;
              int passedDQM2 = ( jetIDFunctorDQM2( *cand2, jetId ) ) ? 1 : 0;
              int passedDQM3 = ( jetIDFunctorDQM3( *cand2, jetId ) ) ? 1 : 0;

              Id = ( passedPURE3 << 5 ) | ( passedPURE2 << 4 ) | ( passedPURE1 << 3 ) |
                ( passedDQM3 << 2 ) | ( passedDQM2 << 1 ) | passedDQM1; 

            }
            break;
          }
        }
        privateData_->uncorrEnergy->push_back(rawEnergy);
        privateData_->Id->push_back(Id);
        privateData_->nHit->push_back(nHit);
        privateData_->nHit90->push_back(nHit90);
        privateData_->fHPD->push_back(fHPD);
        privateData_->covEtaEta->push_back(covEtaEta);
        privateData_->covPhiPhi->push_back(covPhiPhi);
        privateData_->fLS->push_back(fLS);
        privateData_->fOOT->push_back(fOOT);
      } else dumpUncorrEnergy_ = false;


      // if an L2L3corrected jet collection is provided, save also the L2L3corrected energy
      if(L2L3correctedCollection) {
        dumpL2L3CorrEnergy_ = true;
        float L2L3Energy = -999.;
        edm::View<reco::CaloJet>::const_iterator cand2;
        for(cand2=L2L3correctedCollection->begin(); cand2!=L2L3correctedCollection->end(); cand2++) {
          const CaloJet *L2L3correctedRecoJet = dynamic_cast< const CaloJet * > ( &(*cand2) );
          // corrected and uncorrected jets differ only for jet PT 
          if( thisRecoJet->jetArea() == L2L3correctedRecoJet->jetArea() ) {

            L2L3Energy = L2L3correctedRecoJet->energy();
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
  if(saveJetExtras_) treeJetInfo(columnPrefix,columnSuffix);

  if(dumpData) cmstree->dumpData();
	
}









void CmsJetFiller::treeJetInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"emFrac"+colSuffix).c_str(), *privateData_->emFrac, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hadFrac"+colSuffix).c_str(), *privateData_->hadFrac, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"area"+colSuffix).c_str(), *privateData_->area, nCandString.c_str(), 0, "Reco");
  if(saveJetExtras_) {
    cmstree->column((colPrefix+"Id"+colSuffix).c_str(), *privateData_->Id, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"nHit"+colSuffix).c_str(), *privateData_->nHit, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"nHit90"+colSuffix).c_str(), *privateData_->nHit90, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"fHPD"+colSuffix).c_str(), *privateData_->fHPD, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"covEtaEta"+colSuffix).c_str(), *privateData_->covEtaEta, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"covPhiPhi"+colSuffix).c_str(), *privateData_->covPhiPhi, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"fLS"+colSuffix).c_str(), *privateData_->fLS, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"fOOT"+colSuffix).c_str(), *privateData_->fOOT, nCandString.c_str(), 0, "Reco");
  }
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
    //    cmstree->column((colPrefix+"softElectronBJetTags"+colSuffix).c_str(), *privateData_->softElectronBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"softElectronByIP3dBJetTags"+colSuffix).c_str(), *privateData_->softElectronByIP3dBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"softElectronByPtBJetTags"+colSuffix).c_str(), *privateData_->softElectronByPtBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackCountingHighPurBJetTags"+colSuffix).c_str(), *privateData_->trackCountingHighPurBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackCountingHighEffBJetTags"+colSuffix).c_str(), *privateData_->trackCountingHighEffBJetTags, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackCountingVeryHighEffBJetTags"+colSuffix).c_str(), *privateData_->trackCountingVeryHighEffBJetTags, nCandString.c_str(), 0, "Reco");
  }
  if(dumpUncorrEnergy_) {
    cmstree->column((colPrefix+"uncorrEnergy"+colSuffix).c_str(), *privateData_->uncorrEnergy, nCandString.c_str(), 0, "Reco");
  }
  if(dumpL2L3CorrEnergy_) {
    cmstree->column((colPrefix+"L2L3CorrEnergy"+colSuffix).c_str(), *privateData_->L2L3CorrEnergy, nCandString.c_str(), 0, "Reco");
  }
}







void CmsJetFillerData::initialise() {
  initialiseCandidate();
  emFrac = new vector<float>;
  hadFrac = new vector<float>;
  area = new vector<float>;
  Id = new vector<int>;
  nHit = new vector<int>;
  nHit90 = new vector<int>;
  fHPD = new vector<float>;
  covEtaEta = new vector<float>;
  covPhiPhi = new vector<float>;
  fLS = new vector<float>;
  fOOT = new vector<float>;
  combinedSecondaryVertexBJetTags = new vector<float>;
  combinedSecondaryVertexMVABJetTags = new vector<float>;
  jetBProbabilityBJetTags = new vector<float>;
  jetProbabilityBJetTags = new vector<float>;
  simpleSecondaryVertexHighEffBJetTags = new vector<float>;
  simpleSecondaryVertexHighPurBJetTags = new vector<float>;
  softMuonBJetTags = new vector<float>;
  softMuonByIP3dBJetTags = new vector<float>;
  softMuonByPtBJetTags = new vector<float>;
  //  softElectronBJetTags = new vector<float>;
  softElectronByIP3dBJetTags = new vector<float>;
  softElectronByPtBJetTags = new vector<float>;
  trackCountingHighPurBJetTags = new vector<float>;
  trackCountingHighEffBJetTags = new vector<float>;
  trackCountingVeryHighEffBJetTags = new vector<float>;
  uncorrEnergy = new vector<float>;
  L2L3CorrEnergy = new vector<float>;
}




void CmsJetFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();
  emFrac->clear();
  hadFrac->clear();
  area->clear();
  Id->clear();
  nHit->clear();
  nHit90->clear();
  fHPD->clear();
  covEtaEta->clear();
  covPhiPhi->clear();
  fLS->clear();
  fOOT->clear();
  combinedSecondaryVertexBJetTags->clear();
  combinedSecondaryVertexMVABJetTags->clear();
  jetBProbabilityBJetTags->clear();
  jetProbabilityBJetTags->clear();
  simpleSecondaryVertexHighEffBJetTags->clear();
  simpleSecondaryVertexHighPurBJetTags->clear();
  softMuonBJetTags->clear();
  softMuonByIP3dBJetTags->clear();
  softMuonByPtBJetTags->clear();
  //  softElectronBJetTags->clear();
  softElectronByIP3dBJetTags->clear();
  softElectronByPtBJetTags->clear();
  trackCountingHighPurBJetTags->clear();
  trackCountingHighEffBJetTags->clear();
  trackCountingVeryHighEffBJetTags->clear();
  uncorrEnergy->clear();
  L2L3CorrEnergy->clear();
}
