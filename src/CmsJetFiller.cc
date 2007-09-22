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
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchElectronFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsJetFiller.h"

#include "CLHEP/HepMC/GenParticle.h"


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
			   edm::InputTag jetVertexBetaCollection,
			   int maxTracks, int maxMCTracks,
			   bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsJetFillerData)
{
  cmstree=cmsTree;
  jetVertexAlphaCollection_=jetVertexAlphaCollection;
  jetVertexBetaCollection_=jetVertexBetaCollection;

  saveJetExtras_=true;

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();
}

CmsJetFiller::CmsJetFiller(CmsTree *cmsTree, 
			   edm::InputTag jetVertexAlphaCollection,
			   edm::InputTag jetVertexBetaCollection,
			   bool fatTree, 
			   int maxTracks, int maxMCTracks,
			   bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,fatTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsJetFillerData)
{
  cmstree=cmsTree;
  jetVertexAlphaCollection_=jetVertexAlphaCollection;
  jetVertexBetaCollection_=jetVertexBetaCollection;

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
  delete privateData_->beta;

}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out

void CmsJetFiller::saveJetExtras(bool what) { saveJetExtras_=what; }

void CmsJetFiller::writeCollectionToTree(const CandidateCollection *collection,
			   const edm::Event& iEvent, const edm::EventSetup& iSetup,
			   const std::string &columnPrefix, const std::string &columnSuffix,
			   bool dumpData) {
  
  if(hitLimitsMeansNoOutput_ && 
     (int)collection->size() > maxTracks_){
    LogError("CmsJetFiller") << "Track length " << collection->size() 
			     << " is too long for declared max length for tree "
			     << maxTracks_ << " and no output flag is set."
			     << " No tracks written to tuple for this event ";
    return;
  }

  if((int)collection->size() > maxTracks_){
    LogError("CmsJetFiller") << "Track length " << collection->size() 
			     << " is too long for declared max length for tree "
			     << maxTracks_ 
			     << ". Collection will be truncated ";
  }
  
  *(privateData_->ncand) = collection->size();

  LogInfo("CmsJetFiller") << "Filling candidate vectors";
  privateData_->clearTrkVectors();

  typedef std::vector<double> jetVtxQualCollection;
  Handle<jetVtxQualCollection> JV_alpha;
  iEvent.getByLabel(jetVertexAlphaCollection_,JV_alpha);
  Handle<jetVtxQualCollection> JV_beta;
  iEvent.getByLabel(jetVertexBetaCollection_,JV_beta);

  jetVtxQualCollection::const_iterator jetVtxAlphaItr = JV_alpha->begin();
  jetVtxQualCollection::const_iterator jetVtxBetaItr = JV_beta->begin();
  CandidateCollection::const_iterator cand;

  for(cand=collection->begin(); cand!=collection->end(); cand++) {
    // fill basic kinematics
    if(saveCand_) writeCandInfo(&(*cand),iEvent,iSetup);

    // fill jet extra informations
    privateData_->alpha->push_back(*jetVtxAlphaItr);
    privateData_->beta->push_back(*jetVtxBetaItr);
    jetVtxAlphaItr++;
    jetVtxBetaItr++;

    // em, had fractions
    CaloJetRef thisRecoJet = cand->masterClone().castTo<CaloJetRef>();

    if( thisRecoJet.isNonnull() ) { 
      privateData_->emFrac->push_back( thisRecoJet->emEnergyFraction() );
      privateData_->hadFrac->push_back( thisRecoJet->energyFractionHadronic() );
    }
    else {
      privateData_->emFrac->push_back( -1. );
      privateData_->hadFrac->push_back( -1. );
    }

  }
  
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the 
  // tree 
  
  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),collection->size(),0,"Reco");
  
  if(saveJetExtras_) treeJetInfo(columnPrefix,columnSuffix);
  if(dumpData) cmstree->dumpData();

}

void CmsJetFiller::treeJetInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"alpha"+colSuffix).c_str(), *privateData_->alpha, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"beta"+colSuffix).c_str(), *privateData_->beta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"emFrac"+colSuffix).c_str(), *privateData_->emFrac, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hadFrac"+colSuffix).c_str(), *privateData_->hadFrac, nCandString.c_str(), 0, "Reco");

}
void CmsJetFillerData::initialise() {

  alpha = new vector<float>;
  beta = new vector<float>;
  emFrac = new vector<float>;
  hadFrac = new vector<float>;

}

void CmsJetFillerData::clearTrkVectors() {

  alpha->clear();
  beta->clear();
  emFrac->clear();
  hadFrac->clear();

}
