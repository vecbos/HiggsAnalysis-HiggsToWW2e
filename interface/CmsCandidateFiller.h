// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsCandidateFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Emanuele Di Marco
//         Created:  Fri Apr  6 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsCandidateFiller_h
#define CmsCandidateFiller_h

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "CLHEP/HepMC/GenEvent.h"
#include <TTree.h>

using namespace cms;
using namespace edm;
using namespace reco;

struct CmsCandidateFillerData {
  // All the vectors that will store the stuff
  // going into the tree.
  vector<int> *charge;
  vector<float> *energy, *et, *momentum;
  vector<float> *vertexX, *vertexY, *vertexZ;
  vector<float> *theta, *eta, *phi;
  vector<float> *x, *y, *z;
  vector<float> *mass, *mt;
  vector<int> *pdgId;
  vector<int> *nDau;
  vector<int> *d1Index, *d2Index;
  vector<int> *d1pdgId, *d2pdgId;

  vector<int> *mcIndex;
  int *ncand;

public:
  void initialiseCandidate();
  void clearTrkVectorsCandidate();
};


class CmsCandidateFiller {

 public:

  //! Dump everything
  CmsCandidateFiller(CmsTree *, int maxTracks=500, 
		    int maxMCTracks=2000, bool noOutputIfLimitsReached=false );

  //! Dump  everything if fatTree is true and less informations otherwise
  CmsCandidateFiller(CmsTree *, bool fatTree, int maxTracks=500, 
		     int maxMCTracks=2000, bool noOutputIfLimitsReached=false );

  //! Destructor
  virtual ~CmsCandidateFiller();

  //! dump the particle candidate informations
  void saveCand(bool );
  //! do CM truth - reco particle matching, based on delta R
  void doMcMatch(bool );
  //! set the association map between reco - truth particles
  void setMatchMap(edm::InputTag matchMap) { matchMap_ = matchMap;};

  //! write the basic candidate informations for the collection "collection"
  virtual void writeCollectionToTree(edm::InputTag collection,
				     const edm::Event&, const edm::EventSetup&,
				     const std::string &columnPrefix, const std::string &columnSuffix,
				     bool dumpData=false);
  
  //! write a column with indices connecting reco - MC truth particles in the tree
  virtual void writeMcIndicesToTree(edm::InputTag recoCollection,
				    const edm::Event&, const edm::EventSetup&,
				    edm::InputTag mcTruthCollection,
				    const std::string &columnPrefix, const std::string &columnSuffix,
				    bool dumpData=false);


  //! add a collection where to look for daughters of a composite collection 
  virtual void addDaughterCollection(edm::InputTag daugCollectionTag,
				     const edm::Event& iEvent, const edm::EventSetup& iSetup);

 protected:
  

  virtual void writeCandInfo(const Candidate *cand, 
		     const edm::Event&, const edm::EventSetup&);
  virtual void treeCandInfo(const std::string colPrefix, const std::string colSuffix);
  
  virtual void writeMcMatchInfo(const edm::View<reco::Candidate> *, 
			const edm::Event&, const edm::EventSetup&, 
			const edm::View<reco::Candidate> *);
  virtual void treeMcMatchInfo(const std::string colPrefix, const std::string colSuffix);

  // Friends

  CmsCandidateFillerData *privateData_;
  edm::InputTag matchMap_;
  std::vector< const edm::View<reco::Candidate>* > daugCollectionList_;

  CmsTree *cmstree;

  bool saveCand_;
  bool doMcMatch_;

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  std::string *trkIndexName_;
};

#endif // CmsCandidateFiller_h
