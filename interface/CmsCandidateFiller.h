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

  // Constructors

  // Dump everything
  CmsCandidateFiller(CmsTree *, int maxTracks=500, 
		    int maxMCTracks=2000, bool noOutputIfLimitsReached=false );

  // Dump  everything if fatTree is true and less informations otherwise
  CmsCandidateFiller(CmsTree *, bool fatTree, int maxTracks=500, 
		     int maxMCTracks=2000, bool noOutputIfLimitsReached=false );

  // Destructor
  virtual ~CmsCandidateFiller();

  // Modifiers


  void saveCand(bool );
  void doMcMatch(bool );

  void setMatchMap(edm::InputTag matchMap) { matchMap_ = matchMap;};

  // Operators

  // run number and all of that --- to implement

  virtual void writeCollectionToTree(const edm::View<reco::Candidate> *,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);

  virtual void writeMcIndicesToTree(const edm::View<reco::Candidate> *,
			    const edm::Event&, const edm::EventSetup&,
			    const edm::View<reco::Candidate> *,
			    const std::string &columnPrefix, const std::string &columnSuffix,
			    bool dumpData=false);


  // Add daughter list
  virtual void addDaughterList(const edm::View<reco::Candidate> *);

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
