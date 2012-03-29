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

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include <TTree.h>

struct CmsCandidateFillerData {
  // All the vectors that will store the stuff
  // going into the tree.
  vector<int> *charge;
  vector<float> *energy, *et, *momentum;
  vector<float> *vertexX, *vertexY, *vertexZ;
  vector<float> *theta, *eta, *phi;
  vector<float> *x, *y, *z;
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

  //! write the basic candidate informations for the collection "collection"
  virtual void writeCollectionToTree(edm::InputTag collection,
				     const edm::Event&, const edm::EventSetup&,
				     const std::string &columnPrefix, const std::string &columnSuffix,
				     bool dumpData=false);
  
 protected:
  

  virtual void writeCandInfo(const reco::Candidate *cand, 
                             const edm::Event&, const edm::EventSetup&);
  virtual void treeCandInfo(const std::string colPrefix, const std::string colSuffix);
  
  // Friends

  CmsCandidateFillerData *privateData_;

  CmsTree *cmstree;

  bool saveCand_;

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  std::string *trkIndexName_;
};

#endif // CmsCandidateFiller_h
