// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsMetFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Emanuele Di Marco
//         Created:  Fri Apr  6 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsPFCandidateFiller_h
#define CmsPFCandidateFiller_h

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
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include <TTree.h>

struct CmsPFCandidateFillerData : public CmsCandidateFillerData {
  // All the vectors that will store the stuff
  // going into the tree.
  vector<int> *particleId;

public:
  void initialise();
  void clearTrkVectors();
};


class CmsPFCandidateFiller : public CmsCandidateFiller {

 public:

  //! Dump everything
  CmsPFCandidateFiller(CmsTree *, int maxTracks=5000, 
		    int maxMCTracks=5000, bool noOutputIfLimitsReached=false );

  // Destructor
  virtual ~CmsPFCandidateFiller();
    
  //! set the general tracks (reduced collection)
  void setGeneralTracks( edm::InputTag generalTracks) { generalTracks_ = generalTracks; }

  //! write the basic candidate informations for the collection "collection"
  virtual void writeCollectionToTree(edm::InputTag collection,
				     const edm::Event&, const edm::EventSetup&,
				     const std::string &columnPrefix, const std::string &columnSuffix,
				     bool dumpData=false);
  

 protected:
  

  virtual void treePFInfo(const std::string colPrefix, const std::string colSuffix);
  

  // Friends

  CmsPFCandidateFillerData *privateData_;

  CmsTree *cmstree;

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  edm::InputTag generalTracks_;
  edm::Handle< reco::TrackRefVector > h_tracks;

  std::string *trkIndexName_;
};

#endif // CmsPFCandidateFiller_h
