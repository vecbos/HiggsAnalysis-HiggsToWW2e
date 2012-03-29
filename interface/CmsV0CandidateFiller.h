// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsV0CandidateFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Emanuele Di Marco
//         Created:  Fri Apr  6 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsV0CandidateFiller_h
#define CmsV0CandidateFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include <TTree.h>

struct CmsV0CandidateFillerData : public CmsCandidateFillerData {

  vector<int> *dau1Index, *dau2Index;

public:
  void initialise();
  void clearTrkVectors();

};

class CmsV0CandidateFiller : public CmsCandidateFiller {

 public:

  // Constructors

  // Dump everything
  CmsV0CandidateFiller(CmsTree *, 
		int maxTracks=500, int maxMCTracks=2000, 
		bool noOutputIfLimitsReached=false );

  // Dump  everything if fatTree is true and less informations otherwise
  CmsV0CandidateFiller(CmsTree *, 
		bool fatTree, 
		int maxTracks=500, int maxMCTracks=2000, 
		bool noOutputIfLimitsReached=false );
  
  // Destructor
  virtual ~CmsV0CandidateFiller();

  //! write the V0 related informations for the given collection
  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);

  
 private:
  
  void writeTrkInfo(const reco::Candidate *cand, const edm::Event&, const edm::EventSetup&);
  void treeTrkInfo(const std::string &colPrefix, const std::string &colSuffix);

  bool saveTrk_;

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;

  std::string *trkIndexName_;

  CmsV0CandidateFillerData *privateData_;

  CmsTree *cmstree;

};

#endif // CmsV0CandidateFiller_h
