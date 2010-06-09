// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsPFPreIdFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//
//-----------------------------------------------------------------------

#ifndef CmsPFPreIdFiller_h
#define CmsPFPreIdFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include <TTree.h>

struct CmsPFPreIdFillerData {

  int *ncand;
  std::vector<int> *trackIndex;
  std::vector<float> *deltaEtaMatch, *deltaPhiMatch; 
  std::vector<float> *chiEtaMatch, *chiPhiMatch, *chi2Match;
  std::vector<float> *eopMatch;
  std::vector<float> *kfChi2, *kfNHits;
  std::vector<float> *gsfChi2, *chi2Ratio;
  std::vector<bool> *ecalMatching, *psMatching;
  std::vector<bool> *trackFiltered, *preided;

  public:
  void initialise();
  void clearTrkVectors();

};

class CmsPFPreIdFiller {

 public:

  // Constructors

  // Dump everything
  CmsPFPreIdFiller(CmsTree *, 
		      int maxTracks=500, int maxMCTracks=2000, 
		      bool noOutputIfLimitsReached=false );

  // Dump  everything if fatTree is true and less informations otherwise
  CmsPFPreIdFiller(CmsTree *, 
		      bool fatTree, 
		      int maxTracks=500, int maxMCTracks=2000, 
		      bool noOutputIfLimitsReached=false );

  // Destructor
  virtual ~CmsPFPreIdFiller();

  // Operators

  // run number and all of that --- to implement
  void writeCollectionToTree(edm::InputTag collectionTag1, edm::InputTag collectionTag2, 
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);

 private:

  void treePFPreIdInfo(const std::string &colPrefix, const std::string &colSuffix);
  
  // Friends
  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  std::string *trkIndexName_;
  
  CmsPFPreIdFillerData *privateData_;
  
  CmsTree *cmstree;

};

#endif // CmsPFPreIdFiller_h
