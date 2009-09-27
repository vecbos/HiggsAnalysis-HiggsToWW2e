// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsCaloTowerFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Emanuele Di Marco
//         Created:  Fri Apr  6 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsCaloTowerFiller_h
#define CmsCaloTowerFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include <TTree.h>

struct CmsCaloTowerFillerData : public CmsCandidateFillerData {

  std::vector<float> *energy, *x, *y, *z;
  std::vector<int> *CALO, *CaloIndex;
  //std::vector<float> *emFrac, *hadFrac;
  //std::vector<float> *alpha;

public:
  void initialise();
  void clearTrkVectors();

};

class CmsCaloTowerFiller : public CmsCandidateFiller {

 public:

  // Constructors

  // Dump everything
  CmsCaloTowerFiller(CmsTree *, 
		     edm::InputTag hbheLabel,
		     edm::InputTag hoLabel,
		     edm::InputTag hfLabel,
		     std::vector<edm::InputTag> ecalLabels,
		     int maxTracks=500, int maxMCTracks=2000, 
		     bool noOutputIfLimitsReached=false );

  // Dump  everything if fatTree is true and less informations otherwise
  CmsCaloTowerFiller(CmsTree *, 
		     edm::InputTag hbheLabel,
		     edm::InputTag hoLabel,
		     edm::InputTag hfLabel,
		     std::vector<edm::InputTag> ecalLabels,
		     bool fatTree, 
		     int maxTracks=500, int maxMCTracks=2000, 
		     bool noOutputIfLimitsReached=false );
  
  // Destructor
  virtual ~CmsCaloTowerFiller();

  // Modifiers
  void saveCaloTowerExtras(bool );

  // Operators

  // run number and all of that --- to implement

  void writeCollectionToTree(edm::InputTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);

 private:
  
  void writeCaloTowerInfo(const Candidate *cand, const edm::Event&, const edm::EventSetup&);
  void treeCaloTowerInfo(const std::string &colPrefix, const std::string &colSuffix);

  // Friends
  bool saveCaloTowerExtras_;

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  std::string *trkIndexName_;

  CmsCaloTowerFillerData *privateData_;
  edm::InputTag hbheLabel_;
  edm::InputTag hoLabel_;
  edm::InputTag hfLabel_;
  std::vector<edm::InputTag> ecalLabels_;

  CmsTree *cmstree;

};

#endif // CmsCaloTowerFiller_h
