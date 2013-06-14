// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsCalibElectronFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Chiara Rovelli
//
//-----------------------------------------------------------------------

#ifndef CmsCalibElectronFiller_h
#define CmsCalibElectronFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include <TTree.h>

struct CmsCalibElectronFillerData : public CmsCandidateFillerData {

  // All the vectors that will store the stuff
  // going into the tree.
  vector<float> *energyError;
  
public:
  void initialise();
  void clearTrkVectors();
};

class CmsCalibElectronFiller : public CmsCandidateFiller {

 public:

  //! Dump everything
  CmsCalibElectronFiller(CmsTree *, int maxTracks=500,
			 int maxMCTracks=2000, bool noOutputIfLimitsReached=false );


  //! Destructor
  virtual ~CmsCalibElectronFiller();
  
  //! write the electron related informations for the given collection
  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);

 private:
  
  void treeExtraEleInfo(const std::string &colPrefix, const std::string &colSuffix);

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  std::string *trkIndexName_;

  CmsCalibElectronFillerData *privateData_;

  CmsTree *cmstree;
};

#endif // CmsCalibElectronFiller_h
