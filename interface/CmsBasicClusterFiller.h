// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsBasicClusterFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Emanuele Di Marco
//         Created:  Fri Apr  6 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsBasicClusterFiller_h
#define CmsBasicClusterFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include <TTree.h>

struct CmsBasicClusterFillerData {
  // All the vectors that will store the stuff
  // going into the tree.
  vector<int> *nCrystals, *nOverlap3x3;
  vector<float> *energy, *eta, *phi;
  vector<float> *seedEnergy, *eMax, *e3x3, *e5x5;
  int *nSC;

public:
  void initialiseCandidate();
  void clear();
};


class CmsBasicClusterFiller {

public:

  // Constructors

  // Dump everything
  CmsBasicClusterFiller(CmsTree *, int maxSC=500);


  // Destructor
  virtual ~CmsBasicClusterFiller();

  // Operators

  // run number and all of that --- to implement

  virtual void writeCollectionToTree(edm::InputTag collectionTag,
				     const edm::Event&, const edm::EventSetup&,
				     const std::string &columnPrefix, const std::string &columnSuffix,
				     bool dumpData=false);

  //! set the rechits for ECAL barrel (needed for cluster shapes)
  void setEcalBarrelRecHits( edm::InputTag EcalBarrelRecHits ) { EcalBarrelRecHits_ = EcalBarrelRecHits; }
  //! set the rechits for ECAL endcap (needed for cluster shapes)
  void setEcalEndcapRecHits( edm::InputTag EcalEndcapRecHits ) { EcalEndcapRecHits_ = EcalEndcapRecHits; }


protected:
  

  virtual void writeSCInfo(const reco::BasicCluster *cand, 
			   const edm::Event&, const edm::EventSetup&,
                           const EcalRecHitCollection *EBRecHits, const EcalRecHitCollection *EERecHits);
  virtual void treeSCInfo(const std::string colPrefix, const std::string colSuffix);
  

  // Friends

  CmsBasicClusterFillerData *privateData_;

  edm::InputTag EcalBarrelRecHits_;
  edm::InputTag EcalEndcapRecHits_;

  CmsTree *cmstree;

  int maxSC_;
  std::string *trkIndexName_;
};

#endif // CmsBasicClusterFiller_h
