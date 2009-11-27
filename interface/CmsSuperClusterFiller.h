// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsSuperClusterFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Emanuele Di Marco
//         Created:  Fri Apr  6 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsSuperClusterFiller_h
#define CmsSuperClusterFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include <TTree.h>

struct CmsSuperClusterFillerData {
  // All the vectors that will store the stuff
  // going into the tree.
  vector<int>  *nBC, *nCrystals, *iAlgo, *recoFlag;
  vector<float> *rawEnergy, *energy, *seedEnergy, *eta, *theta, *phi, *time, *chi2Prob;
  vector<float> *e3x3, * e5x5, *eMax, *e2x2, *e2nd, *covIEtaIEta, *covIEtaIPhi, *covIPhiIPhi;
  vector<int> *trackIndex;
  vector<float> *deltaR, *deltaPhi, *deltaEta;
  int *nSC;

public:
  void initialiseCandidate();
  void clear();
};


class CmsSuperClusterFiller {

public:

  // Constructors

  // Dump everything
  CmsSuperClusterFiller(CmsTree *, int maxSC=500);


  // Destructor
  virtual ~CmsSuperClusterFiller();

  // Operators

  // run number and all of that --- to implement

  virtual void writeCollectionToTree(edm::InputTag collectionTag,
				     const edm::Event&, const edm::EventSetup&,
				     const std::string &columnPrefix, const std::string &columnSuffix,
				     bool dumpData=false);

  //! set the track collection (to match the superclusters as hand-made electrons)
  void setTracks( edm::InputTag Tracks ) { Tracks_ = Tracks; }
  //! set the rechits for ECAL barrel (needed for cluster shapes)
  void setEcalBarrelRecHits( edm::InputTag EcalBarrelRecHits ) { EcalBarrelRecHits_ = EcalBarrelRecHits; }
  //! set the rechits for ECAL endcap (needed for cluster shapes)
  void setEcalEndcapRecHits( edm::InputTag EcalEndcapRecHits ) { EcalEndcapRecHits_ = EcalEndcapRecHits; }


protected:
  

  virtual void writeSCInfo(const reco::SuperCluster *cand, 
			   const edm::Event&, const edm::EventSetup&,
                           const EcalRecHitCollection *EBRecHits, const EcalRecHitCollection *EERecHits, 
                           const reco::TrackCollection *Tracks);
  virtual void treeSCInfo(const std::string colPrefix, const std::string colSuffix);
  
  // Friends

  CmsSuperClusterFillerData *privateData_;

  CmsTree *cmstree;

  edm::InputTag Tracks_;
  edm::InputTag EcalBarrelRecHits_;
  edm::InputTag EcalEndcapRecHits_;

  int maxSC_;
  std::string *trkIndexName_;
};

#endif // CmsSuperClusterFiller_h
