// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsPhotonFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Emanuele Di Marco
//         Created:  Fri Apr  6 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsPhotonFiller_h
#define CmsPhotonFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"                                                                             
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"             
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include <TTree.h>

struct CmsPhotonFillerData : public CmsCandidateFillerData {

  // All the vectors that will store the stuff
  // going into the tree.
  vector<int> *fiducialFlags, *recoFlags;
  vector<int> *superClusterIndex, *PFsuperClusterIndex;
  vector<float> *hOverE, *hTowOverE, *dr03TkSumPt, *dr03HollowTkSumPt, *dr03EcalRecHitSumEt, *dr03HcalTowerSumEt,
    *dr04TkSumPt, *dr04HollowTkSumPt, *dr04EcalRecHitSumEt, *dr04HcalTowerSumEt;
  vector<float> *chargedHadronIso, *neutralHadronIso, *photonIso;
  vector<int> *hasPixelSeed;
  vector<bool> *hasMatchedConversion;

  vector<float > *dr01chPFIso, *dr01nhPFIso, *dr01phPFIso;
  vector<float > *dr02chPFIso, *dr02nhPFIso, *dr02phPFIso;
  vector<float > *dr03chPFIso, *dr03nhPFIso, *dr03phPFIso;
  vector<float > *dr04chPFIso, *dr04nhPFIso, *dr04phPFIso;
  vector<float > *dr05chPFIso, *dr05nhPFIso, *dr05phPFIso;
  vector<float > *dr06chPFIso, *dr06nhPFIso, *dr06phPFIso;

  void initialise();
  void clearTrkVectors();


};

class CmsPhotonFiller : public CmsCandidateFiller {

 public:

  //! Dump everything
  CmsPhotonFiller(CmsTree *, int maxTracks=500,
		int maxMCTracks=2000, bool noOutputIfLimitsReached=false );

  //! Destructor
  virtual ~CmsPhotonFiller();

  //! write the electron related informations for the given collection
  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);

  //! set the merged ECAL superclusters collection
  void setEcalSuperClusters( edm::InputTag EcalSuperClusters ) { EcalSuperClusters_ = EcalSuperClusters; }
  //! set conversions collection
  void setConversionsProdcer(  edm::InputTag conversionsProducer ) { conversionsProducer_ = conversionsProducer; }

  void setPFCandidates( edm::InputTag pfCands) { pfCandidates_ = pfCands; }
  void setPrimaryVertices( edm::InputTag pv) {primaryVertices_ = pv; }
 private:
  
  void writeEcalInfo(const reco::PhotonRef, const edm::Event&, const edm::EventSetup&, 
                     reco::SuperClusterRef sclusRef, reco::SuperClusterRef pfclusRef);
  void treeEcalInfo(const std::string &colPrefix, const std::string &colSuffix);


  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  std::string *trkIndexName_;

  CmsPhotonFillerData *privateData_;

  edm::InputTag EcalSuperClusters_;
  edm::InputTag conversionsProducer_;
  edm::InputTag pfCandidates_;
  edm::InputTag primaryVertices_;

  edm::Handle<reco::ConversionCollection> hConversions;
  edm::Handle<reco::BeamSpot> bsHandle;
  edm::Handle<reco::PFCandidateCollection> pfCands;
  edm::Handle<reco::VertexCollection> PVs;
  edm::Handle<reco::SuperClusterCollection> hSuperClusters;

  CmsTree *cmstree;


};

#endif // CmsPhotonFiller_h
