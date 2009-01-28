// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsEleIDTreeFiller
//      Simple class for dumping electron ID contents to a ROOT tree
//      
// Original Authors:  Chiara Ilaria Rovelli, Emanuele Di Marco
//          Created:  Fri May  18 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsEleIDTreeFiller_h
#define CmsEleIDTreeFiller_h

#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"

struct CmsEleIDTreeFillerData : public CmsCandidateFillerData {

  vector<int>   *eleClass;
  vector<float> *eleHoE;
  vector<float> *eleNotCorrEoP,    *eleCorrEoP;
  vector<float> *eleNotCorrEoPout, *eleCorrEoPout;
  vector<float> *eleDeltaEtaAtVtx, *eleDeltaEtaAtCalo;
  vector<float> *eleDeltaPhiAtVtx, *eleDeltaPhiAtCalo;
  vector<float> *eleTrackerP;
  vector<float> *minDR03, *minDRveto03, *sumPt03, *sumPtSquared03, *sumN03;
  vector<float> *sumPt04, *sumPt05;
  vector<float> *sumPtPreselection;
  vector<float> *sumHadEt04, *sumEmEt04, *sumHadEt05, *sumEmEt05;
  vector<float> *isoFromDepsTk, *isoFromDepsEcal, *isoFromDepsHcal;
  vector<float> *scBasedEcalSum04, *scBasedEcalSum05, *scHaloBasedEcalSum04, *scHaloBasedEcalSum05;
  vector<float> *eleFullCorrE,     *eleCaloCorrE;
  vector<float> *eleNxtalCorrE,    *eleRawE;
  vector<float> *eleLik;
  vector<bool> *eleIdCutBasedDecision;
  vector<float> *eleTip;

public:
  void initialise();
  void clearTrkVectors();

};


class CmsEleIDTreeFiller : public CmsCandidateFiller {

 public:

  //! Constructor
  CmsEleIDTreeFiller(CmsTree *, int maxTracks=500, 
		     bool noOutputIfLimitsReached=false );

  //! Destructor
  virtual ~CmsEleIDTreeFiller();

  //! write the electron ID starting from a candidate collection
  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);
  
  //! set the rechits for ECAL barrel (needed for cluster shapes)
  void setEcalBarrelRecHits( edm::InputTag EcalBarrelRecHits ) { EcalBarrelRecHits_ = EcalBarrelRecHits; }
  //! set the rechits for ECAL endcap (needed for cluster shapes)
  void setEcalEndcapRecHits( edm::InputTag EcalEndcapRecHits ) { EcalEndcapRecHits_ = EcalEndcapRecHits; }
  //! set the electron ID labels
  void setElectronIdCutsLabel( edm::InputTag electronIdCutsLabel ) { electronIdCutsLabel_ = electronIdCutsLabel; }
  void setElectronIdLikelihoodLabel( edm::InputTag electronIdLikelihoodLabel ) { electronIdLikelihoodLabel_ = electronIdLikelihoodLabel; }
  //! set the tracker isolation producer (egamma)
  void setTkIsolationProducer( edm::InputTag tkIsolationProducer ) { tkIsolationProducer_ = tkIsolationProducer; }
  //! set the HCAL isolation producer with calo towers (egamma)
  void setTowerIsolationProducer( edm::InputTag towerIsolationProducer ) { towerIsolationProducer_ = towerIsolationProducer; }
  //! set the track producer for tracker isolation
  void setTracksProducer( edm::InputTag tracksProducer ) { tracksProducer_ = tracksProducer; }
  //! set the calotower producer for calorimetric isolation
  void setCalotowersProducer( edm::InputTag calotowersProducer ) { calotowersProducer_ = calotowersProducer; }
  
  //! set to false if the column with the block size is set by another object
  void setStandalone(bool );

 private:

  void writeEleInfo(const GsfElectronRef, const edm::Event&, const edm::EventSetup&,
		    const EcalRecHitCollection *EBRecHits,
		    const EcalRecHitCollection *EERecHits,
		    const edm::ValueMap<float> & eIdmapCuts, const edm::ValueMap<float> & eIdmapLikelihood);
  void treeEleInfo(const std::string &colPrefix, const std::string &colSuffix);

  int maxTracks_;
  std::string *trkIndexName_;
  bool standalone_;

  edm::InputTag EcalBarrelRecHits_;
  edm::InputTag EcalEndcapRecHits_;
  
  edm::InputTag electronIdCutsLabel_;
  edm::InputTag electronIdLikelihoodLabel_;
  edm::InputTag tkIsolationProducer_;
  edm::InputTag towerIsolationProducer_;

  edm::InputTag tracksProducer_;
  edm::InputTag calotowersProducer_;

  CmsTree *cmstree;
  CmsEleIDTreeFillerData *privateData_;

  typedef edm::ValueMap<double> isoFromDepositsMap;
  typedef std::vector< edm::Handle<isoFromDepositsMap> > isoContainer;
  isoContainer *eIsoFromDepsValueMap_;
//   edm::Handle< reco::CandViewDoubleAssociations > towerIsolationHandle_;
//   edm::Handle< reco::PMGsfElectronIsoCollection > tkIsolationHandle_;
  edm::Handle<CaloTowerCollection> m_calotowers;
  edm::Handle<reco::TrackCollection> m_tracks;

  const CaloGeometry* caloGeo;
};

#endif // CmsEleIDTreeFiller_h
