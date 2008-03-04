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

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "AnalysisDataFormats/Egamma/interface/ElectronIDAssociation.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

struct CmsEleIDTreeFillerData : public CmsCandidateFillerData {

  vector<int>   *eleClass;
  vector<float> *eleHoE;
  vector<float> *eleNotCorrEoP,    *eleCorrEoP;
  vector<float> *eleNotCorrEoPout, *eleCorrEoPout;
  vector<float> *eleDeltaEtaAtVtx, *eleDeltaEtaAtCalo;
  vector<float> *eleDeltaPhiAtVtx, *eleDeltaPhiAtCalo;
  vector<float> *eleTrackerP;
  vector<float> *eleTrackerIso_minDR, *eleTrackerIso_sumPt;
  vector<float> *eleTrackerIso_minDR_veto;
  vector<float> *eleCaloIso_minDR,    *eleCaloIso_sumPt;
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
  
  //! set the cluster shape association map for ECAL barrel
  void setEcalBarrelClusterShapes( edm::InputTag EcalBarrelClusterShapes ) { EcalBarrelClusterShapes_ = EcalBarrelClusterShapes; }
  //! set the cluster shape association map for ECAL endcap
  void setEcalEndcapClusterShapes( edm::InputTag EcalEndcapClusterShapes ) { EcalEndcapClusterShapes_ = EcalEndcapClusterShapes; }
  //! set the electron ID association map
  void setElectronIdAssociation( edm::InputTag electronIDAssocProducer ) { electronIDAssocProducer_ = electronIDAssocProducer; }
  //! set to false if the column with the block size is set by another object
  void setStandalone(bool );

 private:

  void writeEleInfo(const PixelMatchGsfElectron *, int index, const edm::Event&, const edm::EventSetup&,
		    const reco::BasicClusterShapeAssociationCollection& barrelClShpMap, 
		    const reco::BasicClusterShapeAssociationCollection& endcapClShpMap,
		    const reco::ElectronIDAssociationCollection& eleIdAssoc);
  void treeEleInfo(const std::string &colPrefix, const std::string &colSuffix);

  int maxTracks_;
  std::string *trkIndexName_;
  bool standalone_;

  edm::InputTag EcalBarrelClusterShapes_;
  edm::InputTag EcalEndcapClusterShapes_;
  
  edm::InputTag electronIDAssocProducer_;

  CmsTree *cmstree;
  CmsEleIDTreeFillerData *privateData_;

  edm::Handle<reco::PixelMatchGsfElectronCollection> explicitElectronCollectionHandle_;

  const CaloGeometry* caloGeo;
};

#endif // CmsEleIDTreeFiller_h
