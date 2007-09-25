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
  vector<float> *eleTip;

public:
  void initialise();
  void clearTrkVectors();

};


class CmsEleIDTreeFiller : public CmsCandidateFiller {

 public:
  // Constructors
  CmsEleIDTreeFiller(CmsTree *, int maxTracks=500, 
		     bool noOutputIfLimitsReached=false );

  // Destructor
  virtual ~CmsEleIDTreeFiller();

  // Modifiers
  // write the electron ID starting from a Shallow-Cloned candidate collection
  void writeCollectionToTree(const CandidateCollection *,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);
  
  // write the electron ID starting from a PixelMatchGsfElectronCollection
  void writeCollectionToTree(const PixelMatchGsfElectronCollection *,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);

  // set to false if the column with the block size is set by another object
  void setStandalone(bool );

 private:

  void writeEleInfo(const PixelMatchGsfElectron *, const edm::Event&, const edm::EventSetup&);
  void treeEleInfo(const std::string &colPrefix, const std::string &colSuffix);

  int maxTracks_;
  std::string *trkIndexName_;
  bool standalone_;

  CmsTree *cmstree;
  CmsEleIDTreeFillerData *privateData_;

  const CaloGeometry* caloGeo;
};

#endif // CmsEleIDTreeFiller_h
