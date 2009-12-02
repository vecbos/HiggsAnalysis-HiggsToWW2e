// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description://      Class CmsPFLowElectronFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Daniele Benedetti
//         Created:  Tue Dic 1 17:05:34 CEST 2009
//
//-----------------------------------------------------------------------

#ifndef CmsPFLowElectronFiller_h
#define CmsPFLowElectronFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include <TTree.h>

struct CmsPFLowElectronFillerData : public CmsCandidateFillerData {
  // All the vectors that will store the stuff
  // going into the tree.

  // gsf tracks
  vector<float> *gsf_pxAtInner, *gsf_pyAtInner, *gsf_pzAtInner, *gsf_xAtInner, *gsf_yAtInner, *gsf_zAtInner;
  vector<float> *gsf_pxAtOuter, *gsf_pyAtOuter, *gsf_pzAtOuter, *gsf_xAtOuter, *gsf_yAtOuter, *gsf_zAtOuter;
  vector<float> *gsf_TrackNormalizedChi2;
  vector<float> *gsf_TrackDxy, *gsf_TrackD0, *gsf_TrackDsz, *gsf_TrackDz;
  vector<float> *gsf_TrackDxyError, *gsf_TrackD0Error, *gsf_TrackDszError, *gsf_TrackDzError;
  vector<float> *gsf_TrackValidHits, *gsf_TrackLostHits;
  vector<float> *gsf_TrackVx, *gsf_TrackVy, *gsf_TrackVz; 
  vector<float> *gsf_pxAtInnerMode, *gsf_pyAtInnerMode, *gsf_pzAtInnerMode;
  vector<float> *gsf_charge,*gsf_chargeMode;

  // basic pf candidates info
  vector<float> *MvaOuput;
  vector<float> *PS1Energy,*PS2Energy;
  vector<float> *EcalEnergy,*EcalElectronEnergy;
  

public:
  void initialise();
  void clearTrkVectors();


};

class CmsPFLowElectronFiller : public CmsCandidateFiller {

public:

  //! Dump everything
  CmsPFLowElectronFiller(CmsTree *, int maxTracks=500,
			 int maxMCTracks=2000, bool noOutputIfLimitsReached=false );
  
  //! Dump  everything if fatTree is true and less informations otherwise
  CmsPFLowElectronFiller(CmsTree *, bool fatTree, int maxTracks=500,
			 int maxMCTracks=2000, bool noOutputIfLimitsReached=false );
  

  //! Destructor
  virtual ~CmsPFLowElectronFiller();
  
  void savePFEleGsfTrk(bool what){ savePFEleGsfTrk_=what;};
  void savePFEleBasic(bool what) { savePFEleBasic_=what;};


  //! write the electron related informations for the given collection
  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);
  
 
 
 

private:
  void writePFEleGsfTrkInfo(reco::GsfTrackRef gsfRef);
  void treePFEleGsfTrkInfo(const std::string &colPrefix, const std::string &colSuffix);
  void writePFEleBasicInfo(const reco::PFCandidateRef pflowCandRef);
  void treePFEleBasicInfo(const std::string &colPrefix, const std::string &colSuffix);


  bool savePFEleGsfTrk_;
  bool savePFEleBasic_;
 
  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;




  std::string *trkIndexName_;

  CmsPFLowElectronFillerData *privateData_;
  edm::InputTag matchMap_;

  CmsTree *cmstree;


};
#endif // CmsPFLowElectronFiller_h
