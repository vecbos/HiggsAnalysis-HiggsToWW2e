// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description://      Class CmsPFlowElectronFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Daniele Benedetti
//         Created:  Tue Dic 1 17:05:34 CEST 2009
//
//-----------------------------------------------------------------------

#ifndef CmsPFlowElectronFiller_h
#define CmsPFlowElectronFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include <TTree.h>

struct CmsPFlowElectronFillerData : public CmsCandidateFillerData {
  // All the vectors that will store the stuff
  // going into the tree.

  // tracks
  vector<int> *trackIndex, *gsfTrackIndex;
  vector<int> *convTrackIndex;
  vector<float> *convDist, *convDcot, *convRadius;
  
  // id a la egamma 
  vector<float> *deltaEtaAtVtxEg;
  vector<float> *deltaPhiAtVtxEg;

  // basic pf candidates info
  vector<float> *MvaOutput;
  vector<float> *PS1Energy,*PS2Energy;
  vector<float> *EcalEnergy,*RawEcalEnergy;  
  vector<float> *HcalEnergy,*RawHcalEnergy;
  vector<float> *PositionAtEcalX;
  vector<float> *PositionAtEcalY;
  vector<float> *PositionAtEcalZ;
  vector<float> *chIso03veto,      *chIso04veto,      *chIso05veto;
  vector<float> *chIso03noVeto,    *chIso04noVeto,    *chIso05noVeto;
  vector<float> *nhIso03veto,      *nhIso04veto,      *nhIso05veto;
  vector<float> *nhIso03noVeto,    *nhIso04noVeto,    *nhIso05noVeto;
  vector<float> *phIso03veto,      *phIso04veto,      *phIso05veto;
  vector<float> *phIso03noVeto,    *phIso04noVeto,    *phIso05noVeto;
  vector<float> *chIso03vetoNVC,   *chIso04vetoNVC,   *chIso05vetoNVC;
  vector<float> *chIso03noVetoNVC, *chIso04noVetoNVC, *chIso05noVetoNVC;
  vector<float> *nhIso03vetoTJ,    *nhIso04vetoTJ,    *nhIso05vetoTJ;
  vector<float> *phIso03vetoTJ,    *phIso04vetoTJ,    *phIso05vetoTJ;

public:
  void initialise();
  void clearTrkVectors();


};

class CmsPFlowElectronFiller : public CmsCandidateFiller {

public:

  //! Dump everything
  CmsPFlowElectronFiller(CmsTree *, int maxTracks=500,
			 int maxMCTracks=2000, bool noOutputIfLimitsReached=false );
  
  //! Dump  everything if fatTree is true and less informations otherwise
  CmsPFlowElectronFiller(CmsTree *, bool fatTree, int maxTracks=500,
			 int maxMCTracks=2000, bool noOutputIfLimitsReached=false );
  

  //! Destructor
  virtual ~CmsPFlowElectronFiller();
  
  void savePFEleBasic(bool what)  { savePFEleBasic_=what;};
  void savePFEleIsoDep(bool what) { savePFEleIsoDep_=what;};

  //! write the electron related informations for the given collection
  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);
  

  //! set the general tracks 
  void setGeneralTracks( edm::InputTag generalTracks) { generalTracks_ = generalTracks; }

 
 
 

private:
  void writePFEleTrkInfo(const reco::PFCandidateRef, reco::GsfTrackRef gsfRef, reco::TrackRef kfTrackRef, const edm::EventSetup& iSetup );
  void treePFEleTrkInfo(const std::string &colPrefix, const std::string &colSuffix);

  void writePFEleBasicInfo(const reco::PFCandidateRef pflowCandRef);
  void treePFEleBasicInfo(const std::string &colPrefix, const std::string &colSuffix);

  void treePFEleIsoInfo(const std::string &colPrefix, const std::string &colSuffix);

  void writePFEleEgammaIDInfo(const reco::PFCandidateRef pflowCandRef);
  void treePFEleEgammaIDInfo(const std::string &colPrefix, const std::string &colSuffix);

  void getConversionInfo(reco::PFCandidateRef pflowCandRef,
                         const edm::Handle<reco::TrackCollection>& track_h,
                         const double bFieldAtOrigin);

  bool savePFEleBasic_;
  bool savePFEleIsoDep_;

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  float theCFDist_, theCFDcot_, theCFRconv_;
  int theCFTrackIndex_;
  edm::InputTag generalTracks_;
  edm::Handle< reco::TrackCollection > h_tracks;

  std::string *trkIndexName_;

  CmsPFlowElectronFillerData *privateData_;

  CmsTree *cmstree;


};
#endif // CmsPFlowElectronFiller_h
