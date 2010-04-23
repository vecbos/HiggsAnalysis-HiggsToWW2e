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

  // basic pf candidates info
  vector<float> *MvaOutput;
  vector<float> *PS1Energy,*PS2Energy;
  vector<float> *EcalEnergy,*RawEcalEnergy;  
  vector<float> *HcalEnergy,*RawHcalEnergy;
  vector<float> *PositionAtEcalX;
  vector<float> *PositionAtEcalY;
  vector<float> *PositionAtEcalZ;

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
  
  void savePFEleTrk(bool what)   { savePFEleTrk_=what;};
  void savePFEleBasic(bool what) { savePFEleBasic_=what;};

  //! write the electron related informations for the given collection
  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);
  
 
 
 

private:
  void writePFEleTrkInfo(reco::GsfTrackRef gsfRef,reco::TrackRef kfTrackRef);
  void treePFEleTrkInfo(const std::string &colPrefix, const std::string &colSuffix);
  void writePFEleBasicInfo(const reco::PFCandidateRef pflowCandRef);
  void treePFEleBasicInfo(const std::string &colPrefix, const std::string &colSuffix);

  bool savePFEleTrk_;
  bool savePFEleBasic_;
  
  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;


  std::string *trkIndexName_;

  CmsPFlowElectronFillerData *privateData_;
  edm::InputTag matchMap_;

  CmsTree *cmstree;


};
#endif // CmsPFlowElectronFiller_h
