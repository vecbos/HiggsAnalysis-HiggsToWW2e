// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsGsfTrackFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
//-----------------------------------------------------------------------

#ifndef CmsGsfTrackFiller_h
#define CmsGsfTrackFiller_h

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTrackFiller.h"

struct CmsGsfTrackFillerData : public CmsTrackFillerData {

  vector<int> *chargeMode;
  vector<float> *pxMode, *pyMode, *pzMode;

public:
  void initialiseGsf();
  void clearTrkVectorsGsf();

};

class CmsGsfTrackFiller : public CmsTrackFiller {

 public:

  // Constructors

  // Dump everything
  CmsGsfTrackFiller(CmsTree *, 
                    edm::InputTag vertexCollection,
                    int maxTracks=500, int maxMCTracks=2000, 
                    bool noOutputIfLimitsReached=false );
  
  // Dump  everything if fatTree is true and less informations otherwise
  CmsGsfTrackFiller(CmsTree *, 
                    edm::InputTag vertexCollection,
                    bool fatTree, 
                    int maxTracks=500, int maxMCTracks=2000, 
                    bool noOutputIfLimitsReached=false, bool vtxtrack=false);
  
  // Destructor
  virtual ~CmsGsfTrackFiller();

  /// write the muon related informations for the given collection
  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);

 private:
  
  void writeGsfTrkInfo(reco::GsfTrackRef trackRef);
  void treeGsfTrkInfo(const std::string &colPrefix, const std::string &colSuffix);

  CmsGsfTrackFillerData *privateData_;

  CmsTree *cmstree;
  
};

#endif // CmsGsfTrackFiller_h

