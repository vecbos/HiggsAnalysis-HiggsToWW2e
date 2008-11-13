// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      ZjetsAnalysis/ZllProducer
// Description:
//      Class CsmTrackFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
//-----------------------------------------------------------------------

#ifndef CmsTrackFiller_h
#define CmsTrackFiller_h

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "DataFormats/Math/interface/Point3D.h"
#include <TTree.h>

using namespace cms;
using namespace edm;
using namespace reco;

struct CmsTrackFillerData : public CmsCandidateFillerData {

  vector<int> *vtxIndex;
  vector<float> *vtxWeight;

  vector<float> *pxAtInner, *pyAtInner, *pzAtInner, *xAtInner, *yAtInner, *zAtInner;
  vector<float> *pxAtOuter, *pyAtOuter, *pzAtOuter, *xAtOuter, *yAtOuter, *zAtOuter;
  vector<float> *trackNormalizedChi2;
  vector<float> *trackDxy, *trackD0, *trackDsz, *trackDz;
  vector<float> *trackDxyError, *trackD0Error, *trackDszError, *trackDzError;
  vector<float> *trackValidHits, *trackLostHits;
  vector<float> *trackVx, *trackVy, *trackVz; 

  vector<float> *charge, *pterr, *recHitsSize;
  vector<float> *trackDxyPV, *trackD0PV, *trackDszPV, *trackDzPV;
  vector<float> *trackDxyErrorPV, *trackD0ErrorPV, *trackDszErrorPV, *trackDzErrorPV;

public:
  void initialise();
  void clearTrkVectors();

};

class CmsTrackFiller : public CmsCandidateFiller {

 public:

  // Constructors

  // Dump everything
  CmsTrackFiller(CmsTree *, 
		 edm::InputTag vertexCollection,
		 int maxTracks=500, int maxMCTracks=2000, 
		 bool noOutputIfLimitsReached=false );
  
  // Dump  everything if fatTree is true and less informations otherwise
  CmsTrackFiller(CmsTree *, 
		 edm::InputTag vertexCollection,
		 bool fatTree, 
		 int maxTracks=500, int maxMCTracks=2000, 
		 bool noOutputIfLimitsReached=false, bool vtxtrack=false);
  
  // Destructor
  virtual ~CmsTrackFiller();

  /// dump tracking related variables
  void saveTrk(bool );
  /// dump track-extras related variables
  void saveFatTrk(bool );
  /// dump Vtx related variables
  void saveVtxTrk(bool );
  /// Find Primary Vertex
  void findPrimaryVertex(const edm::Event& iEvent);

  // Operators

  /// write the muon related informations for the given collection
  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);

  

 private:
  
  void writeTrkInfo(const Candidate *cand, const edm::Event&, const edm::EventSetup&, TrackRef trkRef);
  void treeTrkInfo(const std::string &colPrefix, const std::string &colSuffix);

  bool saveTrk_;
  bool saveFatTrk_;
  bool saveVtxTrk_;

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  std::string *trkIndexName_;

  CmsTrackFillerData *privateData_;
  edm::InputTag matchMap_;

  CmsTree *cmstree;

  edm::InputTag vertexCollection_;

  // Primary Vertex in point format
  float x0, y0, z0;
};

#endif // CmsTrackFiller_h
