// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CsmTrackFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
//-----------------------------------------------------------------------

#ifndef CmsTrackFiller_h
#define CmsTrackFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "boost/mpl/vector.hpp" 
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include <TTree.h>

struct CmsTrackFillerData {

  int *ncand;
  vector<int> *vtxIndex;
  vector<float> *vtxWeight;

  vector<float> *pxAtInner, *pyAtInner, *pzAtInner, *xAtInner, *yAtInner, *zAtInner;
  vector<float> *pxAtOuter, *pyAtOuter, *pzAtOuter, *xAtOuter, *yAtOuter, *zAtOuter;
  vector<float> *px, *py, *pz;
  vector<float> *trackNormalizedChi2;
  vector<int> *qualityMask;
  vector<float> *trackValidHits, *trackLostHits;
  vector<float> *trackVx, *trackVy, *trackVz; 

  vector<float> *charge, *pterr, *recHitsSize;
  vector<float> *impactPar3D, *impactPar3DError, *transvImpactPar, *transvImpactParError;
  vector<float> *impactPar3DBiased, *impactPar3DBiasedError, *transvImpactParBiased, *transvImpactParBiasedError;

  vector<float> *truncatedDeDx, *truncatedDeDxError, *truncatedDeDxNoM;
  vector<float> *medianDeDx, *medianDeDxError, *medianDeDxNoM;
  vector<float> *harmonic2DeDx, *harmonic2DeDxError, *harmonic2DeDxNoM;
  vector<int> *pixelHits, *expInnerLayers, *trackerLayersWithMeasurement;
  vector<float> *d0,*d0Error,*dz,*dzError;
  vector<int> *numberOfValidPixelBarrelHits, *numberOfValidPixelEndcapHits;
  vector<int> *numberOfValidStripTIBHits, *numberOfValidStripTIDHits, *numberOfValidStripTOBHits, *numberOfValidStripTECHits;
  vector<int> *numberOfValidMuonHits;

public:
  void initialise();
  void clearTrkVectors();

};

class CmsTrackFiller {

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

  //! set the tracker isolation producer
  void setRefittedTracksForDeDxProducer( edm::InputTag tracksTag ) { refittedTracksForDeDxTag_ = tracksTag; }

  /// dump tracking related variables
  void saveTrk(bool );
  /// are tracks Gsf
  void isGsf(bool );
  /// dump track-extras related variables
  void saveFatTrk(bool );
  /// dump Vtx related variables
  void saveVtxTrk(bool );
  /// Find Primary Vertex
  void findHardestPrimaryVertex(const edm::Event& iEvent);
  /// Save the dEdX: needs the right sequence to be run in the cms path
  void saveDeDx(bool );

  // Operators

  /// write the muon related informations for the given collection
  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);

  

 protected:
  
  void writeTrkInfo(edm::RefToBase<reco::Track> trackRef, const edm::ESHandle<MagneticField>& magfield, const edm::ESHandle<GlobalTrackingGeometry>& theTrackingGeometry);
  void writeDeDxInfo(edm::RefToBase<reco::Track> refittedTrack);
  void treeTrkInfo(const std::string &colPrefix, const std::string &colSuffix);
  void treeDeDxInfo(const std::string &colPrefix, const std::string &colSuffix);
  bool hasValidHitInNthPixelBarrel(uint32_t nlayer, reco::HitPattern pattern);
  bool hasValidHitInNthPixelEndcap(uint32_t nlayer, reco::HitPattern pattern);

  bool saveTrk_;
  bool isGsf_;
  bool saveFatTrk_;
  bool saveVtxTrk_;
  bool saveDeDx_;

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  std::string *trkIndexName_;

  CmsTrackFillerData *privateData_;

  CmsTree *cmstree;

  edm::InputTag vertexCollection_;
  edm::InputTag refittedTracksForDeDxTag_;

  edm::Handle< edm::View<reco::Track> > refittedTracksForDeDx_;
  edm::Handle< reco::DeDxDataValueMap >  truncatedEnergyLoss_,  medianEnergyLoss_, harmonic2EnergyLoss_;
  edm::Handle<reco::VertexCollection> primaryVertex_;
  edm::ESHandle<TransientTrackBuilder> trackBuilder_;

  // number of 32 bit integers to store the full pattern
  const static unsigned short PatternSize = 25;

  // number of bits used for each hit
  const static unsigned short HitSize = 11;    

  // Primary Vertex in point format
  reco::Vertex bestPrimaryVertex_;
  float x0, y0, z0;
};

#endif // CmsTrackFiller_h
