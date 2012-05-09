// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CsmConversionFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
//-----------------------------------------------------------------------

#ifndef CmsConversionFiller_h
#define CmsConversionFiller_h

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

struct CmsConversionFillerData {

  int *ncand;
  vector<float> *pxPair, *pyPair, *pzPair;
  vector<float> *pxRefittedPair, *pyRefittedPair, *pzRefittedPair;
  vector<float> *etaRefittedPair, *phiRefittedPair, *ptRefittedPair, *energyRefittedPair;
  vector<float> *eOverPRefittedPair,*zOfPVFromTracks;

  vector<float> *xVtx,*yVtx,*zVtx,*chi2Vtx,*chi2ProbVtx;
  vector<float> *mvaOutVtx;
  vector<int>   *nTracksVtx,*isValidVtx;
  vector<float>   *trk1Dz,*trk1DzError,*trk1Charge,*trk1Algo,*trk1D0,*trk1Pout,*trk1Pin,*trk1Pt;
  vector<float>   *trk2Dz,*trk2DzError,*trk2Charge,*trk2Algo,*trk2D0,*trk2Pout,*trk2Pin,*trk2Pt;

public:
  void initialise();
  void clearConvVectors();

};

class CmsConversionFiller {

 public:

  // Constructors

  // Dump everything
  CmsConversionFiller(CmsTree *, 
		      //edm::InputTag conversionCollection,
		      int maxConv=500, int maxMCConv=2000, 
		      bool noOutputIfLimitsReached=false );
  
  // Destructor
  virtual ~CmsConversionFiller();

  /// dump tracking related variables
  // Operators

  /// write the muon related informations for the given collection
  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);

  

 protected:
  
  void writeConvInfo(edm::RefToBase<reco::Conversion> trackRef, const edm::ESHandle<MagneticField>& magfield, const edm::ESHandle<GlobalTrackingGeometry>& theTrackingGeometry);
  void treeConvInfo(const std::string &colPrefix, const std::string &colSuffix);

  bool hitLimitsMeansNoOutput_;
  int maxConv_;
  int maxMCConv_;

  std::string *convIndexName_;

  CmsConversionFillerData *privateData_;

  CmsTree *cmstree;

  edm::InputTag conversionCollection_;

  // number of 32 bit integers to store the full pattern
  const static unsigned short PatternSize = 25;

  // number of bits used for each hit
  const static unsigned short HitSize = 11;    

};

#endif // CmsConversionFiller_h
