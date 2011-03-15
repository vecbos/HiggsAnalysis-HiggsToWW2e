// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsElectronFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Emanuele Di Marco
//         Created:  Fri Apr  6 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsHcalNoiseFiller_h
#define CmsHcalNoiseFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "TrackingTools/TrackAssociator/interface/CachedTrajectory.h"
#include "TrackingTools/TrackAssociator/interface/CaloDetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/EcalDetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/MuonDetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/HcalDetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/HODetIdAssociator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include <TTree.h>

struct CmsHcalNoiseFillerData
{
public:
  void initialise() {}
};

class CmsHcalNoiseFiller
{

 public:

  // Constructors

  // Dump everything
  CmsHcalNoiseFiller(CmsTree *);

  // Dump  everything if fatTree is true and less informations otherwise
  CmsHcalNoiseFiller(CmsTree *, bool fatTree);
  
  // Destructor
  virtual ~CmsHcalNoiseFiller();

  // Operators

  //! write the hcal noise related informations for the given collection
  void writeHcalNoiseSummaryToTree(edm::InputTag rechitHBHETag,
              edm::InputTag rechitHFTag,
              edm::InputTag noiseSummaryTag,
			     const edm::Event &iEvent, const edm::EventSetup &iSetup);

  
 private:
  
  CmsHcalNoiseFillerData *privateData_;

  CmsTree *cmstree;

};

#endif // CmsHcalNoiseFiller_h
