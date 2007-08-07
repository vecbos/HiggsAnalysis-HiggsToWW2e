// -*- C++ -*-
//
// Package:    HtoWWElectrons
// Class:      hwwEleCaloIsolation
// 
/*
   Description: <one line class summary>
   Electron isolation using tracker info

   Implementation:
 
*/
//
// Original Author:  Chiara Rovelli
//
//


#ifndef hwwEleCaloIsolation_h
#define hwwEleCaloIsolation_h

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

using namespace edm;
using namespace std;
using namespace reco;

class hwwEleCaloIsolation{
 public:
  
  //constructors
  hwwEleCaloIsolation();
  hwwEleCaloIsolation(const PixelMatchGsfElectron *gsfEle, const HBHERecHitCollection trackColl, edm::ESHandle<CaloGeometry>);
  
  //methods
  void setExtRadius (float extRadius);
  
  float getEtHcal () const;
  bool isIsolated (float etCut = 0.05) const;

  //destructor 
  ~hwwEleCaloIsolation();
  
 private:

  const PixelMatchGsfElectron* _myGsfEle;  	     
  const HBHERecHitCollection _hcalrh;
  edm::ESHandle<CaloGeometry> _caloGeo;
  
  float _extRadius;
};

#endif
