// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HtoWWTreeDumper
// Description:
//      Class CmsGenInfoFiller
//      Simple class for dumping RECO (or AOD) contents to an ntuple
//      
// Original Author:  Alessio Ghezzi, Pietro Govoni
//         Created:  Fri Apr  6 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsGenInfoFiller_h
#define CmsGenInfoFiller_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "SimDataFormats/HepMCProduct/interface/GenInfoProduct.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"


#include "CLHEP/HepMC/GenEvent.h"
#include <TTree.h>

class CmsGenInfoFillerData;

class CmsGenInfoFiller {

 public:

  //! Constructors
  CmsGenInfoFiller(CmsTree *);

  //! Destructor
  virtual ~CmsGenInfoFiller();

  // to be used also in AODSIM 
  void
  writeGenInfoToTree (double processID, double ptHat, double genFilterEff,  double genXsec,  double weight, double AlpgenID=0);

  void 
  writeGenInfoToTree (edm::Handle<GenInfoProduct> & gi, edm::Handle<HepMCProduct>& mc, double weight );

 private:
  
  CmsGenInfoFillerData *privateData_;

};

#endif // CmsGenInfoFiller_h
