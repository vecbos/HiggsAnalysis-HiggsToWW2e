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
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"

#include <TTree.h>

struct CmsGenInfoFillerData {
  
  float *ptHat, *processID, *weight;
  float *alphaQCD, *alphaQED;

public:
  void initialise();

};

class CmsGenInfoFiller {

 public:

  //! Constructors
  CmsGenInfoFiller(CmsTree *);

  //! Destructor
  virtual ~CmsGenInfoFiller();

  void writeGenInfoToTree ( edm::Handle<GenEventInfoProduct> & gei);

 private:
  
  CmsTree *cmstree;
  CmsGenInfoFillerData *privateData_;

};

#endif // CmsGenInfoFiller_h
