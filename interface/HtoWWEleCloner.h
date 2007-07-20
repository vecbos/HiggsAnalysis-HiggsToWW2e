#ifndef HtoWWEleCloner_h
#define HtoWWEleCloner_h

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h" 
#include "DataFormats/Common/interface/EDProduct.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include <TFile.h>
#include <TH1.h>

// class declaration
class HtoWWEleCloner : public edm::EDProducer
{
 public:
  
  explicit HtoWWEleCloner(const edm::ParameterSet&);
  
  virtual ~HtoWWEleCloner();
  
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;


 private:

  edm::InputTag electronProducer_ ;
  std::string hwwEleCandidCollection_ ;

 };
  
#endif

 
