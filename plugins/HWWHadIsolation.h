#ifndef HWWHadIsolation_h
#define HWWHadIsolation_h

/** \class HWWHadIsolation
 *
 * Object selector perform electron track isolation selection
 *
 */  
 
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/EventPrincipal.h" 
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/Common/interface/Handle.h" 
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"

class HWWHadIsolation {

 public:
  
  HWWHadIsolation(const edm::ParameterSet& conf);
  
  ~HWWHadIsolation();
  
  typedef reco::GsfElectronCollection collection;
  typedef std::vector<const reco::GsfElectron *> container;
  typedef container::const_iterator const_iterator;

  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
 
  void select( edm::Handle<reco::GsfElectronCollection>, 
               const edm::Event&, const edm::EventSetup& );

 private:
 
  container selected_;
  
  std::string hcalrhitsLabel_;
  double radiusConeExt_;
  double radiusConeInt_;
  double eTMin_;
  double cut_;
  
  edm::ESHandle<CaloGeometry> theCaloGeom_;

};
  
#endif
 


