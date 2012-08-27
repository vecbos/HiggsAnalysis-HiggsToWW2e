#ifndef ElectronAndPhotonSuperClusterProducer_h
#define ElectronAndPhotonSuperClusterProducer_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include <functional>
#include <vector>
#include <map>

// N.B. the reason to make a producer and not a filter is that 
// we cannot merge two RefVectors coming from 2 different collections into one (EB + EE)
class ElectronAndPhotonSuperClusterProducer : public edm::EDProducer{

 public:

  //! ctor
  explicit ElectronAndPhotonSuperClusterProducer (const edm::ParameterSet& conf) ;
  //!dtor
  ~ElectronAndPhotonSuperClusterProducer () {};

 private:

  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  
  virtual void beginRun(edm::Run&, edm::EventSetup const&);
  virtual void endRun(edm::Run&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  //! the configurable inputs
  edm::InputTag m_ElectronLabel;
  edm::InputTag m_PhotonLabel;
  edm::InputTag m_SuperClusterLabel;
  bool m_includePhotonSuperClusters;

};  

#endif

