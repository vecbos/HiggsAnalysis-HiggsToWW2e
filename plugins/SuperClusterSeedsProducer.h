#ifndef SuperClusterSeedFilter_h
#define SuperClusterSeedFilter_h

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


class SuperClusterSeedsProducer : public edm::EDProducer{

 public:

  //! ctor
  explicit SuperClusterSeedsProducer (const edm::ParameterSet& conf) ;
  //!dtor
  ~SuperClusterSeedsProducer () {};

 private:

  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  
  virtual void beginRun(edm::Run&, edm::EventSetup const&);
  virtual void endRun(edm::Run&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  //! the configurable inputs
  edm::InputTag m_BasicClusterLabel;
  edm::InputTag m_SuperClusterLabel;

};  

#endif

