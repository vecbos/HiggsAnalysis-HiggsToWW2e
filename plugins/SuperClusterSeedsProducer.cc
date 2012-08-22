#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include "HiggsAnalysis/HiggsToWW2e/plugins/SuperClusterSeedsProducer.h"

#include <iostream>
#include <algorithm>

using namespace reco;

SuperClusterSeedsProducer::SuperClusterSeedsProducer (const edm::ParameterSet& conf)
{
  m_BasicClusterLabel = conf.getParameter<edm::InputTag>("BasicClusterLabel");
  m_SuperClusterLabel = conf.getParameter<edm::InputTag>("SuperClusterLabel");
  produces<BasicClusterCollection>();
}


// ------------------------------------------------------------


void 
SuperClusterSeedsProducer::produce (edm::Event& iEvent,  const edm::EventSetup& evtSetup ) 
{

  std::auto_ptr<BasicClusterCollection> selected( new BasicClusterCollection );  

  // get the BasicClusterCollection
  edm::Handle<reco::BasicClusterCollection> basicclusters;
  try { iEvent.getByLabel(m_BasicClusterLabel, basicclusters); }
  catch ( cms::Exception& ex ) { edm::LogWarning("SuperClusterSeedsProducer") << "Can't get BC Collection: " << m_BasicClusterLabel; }

  // get the SuperClusterCollection
  edm::Handle<reco::SuperClusterCollection> superclusters;
  try { iEvent.getByLabel(m_SuperClusterLabel, superclusters); }
  catch ( cms::Exception& ex ) { edm::LogWarning("SuperClusterSeedsProducer") << "Can't get SC Collection: " << m_SuperClusterLabel; }

  edm::LogInfo("SuperClusterSeedsProducer") << "SuperClusterSeedsProducer starting, original collection size = " 
                                            << basicclusters->size();

  reco::BasicClusterCollection::const_iterator it = basicclusters->begin();
  for(unsigned int ibclu = 0; ibclu < basicclusters->size(); ibclu++) {

    reco::BasicClusterRef bcRef(basicclusters,ibclu);
    double bcEnergy = bcRef->energy();
    const math::XYZPoint &bcPosition = bcRef->position();

    // look if it is a seed of a supercluster
    for(unsigned int isclu = 0; isclu < superclusters->size(); isclu++) {
      const reco::SuperClusterRef scRef(superclusters,isclu);
      edm::Ptr<reco::CaloCluster> theSeed = scRef->seed();

      double seedEnergy = theSeed->energy();
      const math::XYZPoint &seedPosition = theSeed->position();
      
      if(bcEnergy==seedEnergy && bcPosition==seedPosition) {
        selected->push_back(*it);
      }
    }
    it++;
  }

  edm::LogInfo("SuperClusterSeedsProducer") << "SuperClusterSeedsProducer ending, final collection size = " 
                                            << selected->size();
  
  iEvent.put( selected );

}       

// ------------ method called once each job just before starting event loop  ------------
void 
SuperClusterSeedsProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SuperClusterSeedsProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
SuperClusterSeedsProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
SuperClusterSeedsProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
SuperClusterSeedsProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
SuperClusterSeedsProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}
