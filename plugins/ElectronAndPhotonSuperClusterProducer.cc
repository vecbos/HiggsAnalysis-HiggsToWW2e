#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "HiggsAnalysis/HiggsToWW2e/plugins/ElectronAndPhotonSuperClusterProducer.h"

#include <iostream>
#include <algorithm>

using namespace reco;
using namespace std;

ElectronAndPhotonSuperClusterProducer::ElectronAndPhotonSuperClusterProducer (const edm::ParameterSet& conf)
{
  m_ElectronLabel = conf.getParameter<edm::InputTag>("ElectronLabel");
  m_PhotonLabel = conf.getParameter<edm::InputTag>("PhotonLabel");
  m_SuperClusterLabel = conf.getParameter<edm::InputTag>("SuperClusterLabel");
  m_includePhotonSuperClusters = conf.getParameter<bool>("includePhotonSuperClusters");
  produces<SuperClusterCollection>();
}


// ------------------------------------------------------------


void 
ElectronAndPhotonSuperClusterProducer::produce (edm::Event& iEvent,  const edm::EventSetup& evtSetup ) 
{

  std::auto_ptr<SuperClusterCollection> selected( new SuperClusterCollection );  

  // get the original SuperClusterCollection
  edm::Handle<SuperClusterCollection> superclusters;
  try { iEvent.getByLabel(m_SuperClusterLabel, superclusters); }
  catch ( cms::Exception& ex ) { edm::LogWarning("ElectronAndPhotonSuperClusterProducer") << "Can't get SC Collection: " << m_SuperClusterLabel; }

  // get the ElectronCollection
  edm::Handle<GsfElectronCollection> electrons;
  try { iEvent.getByLabel(m_ElectronLabel, electrons); }
  catch ( cms::Exception& ex ) { edm::LogWarning("ElectronAndPhotonSuperClusterProducer") << "Can't get electron Collection: " << m_ElectronLabel; }

  // get the PhotonCollection
  edm::Handle<PhotonCollection> photons;
  if(m_includePhotonSuperClusters) {
    try { iEvent.getByLabel(m_PhotonLabel, photons); }
    catch ( cms::Exception& ex ) { edm::LogWarning("ElectronAndPhotonSuperClusterProducer") << "Can't get photon Collection: " << m_PhotonLabel; }
  }

  edm::LogInfo("ElectronAndPhotonSuperClusterProducer") << "ElectronAndPhotonSuperClusterProducer starting, original collection size = " 
                                                        << electrons->size();

  reco::SuperClusterCollection::const_iterator it = superclusters->begin();
  for(unsigned int isclu = 0; isclu < superclusters->size(); isclu++) {
    reco::SuperClusterRef scRef(superclusters,isclu);
    bool linked=false;

    // look if it is linked to an electron
    for(unsigned int iele = 0; iele < electrons->size(); iele++) {
      const reco::GsfElectronRef electronRef(electrons,iele);
      if ( !(electronRef.isNull()) && electronRef->ecalDrivenSeed() ) {
        reco::SuperClusterRef eleSc = electronRef->superCluster();
        if ( eleSc.isNonnull() ) { 
          // match by energy and position, since the input is a merged EB+EE collection where the ref is lost
          if( eleSc->energy()==scRef->energy() && eleSc->position()==scRef->position() ) {
            selected->push_back(*it);
            linked=true;
            break;
          }
        }
      }
    }

    if(m_includePhotonSuperClusters) {
      // look if it is linked to a photon
      for(unsigned int ipho = 0; ipho < photons->size(); ipho++) {
        const reco::PhotonRef photonRef(photons,ipho);
        if ( !(photonRef.isNull()) ) {
          reco::SuperClusterRef phoSc = photonRef->superCluster();
          if ( phoSc.isNonnull() ) { 
            // match by energy and position, since the input is a merged EB+EE collection where the ref is lost
            if( phoSc->energy()==scRef->energy() && phoSc->position()==scRef->position() && !linked) {
              selected->push_back(*it);
              break;
            }
          }
        }
      }
    }

    it++;
  }
  
  edm::LogInfo("ElectronAndPhotonSuperClusterProducer") << "ElectronAndPhotonSuperClusterProducer ending, final collection size = " 
                                               << selected->size();

  iEvent.put( selected );
  
}      

// ------------ method called once each job just before starting event loop  ------------
void 
ElectronAndPhotonSuperClusterProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElectronAndPhotonSuperClusterProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
ElectronAndPhotonSuperClusterProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ElectronAndPhotonSuperClusterProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ElectronAndPhotonSuperClusterProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ElectronAndPhotonSuperClusterProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}
