// system include files
#include <memory>
#include <cmath>
#include <algorithm>

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/METReco/interface/CommonMETData.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "RecoMET/METAlgorithms/interface/PFSpecificAlgo.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "HiggsAnalysis/HiggsToWW2e/plugins/ChargedPFMetProducer.h"

ChargedPFMetProducer::ChargedPFMetProducer(const edm::ParameterSet& iConfig):
  collectionTag_(iConfig.getParameter<edm::InputTag>("collectionTag")) {

  produces<reco::PFMETCollection>();
}

ChargedPFMetProducer::~ChargedPFMetProducer() { }
void ChargedPFMetProducer::beginJob() { }
void ChargedPFMetProducer::endJob() { } 

void ChargedPFMetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<edm::View<reco::Candidate> > reducedPFCands;
  iEvent.getByLabel(collectionTag_,  reducedPFCands);

  reco::Candidate::LorentzVector totalP4;
  float sumet = 0.0;
  edm::View<reco::Candidate>::const_iterator it;

  for(it= reducedPFCands->begin(); it!=reducedPFCands->end(); ++it) {
    totalP4 += it->p4();
    sumet += it->pt();
  }

  reco::Candidate::LorentzVector invertedP4(-totalP4);

  CommonMETData output;
  output.mex = invertedP4.px();
  output.mey = invertedP4.py();
  output.mez = invertedP4.pz();
  output.met = invertedP4.pt();
  output.sumet = sumet;
  output.phi = atan2(invertedP4.py(),invertedP4.px());
  PFSpecificAlgo pf;
  std::auto_ptr<reco::PFMETCollection> pfmetcoll;
  pfmetcoll.reset (new reco::PFMETCollection);
  pfmetcoll->push_back( pf.addInfo(reducedPFCands, output) );

  // and finally put it in the event
  iEvent.put( pfmetcoll );

}
