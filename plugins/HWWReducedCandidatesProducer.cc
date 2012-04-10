// -*- C++ -*-
//
// Package:    HWWReducedCandidatesProducer
// Class:      HWWReducedCandidatesProducer
// 
/**\class HWWReducedCandidatesProducer HWWReducedCandidatesProducer.cc WWAnalysis/HWWReducedCandidatesProducer/src/HWWReducedCandidatesProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Mon Jun 13 00:56:48 CEST 2011
// $Id: HWWReducedCandidatesProducer.cc,v 1.1 2011/06/13 13:42:31 mangano Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/GsfTrackReco/interface/GsfTrack.h>

//
// class declaration
//

class HWWReducedCandidatesProducer : public edm::EDProducer {
   public:
      explicit HWWReducedCandidatesProducer(const edm::ParameterSet&);
      ~HWWReducedCandidatesProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
      edm::InputTag srcCands_;
      edm::InputTag srcVertices_;
      double dz_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
HWWReducedCandidatesProducer::HWWReducedCandidatesProducer(const edm::ParameterSet& iConfig):
  srcCands_( iConfig.getParameter<edm::InputTag>("srcCands") ),
  srcVertices_( iConfig.getParameter<edm::InputTag>("srcVertices") ),
  dz_( iConfig.getParameter<double>("dz") )
{
  using namespace reco;
  produces<CandidateCollection>(); 
}


HWWReducedCandidatesProducer::~HWWReducedCandidatesProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
HWWReducedCandidatesProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco; 
  std::auto_ptr<CandidateCollection> comp( new CandidateCollection );

  Handle<PFCandidateCollection>   cands;
  iEvent.getByLabel( srcCands_, cands );

  edm::Handle<vector<reco::Vertex> >  vtxH;
  iEvent.getByLabel( srcVertices_ ,vtxH);

  reco::Vertex vtx = vtxH->front();

  for( PFCandidateCollection::const_iterator c = cands->begin(); c != cands->end(); ++c ) {
     if(c->charge()==0) continue;
     double dz(99999.);
     if(c->trackRef().isNonnull()) 
       dz = fabs(c->trackRef()->dz(vtx.position()));
     else if(c->gsfTrackRef().isNonnull())
       dz = fabs(c->gsfTrackRef()->dz(vtx.position()));
     else if(c->muonRef().isNonnull() && c->muonRef()->innerTrack().isNonnull())
       dz = fabs(c->muonRef()->innerTrack()->dz(vtx.position()));
     else {
       cout << "WARNING: found charged PF candidate without any track ref" << endl;
       continue;
     }
     if(dz>dz_) continue;
     std::auto_ptr<Candidate> cand( new LeafCandidate( * c ) );
     comp->push_back( cand.release() );
   }

   iEvent.put( comp );
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
HWWReducedCandidatesProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HWWReducedCandidatesProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
HWWReducedCandidatesProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
HWWReducedCandidatesProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HWWReducedCandidatesProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HWWReducedCandidatesProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HWWReducedCandidatesProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HWWReducedCandidatesProducer);
