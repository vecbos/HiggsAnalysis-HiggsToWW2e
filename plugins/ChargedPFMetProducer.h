// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

class ChargedPFMetProducer : public edm::EDProducer {
    public:
        explicit ChargedPFMetProducer(const edm::ParameterSet&);
        ~ChargedPFMetProducer();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag collectionTag_;    
        edm::InputTag vertexTag_;
        double dzCut_;

};
