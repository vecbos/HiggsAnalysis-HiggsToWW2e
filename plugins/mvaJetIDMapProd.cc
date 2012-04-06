#include <memory>
#include "Math/VectorUtil.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "CMG/JetIDAnalysis/interface/PileupJetIdNtupleAlgo.h"

#include "FWCore/Utilities/interface/Exception.h"


#include "DataFormats/Common/interface/ValueMap.h"
#include "TMath.h"

class mvaJetIDMapProd : public edm::EDProducer {
public:
  explicit mvaJetIDMapProd(const edm::ParameterSet&);
  ~mvaJetIDMapProd();
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  edm::InputTag vtxLabel_;
  edm::InputTag uncorrJetLabel_;
  edm::InputTag corrJetLabel_;
  PileupJetIdNtupleAlgo *puIdAlgo_;
};

using namespace edm;
using namespace reco;
using namespace std;

mvaJetIDMapProd::mvaJetIDMapProd(const edm::ParameterSet& iConfig) :
  
  vtxLabel_(iConfig.getUntrackedParameter<edm::InputTag>("vtxLabel")),
  uncorrJetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("uncorrJetLabel")),
  corrJetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("corrJetLabel")) {
  
  puIdAlgo_ = new PileupJetIdNtupleAlgo(iConfig.getParameter<edm::ParameterSet>("puJetIDAlgo"));
  
  produces<edm::ValueMap<float> >().setBranchAlias("mvaJetID");
}

void mvaJetIDMapProd::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // uncorrected PF jets
  edm::Handle<edm::View<reco::Candidate> > uncorrPFJetCollectionHandle;
  try { iEvent.getByLabel(uncorrJetLabel_,uncorrPFJetCollectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("NoPUMetProducer") << "Can't get uncorrected jets"; }
  const edm::View<reco::Candidate> *uncorrPFJetColl = uncorrPFJetCollectionHandle.product();

  // fully corrected PF jets  
  edm::Handle< edm::View<reco::Candidate> > corrPFJetCollectionHandle;
  try { iEvent.getByLabel(corrJetLabel_,corrPFJetCollectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("NoPUMetProducer") << "Can't get corrected jets"; }
  const edm::View<reco::Candidate> *corrPFJetColl = corrPFJetCollectionHandle.product();

  // vertices 
  edm::Handle<reco::VertexCollection> primaryVertex;
  try { iEvent.getByLabel(vtxLabel_, primaryVertex); }
  catch(cms::Exception& ex ) {edm::LogWarning("NoPUMetProducer") << "Can't get vertices"; }
  reco::VertexCollection vertexCollection = *(primaryVertex.product());
  reco::VertexCollection::const_iterator vMax = primaryVertex->begin();
  reco::Vertex pv;
  if (primaryVertex->size()>0) pv = *vMax;

  // ----------------------------------------
  std::vector<float> mvaidV;
  std::auto_ptr<edm::ValueMap<float> > mvaidM(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler mvaidF(*mvaidM);

  bool isMatched;
  for(int index2 = 0; index2 < (int)corrPFJetColl->size(); index2++) {   // corrected jets collection
    const Candidate *corrCand   = &(corrPFJetColl->at(index2));
    const PFJet     *pCorrPFJet = dynamic_cast< const PFJet * > ( &(*corrCand) );    

    isMatched = false;
    for(int index = 0; index < (int)uncorrPFJetColl->size(); index++) {      // uncorrected jets collection
      const Candidate *uncorrCand   = &(uncorrPFJetColl->at(index));
      const PFJet     *pUncorrPFJet = dynamic_cast< const PFJet * > ( &(*uncorrCand) );    
      
      if(  pUncorrPFJet->jetArea() == pCorrPFJet->jetArea() ) {      // to match corrected and uncorrected jets
	if(  fabs(pUncorrPFJet->eta() - pCorrPFJet->eta())<0.01 ) {  // to match corrected and uncorrected jets

	  double lJec = pCorrPFJet->pt()/pUncorrPFJet->pt();
	  PileupJetIdentifier lPUJetId =  puIdAlgo_->computeIdVariables(pUncorrPFJet,lJec,&pv,vertexCollection,true);
	  mvaidV.push_back(lPUJetId.mva());

	  isMatched = true;
	  break ;
	}
      }
    }
    if (!isMatched) { mvaidV.push_back(-999.); }
  }
	
  mvaidF.insert(corrPFJetCollectionHandle,mvaidV.begin(),mvaidV.end());

  mvaidF.fill();
  iEvent.put(mvaidM);
}

mvaJetIDMapProd::~mvaJetIDMapProd() { 

  delete puIdAlgo_;
}
void mvaJetIDMapProd::beginJob() { }
void mvaJetIDMapProd::endJob() { }

DEFINE_FWK_MODULE(mvaJetIDMapProd);
