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
#include "CMGTools/External/interface/PileupJetIdAlgo.h"

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
  PileupJetIdAlgo *puIdAlgo_;
};

using namespace edm;
using namespace reco;
using namespace std;

mvaJetIDMapProd::mvaJetIDMapProd(const edm::ParameterSet& iConfig) :
  
  vtxLabel_(iConfig.getUntrackedParameter<edm::InputTag>("vtxLabel")),
  uncorrJetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("uncorrJetLabel")),
  corrJetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("corrJetLabel")) {
  
  puIdAlgo_ = new PileupJetIdAlgo(iConfig.getParameter<edm::ParameterSet>("puJetIDAlgo"));
  
  // MVA output
  produces<edm::ValueMap<float> >("mva");
  
  // MVA inputs
  produces<edm::ValueMap<float> >("nCharged");
  produces<edm::ValueMap<float> >("nNeutrals");
  produces<edm::ValueMap<float> >("dZ");
  produces<edm::ValueMap<float> >("nParticles");
  produces<edm::ValueMap<float> >("dR2Mean");
  produces<edm::ValueMap<float> >("dRMean");
  produces<edm::ValueMap<float> >("frac01");
  produces<edm::ValueMap<float> >("frac02");
  produces<edm::ValueMap<float> >("frac03");
  produces<edm::ValueMap<float> >("frac04");
  produces<edm::ValueMap<float> >("frac05");
  produces<edm::ValueMap<float> >("beta");
  produces<edm::ValueMap<float> >("betastar");
  produces<edm::ValueMap<float> >("betastarclassic");
}

void mvaJetIDMapProd::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // uncorrected PF jets
  edm::Handle<edm::View<reco::Candidate> > uncorrPFJetCollectionHandle;
  try { iEvent.getByLabel(uncorrJetLabel_,uncorrPFJetCollectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("mvaJetIdMapProducer") << "Can't get uncorrected jets"; }
  const edm::View<reco::Candidate> *uncorrPFJetColl = uncorrPFJetCollectionHandle.product();

  // fully corrected PF jets  
  edm::Handle< edm::View<reco::Candidate> > corrPFJetCollectionHandle;
  try { iEvent.getByLabel(corrJetLabel_,corrPFJetCollectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("mvaJetIdMapProducer") << "Can't get corrected jets"; }
  const edm::View<reco::Candidate> *corrPFJetColl = corrPFJetCollectionHandle.product();

  // vertices 
  edm::Handle<reco::VertexCollection> primaryVertex;
  try { iEvent.getByLabel(vtxLabel_, primaryVertex); }
  catch(cms::Exception& ex ) {edm::LogWarning("mvaJetIdMapProducer") << "Can't get vertices"; }
  reco::VertexCollection vertexCollection = *(primaryVertex.product());
  reco::VertexCollection::const_iterator vMax = primaryVertex->begin();
  reco::Vertex pv;
  if (primaryVertex->size()>0) pv = *vMax;

  // ----------------------------------------
  // MVA output producer
  std::vector<float> mvaidV;
  std::auto_ptr<edm::ValueMap<float> > mvaidM(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler mvaidF(*mvaidM);

  // MVA inputs producer
  std::vector<float> mvaidV1, mvaidV2, mvaidV3,  mvaidV4,  mvaidV5,  mvaidV6,  mvaidV7;
  std::vector<float> mvaidV8, mvaidV9, mvaidV10, mvaidV11, mvaidV12, mvaidV13, mvaidV14;
  std::auto_ptr<edm::ValueMap<float> > mvaidM1(new edm::ValueMap<float> ());
  std::auto_ptr<edm::ValueMap<float> > mvaidM2(new edm::ValueMap<float> ());
  std::auto_ptr<edm::ValueMap<float> > mvaidM3(new edm::ValueMap<float> ());
  std::auto_ptr<edm::ValueMap<float> > mvaidM4(new edm::ValueMap<float> ());
  std::auto_ptr<edm::ValueMap<float> > mvaidM5(new edm::ValueMap<float> ());
  std::auto_ptr<edm::ValueMap<float> > mvaidM6(new edm::ValueMap<float> ());
  std::auto_ptr<edm::ValueMap<float> > mvaidM7(new edm::ValueMap<float> ());
  std::auto_ptr<edm::ValueMap<float> > mvaidM8(new edm::ValueMap<float> ());
  std::auto_ptr<edm::ValueMap<float> > mvaidM9(new edm::ValueMap<float> ());
  std::auto_ptr<edm::ValueMap<float> > mvaidM10(new edm::ValueMap<float> ());
  std::auto_ptr<edm::ValueMap<float> > mvaidM11(new edm::ValueMap<float> ());
  std::auto_ptr<edm::ValueMap<float> > mvaidM12(new edm::ValueMap<float> ());
  std::auto_ptr<edm::ValueMap<float> > mvaidM13(new edm::ValueMap<float> ());
  std::auto_ptr<edm::ValueMap<float> > mvaidM14(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler mvaidF1(*mvaidM1);
  edm::ValueMap<float>::Filler mvaidF2(*mvaidM2);
  edm::ValueMap<float>::Filler mvaidF3(*mvaidM3);
  edm::ValueMap<float>::Filler mvaidF4(*mvaidM4);
  edm::ValueMap<float>::Filler mvaidF5(*mvaidM5);
  edm::ValueMap<float>::Filler mvaidF6(*mvaidM6);
  edm::ValueMap<float>::Filler mvaidF7(*mvaidM7);
  edm::ValueMap<float>::Filler mvaidF8(*mvaidM8);
  edm::ValueMap<float>::Filler mvaidF9(*mvaidM9);
  edm::ValueMap<float>::Filler mvaidF10(*mvaidM10);
  edm::ValueMap<float>::Filler mvaidF11(*mvaidM11);
  edm::ValueMap<float>::Filler mvaidF12(*mvaidM12);
  edm::ValueMap<float>::Filler mvaidF13(*mvaidM13);
  edm::ValueMap<float>::Filler mvaidF14(*mvaidM14);

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
	  PileupJetIdentifier lPUJetId =  puIdAlgo_->computeIdVariables(pUncorrPFJet,1.,&pv,vertexCollection,true);
	  mvaidV.push_back(lPUJetId.mva());

	  mvaidV1.push_back(lPUJetId.nCharged());
	  mvaidV2.push_back(lPUJetId.nNeutrals());
	  mvaidV3.push_back(lPUJetId.dZ());
	  mvaidV4.push_back(lPUJetId.nParticles());
	  mvaidV5.push_back(lPUJetId.dR2Mean());
	  mvaidV6.push_back(lPUJetId.dRMean());
	  mvaidV7.push_back(lPUJetId.frac01());
	  mvaidV8.push_back(lPUJetId.frac02());
	  mvaidV9.push_back(lPUJetId.frac03());
	  mvaidV10.push_back(lPUJetId.frac04());
	  mvaidV11.push_back(lPUJetId.frac05());
	  mvaidV12.push_back(lPUJetId.beta());
	  mvaidV13.push_back(lPUJetId.betaStar());
	  mvaidV14.push_back(lPUJetId.betaStarClassic());

	  isMatched = true;
	  break ;
	}
      }
    }
    if (!isMatched) { mvaidV.push_back(-999.); }
    if (!isMatched) {
      mvaidV1.push_back(-999.); 
      mvaidV2.push_back(-999.); 
      mvaidV3.push_back(-999.); 
      mvaidV4.push_back(-999.); 
      mvaidV5.push_back(-999.); 
      mvaidV6.push_back(-999.); 
      mvaidV7.push_back(-999.); 
      mvaidV8.push_back(-999.); 
      mvaidV9.push_back(-999.); 
      mvaidV10.push_back(-999.); 
      mvaidV11.push_back(-999.); 
      mvaidV12.push_back(-999.); 
      mvaidV13.push_back(-999.); 
      mvaidV14.push_back(-999.); 
    }
  }
	
  // MVA output
  mvaidF.insert(corrPFJetCollectionHandle,mvaidV.begin(),mvaidV.end());
  mvaidF.fill();

  // MVA inputs
  mvaidF1.insert(corrPFJetCollectionHandle,mvaidV1.begin(),mvaidV1.end());
  mvaidF2.insert(corrPFJetCollectionHandle,mvaidV2.begin(),mvaidV2.end());
  mvaidF3.insert(corrPFJetCollectionHandle,mvaidV3.begin(),mvaidV3.end());
  mvaidF4.insert(corrPFJetCollectionHandle,mvaidV4.begin(),mvaidV4.end());
  mvaidF5.insert(corrPFJetCollectionHandle,mvaidV5.begin(),mvaidV5.end());
  mvaidF6.insert(corrPFJetCollectionHandle,mvaidV6.begin(),mvaidV6.end());
  mvaidF7.insert(corrPFJetCollectionHandle,mvaidV7.begin(),mvaidV7.end());
  mvaidF8.insert(corrPFJetCollectionHandle,mvaidV8.begin(),mvaidV8.end());
  mvaidF9.insert(corrPFJetCollectionHandle,mvaidV9.begin(),mvaidV9.end());
  mvaidF10.insert(corrPFJetCollectionHandle,mvaidV10.begin(),mvaidV10.end());
  mvaidF11.insert(corrPFJetCollectionHandle,mvaidV11.begin(),mvaidV11.end());
  mvaidF12.insert(corrPFJetCollectionHandle,mvaidV12.begin(),mvaidV12.end());
  mvaidF13.insert(corrPFJetCollectionHandle,mvaidV13.begin(),mvaidV13.end());
  mvaidF14.insert(corrPFJetCollectionHandle,mvaidV14.begin(),mvaidV14.end());
  //
  mvaidF1.fill();
  mvaidF2.fill();
  mvaidF3.fill();
  mvaidF4.fill();
  mvaidF5.fill();
  mvaidF6.fill();
  mvaidF7.fill();
  mvaidF8.fill();
  mvaidF9.fill();
  mvaidF10.fill();
  mvaidF11.fill();
  mvaidF12.fill();
  mvaidF13.fill();
  mvaidF14.fill();


  // putting everything into the event
  iEvent.put(mvaidM,"mva");
  iEvent.put(mvaidM1,"nCharged");
  iEvent.put(mvaidM2,"nNeutrals");
  iEvent.put(mvaidM3,"dZ");
  iEvent.put(mvaidM4,"nParticles");
  iEvent.put(mvaidM5,"dR2Mean");
  iEvent.put(mvaidM6,"dRMean");
  iEvent.put(mvaidM7,"frac01");
  iEvent.put(mvaidM8,"frac02");
  iEvent.put(mvaidM9,"frac03");
  iEvent.put(mvaidM10,"frac04");
  iEvent.put(mvaidM11,"frac05");
  iEvent.put(mvaidM12,"beta");
  iEvent.put(mvaidM13,"betastar");
  iEvent.put(mvaidM14,"betastarclassic");
}

mvaJetIDMapProd::~mvaJetIDMapProd() { 

  delete puIdAlgo_;
}
void mvaJetIDMapProd::beginJob() { }
void mvaJetIDMapProd::endJob() { }

DEFINE_FWK_MODULE(mvaJetIDMapProd);
