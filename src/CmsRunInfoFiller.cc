#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsRunInfoFiller.h"


using namespace edm;
using namespace std;

CmsRunInfoFiller::CmsRunInfoFiller(CmsTree *cmsTree, bool isMC): 
  privateData_(new CmsRunInfoFillerData) {

  cmstree=cmsTree;
  privateData_->setMC(isMC);
  privateData_->initialise();
  isMC_=isMC;

}

CmsRunInfoFiller::~CmsRunInfoFiller() {

  // delete here the vector ptr's
  delete privateData_->run;
  delete privateData_->event;
  delete privateData_->lumisection;
  delete privateData_->bx;
  delete privateData_->orbit;
  delete privateData_->nl1Global;
  delete privateData_->nl1Technical;
  delete privateData_->l1Global;
  delete privateData_->l1Technical;
  if(isMC_) {
    delete privateData_->nBX;
    delete privateData_->nInteractions;
    delete privateData_->bxPU;
  }
  delete privateData_->beamSpotX;
  delete privateData_->beamSpotY;
  delete privateData_->beamSpotZ;
  delete privateData_;
}

void CmsRunInfoFiller::writeRunInfoToTree(const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                          bool dumpData) {

  privateData_->clearVectors();

  *(privateData_->run) = iEvent.id().run();
  *(privateData_->event) = iEvent.id().event();

  if( iEvent.isRealData() ) {
    *(privateData_->lumisection) = iEvent.luminosityBlock();
    *(privateData_->bx) = iEvent.bunchCrossing();
    *(privateData_->orbit) = iEvent.orbitNumber();
  } else {
    *(privateData_->lumisection) = -1;
    *(privateData_->bx) = -1;
    *(privateData_->orbit) = -1;
  }

  Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel( "offlineBeamSpot", beamSpotHandle);  
  if(beamSpotHandle.isValid()){
    *(privateData_->beamSpotX) = beamSpotHandle->position().X();
    *(privateData_->beamSpotY) = beamSpotHandle->position().Y();
    *(privateData_->beamSpotZ) = beamSpotHandle->position().Z();
  }else{
    *(privateData_->beamSpotX) = 0.;
    *(privateData_->beamSpotY) = 0.;
    *(privateData_->beamSpotZ) = 0.;
  }
  Handle<L1GlobalTriggerReadoutRecord> gtReadoutRecord;
  iEvent.getByLabel( "gtDigis", gtReadoutRecord );

  // Technical L1 trigger bits
  const std::vector<bool>& gtTechTrigWord = gtReadoutRecord->technicalTriggerWord();

  *(privateData_->nl1Technical) = gtTechTrigWord.size();
  int techBlockSize = gtTechTrigWord.size();

   // use words in 32 bits: use 30/32 bits per word
  int nTechWords = (techBlockSize-1)/30+1;
  privateData_->l1Technical->resize(nTechWords);
  for(int block=0; block<nTechWords; block++) (*privateData_->l1Technical)[block]=0;
  cmstree->column("nl1Technical",nTechWords,0,"L1T");
  
  for (unsigned int iBit=0;iBit<gtTechTrigWord.size();++iBit) {
    int block = iBit/30;
    int pos = iBit%30;
    int passed = ( gtTechTrigWord[iBit] ) ? 1 : 0;
    (*privateData_->l1Technical)[block] |= (passed << pos);
  }

  cmstree->column("l1Technical", *privateData_->l1Technical, "nl1Technical", 0, "L1T");
  
  // Physics L1 trigger bits
  const std::vector<bool>& gtDecisionTrigWord = gtReadoutRecord->decisionWord();

  *(privateData_->nl1Global) = gtDecisionTrigWord.size();
  int gtBlockSize = gtDecisionTrigWord.size();

  // use words in 32 bits: use 30/32 bits per word
  int nGTWords = (gtBlockSize-1)/30+1;
  privateData_->l1Global->resize(nGTWords);
  for(int block=0; block<nGTWords; block++) (*privateData_->l1Global)[block]=0;

  cmstree->column("nl1Global",nGTWords,0,"L1T");

  for (unsigned int iBit=0;iBit<gtDecisionTrigWord.size();++iBit) {
    int block = iBit/30;
    int pos = iBit%30;
    int passed = ( gtDecisionTrigWord[iBit] ) ? 1 : 0;
    (*privateData_->l1Global)[block] |= (passed << pos);
  }

  cmstree->column("l1Global", *privateData_->l1Global, "nl1Global", 0, "L1T");

  if(isMC_) {

    // Pile Up informations
    Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel("addPileupInfo", PupInfo);

    *(privateData_->nBX) = PupInfo->size();

    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI < PupInfo->end(); ++PVI) {
      privateData_->nInteractions->push_back(PVI->getPU_NumInteractions());
      privateData_->bxPU->push_back(PVI->getBunchCrossing());
    }

  }

  // get the rho parameter from Fastjet computation - isolation
  edm::Handle<double> rhoH, sigmaH;
  float rho;
  if( iEvent.getByLabel(edm::InputTag("kt6PFJetsForIsolation","rho"),rhoH) )
    rho = *rhoH;
  else
    rho = 0.;
  float sigma;
  if( iEvent.getByLabel(edm::InputTag("kt6PFJetsForIsolation","sigma"),sigmaH) )
    sigma = *sigmaH;
  else
    sigma = 0.;
  cmstree->column("rhoFastjet", rho, float(0.), "Iso");
  cmstree->column("sigmaFastjet", sigma, float(0.), "Iso");

  // get the rho parameter from Fastjet computation - jets
  edm::Handle<double> rhoHjets;
  float rhoJets;
  if( iEvent.getByLabel(edm::InputTag("kt6PFJets","rho"),rhoHjets) )
    rhoJets = *rhoHjets;
  else
    rhoJets = 0.;
  cmstree->column("rhoJetsFastJet", rhoJets, float(0.), "Jets");


  edm::Handle<double> rhoNoPuAllH;
  float rhoJets_nopu;
  if (iEvent.getByLabel(edm::InputTag("kt6PFJetsNoPU", "rho"), rhoNoPuAllH))
    rhoJets_nopu = *rhoNoPuAllH;
  else
    rhoJets_nopu = 0.;
  cmstree->column("rhoJetsFastJet_nopu", rhoJets_nopu, float(0.), "Jets");

//   edm::Handle<double> rhoNoPuH;
//   float rho_nopu;
//   if (iEvent.getByLabel(edm::InputTag("kt6PFJetsNoPUForIsolation", "rho"), rhoNoPuH))
//     rho_nopu = *rhoNoPuH;
//   else
//     rho_nopu = 0.;
//   cmstree->column("rhoFastJet_nopu", rho_nopu, float(0.), "Iso");




  // log errors (tracker failures)
  if(!isMC_) {
    edm::Handle< bool > tooManySeedsH;
    iEvent.getByLabel(edm::InputTag("tooManySeeds:TrackerLogError"), tooManySeedsH);
    bool tooManySeeds = *tooManySeedsH;
    
    edm::Handle< bool > tooManyClustersH;
    iEvent.getByLabel(edm::InputTag("tooManyClusters:TrackerLogError"), tooManyClustersH);
    bool tooManyClusters = *tooManyClustersH;
    
    int trackerFailures = ( tooManySeeds << 1 ) | tooManyClusters;
    
    cmstree->column("tooManyTrackerFailures", trackerFailures, 0, "LogError");
  }

  treeRunInfo();

  if(dumpData) cmstree->dumpData();

}

void CmsRunInfoFiller::treeRunInfo() {

  cmstree->column("runNumber", *privateData_->run, 0, "L1T");
  cmstree->column("eventNumber", *privateData_->event, uint64_t(0), "L1T");
  cmstree->column("lumiBlock", *privateData_->lumisection, 0, "L1T");
  cmstree->column("bunchCrossing", *privateData_->bx, 0, "L1T");
  cmstree->column("orbitNumber", *privateData_->orbit, 0, "L1T");

  cmstree->column("beamSpotX",*privateData_->beamSpotX, 0., "L1T");
  cmstree->column("beamSpotY",*privateData_->beamSpotY, 0., "L1T");
  cmstree->column("beamSpotZ",*privateData_->beamSpotZ, 0., "L1T");

  if(isMC_) {
    cmstree->column("nBX", *(privateData_->nBX), 0, "Sim");
    cmstree->column("nPU", *privateData_->nInteractions, "nBX", 0, "Sim");
    cmstree->column("bxPU", *privateData_->bxPU, "nBX", 0, "Sim");
  }

}

void CmsRunInfoFillerData::initialise() {

  run = new int;
  event = new uint64_t;
  lumisection = new int;
  bx = new int;
  orbit = new int;
  
  nl1Global = new int;
  nl1Technical = new int;
  l1Global = new vector<int>;
  l1Technical = new vector<int>;

  if(isMC_) {
    nBX = new int;
    nInteractions = new vector<int>;
    bxPU = new vector<int>;
  }

  beamSpotX = new double;
  beamSpotY = new double;
  beamSpotZ = new double;

}

void CmsRunInfoFillerData::clearVectors() {
  if(isMC_) {
    nInteractions->clear();
    bxPU->clear();
  }
}
