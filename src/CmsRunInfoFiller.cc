#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsRunInfoFiller.h"

using namespace edm;
using namespace std;

CmsRunInfoFiller::CmsRunInfoFiller(CmsTree *cmsTree, bool isMC): 
  privateData_(new CmsRunInfoFillerData) {

  cmstree=cmsTree;
  privateData_->setMC(isMC_);
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
    delete privateData_->nInteractions;
    delete privateData_->PUzposition;
    delete privateData_->PUsumpTLowpT;
    delete privateData_->PUsumpTHighpT;
    delete privateData_->PUnTracksLowpT;
    delete privateData_->PUnTracksHighpT;
  }
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
    edm::Handle<PileupSummaryInfo> PupInfo;
    iEvent.getByLabel("addPileupInfo", PupInfo);
    
    *(privateData_->nInteractions) = PupInfo->getPU_NumInteractions();
    *(privateData_->PUzposition) = PupInfo->getPU_zpositions();
    *(privateData_->PUsumpTLowpT) = PupInfo->getPU_sumpT_lowpT();
    *(privateData_->PUsumpTHighpT) = PupInfo->getPU_sumpT_highpT();
    *(privateData_->PUnTracksLowpT) = PupInfo->getPU_ntrks_lowpT();
    *(privateData_->PUnTracksHighpT) = PupInfo->getPU_ntrks_highpT();
  }

  // get the rho parameter from Fastjet computation 
  edm::Handle<double> rhoH;
  iEvent.getByLabel(edm::InputTag("kt6PFJetsForIsolation","rho"),rhoH);
  float rho = *rhoH;
  cmstree->column("rhoFastjet", rho, float(0.), "Iso");
  
  treeRunInfo();

  if(dumpData) cmstree->dumpData();

}

void CmsRunInfoFiller::treeRunInfo() {

  cmstree->column("runNumber", *privateData_->run, 0, "L1T");
  cmstree->column("eventNumber", *privateData_->event, 0, "L1T");
  cmstree->column("lumiBlock", *privateData_->lumisection, 0, "L1T");
  cmstree->column("bunchCrossing", *privateData_->bx, 0, "L1T");
  cmstree->column("orbitNumber", *privateData_->orbit, 0, "L1T");

  if(isMC_) {
    cmstree->column("nPU", *privateData_->nInteractions, 0, "Sim");
    cmstree->column("zpositionPU", *privateData_->PUzposition, "nPU", 0, "Sim");
    cmstree->column("sumpTLowpTPU", *privateData_->PUsumpTLowpT, "nPU", 0, "Sim");
    cmstree->column("sumpTHighpTPU", *privateData_->PUsumpTHighpT, "nPU", 0, "Sim");
    cmstree->column("nTracksLowpTPU", *privateData_->PUnTracksLowpT, "nPU", 0, "Sim");
    cmstree->column("nTracksHighpTPU", *privateData_->PUnTracksHighpT, "nPU", 0, "Sim");
  }

}

void CmsRunInfoFillerData::initialise() {

  run = new int;
  event = new int;
  lumisection = new int;
  bx = new int;
  orbit = new int;
  
  nl1Global = new int;
  nl1Technical = new int;
  l1Global = new vector<int>;
  l1Technical = new vector<int>;

  if(isMC_) {
    nInteractions = new int;
    PUzposition = new vector<float>;
    PUsumpTLowpT = new vector<float>;
    PUsumpTHighpT = new vector<float>;
    PUnTracksLowpT = new vector<int>;
    PUnTracksHighpT = new vector<int>;
  }

}

void CmsRunInfoFillerData::clearVectors() {
  if(isMC_) {
    PUzposition->clear();
    PUsumpTLowpT->clear();
    PUsumpTHighpT->clear();
    PUnTracksLowpT->clear();
    PUnTracksHighpT->clear();
  }
}
