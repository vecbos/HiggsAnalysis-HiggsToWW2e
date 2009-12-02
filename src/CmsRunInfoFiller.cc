#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsRunInfoFiller.h"

using namespace edm;
using namespace std;

CmsRunInfoFiller::CmsRunInfoFiller(CmsTree *cmsTree): 
  privateData_(new CmsRunInfoFillerData) {

  cmstree=cmsTree;
  privateData_->initialise();

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

}

void CmsRunInfoFiller::writeRunInfoToTree(const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                          bool dumpData) {

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

  const std::vector<bool>& gtTechTrigWord = gtReadoutRecord->technicalTriggerWord();

  *(privateData_->nl1Technical) = gtTechTrigWord.size();
  int techBlockSize = gtTechTrigWord.size();
  cmstree->column("nl1Technical",techBlockSize,0,"L1T");
  
  for (unsigned int iBit=0;iBit<gtTechTrigWord.size();++iBit) {
    privateData_->l1Technical->push_back(gtTechTrigWord[iBit]);
  }

  cmstree->column("l1Technical", *privateData_->l1Technical, "nl1Technical", 0, "L1T");
  
  const std::vector<bool>& gtDecisionTrigWord = gtReadoutRecord->decisionWord();

  *(privateData_->nl1Global) = gtDecisionTrigWord.size();
  int gtBlockSize = gtDecisionTrigWord.size();
  cmstree->column("nl1Global",gtBlockSize,0,"L1T");

  for (unsigned int iBit=0;iBit<gtDecisionTrigWord.size();++iBit) {
    privateData_->l1Global->push_back(gtDecisionTrigWord[iBit]);
  }

  cmstree->column("l1Global", *privateData_->l1Global, "nl1Global", 0, "L1T");

  treeRunInfo();

  if(dumpData) cmstree->dumpData();

}

void CmsRunInfoFiller::treeRunInfo() {

  cmstree->column("runNumber", *privateData_->run, 0, "L1T");
  cmstree->column("eventNumber", *privateData_->event, 0, "L1T");
  cmstree->column("lumiBlock", *privateData_->lumisection, 0, "L1T");
  cmstree->column("bunchCrossing", *privateData_->bx, 0, "L1T");
  cmstree->column("orbitNumber", *privateData_->orbit, 0, "L1T");
  
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

}

