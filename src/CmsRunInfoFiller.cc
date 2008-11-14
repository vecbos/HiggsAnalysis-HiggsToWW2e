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

  treeRunInfo();

  if(dumpData) cmstree->dumpData();

}

void CmsRunInfoFiller::treeRunInfo() {

  cmstree->column("runNumber", *privateData_->run, 0, "Reco");
  cmstree->column("eventNumber", *privateData_->event, 0, "Reco");
  cmstree->column("lumiBlock", *privateData_->lumisection, 0, "Reco");
  cmstree->column("bunchCrossing", *privateData_->bx, 0, "Reco");
  cmstree->column("orbitNumber", *privateData_->orbit, 0, "Reco");

}

void CmsRunInfoFillerData::initialise() {

  run = new int;
  event = new int;
  lumisection = new int;
  bx = new int;
  orbit = new int;

}

