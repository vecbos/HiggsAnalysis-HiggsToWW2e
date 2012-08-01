#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPdfWeightFiller.h"

using namespace edm;
using namespace std;

CmsPdfWeightFiller::CmsPdfWeightFiller(CmsTree *cmsTree) {
  cmstree=cmsTree;
  trkIndexName_ = new std::string("n");
}

CmsPdfWeightFiller::~CmsPdfWeightFiller() {
}

void CmsPdfWeightFiller::writePdfWeightToTree(edm::InputTag pdfSet, const edm::Event& iEvent, const edm::EventSetup& iSetup,
					      const std::string &columnPrefix, const std::string &columnSuffix,
					      bool dumpData) {

  std::string label = pdfSet.label();
  if(label.size()>0) {
    edm::Handle<std::vector<double> > weightHandle;
    iEvent.getByLabel(pdfSet, weightHandle);
    
    std::vector<double> myweights = (*weightHandle);
    
    std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix;
    cmstree->column(nCandString.c_str(),int(myweights.size()),0,"Reco");
    cmstree->column((columnPrefix+"w"+columnSuffix).c_str(), myweights, nCandString.c_str(), 0, "Reco");
  }
  if(dumpData) cmstree->dumpData();

}



