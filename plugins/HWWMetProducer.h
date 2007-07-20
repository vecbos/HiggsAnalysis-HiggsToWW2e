// -*- C++ -*-
// 
// Original Author:  Emanuele Di Marco
//         Created:  Mon Apr  2 17:28:56 CEST 2007
// $Id$
//
//

#ifndef HWWMetProducer_h
#define HWWMetProducer_h


#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TFile.h>
#include <TH1.h>

class HWWMetProducer: public edm::EDProducer {
public:
  explicit HWWMetProducer(const edm::ParameterSet&);
  explicit HWWMetProducer();
  virtual ~HWWMetProducer();
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void BookHistos();

private:
  // ---------- member data ---------------------------
  edm::InputTag genMetsTag_;
  edm::InputTag recoMetsTag_;
  double min_Pt_;
  bool produceGenMET_;
  std::string genMetsOut_;
  std::string recoMetsOut_;
  std::string genMetCandsOut_;
  std::string recoMetCandsOut_;
  std::string outputFile_;
  // ---------- monitoring histograms -----------------
  TFile *dataFile_;
  TH1F* hRecoMET_phi, *hRecoMET_MET, *hRecoMET_sig;
  TH1F* hGenMET_phi, *hGenMET_MET, *hGenMET_sig;
};
#endif 
