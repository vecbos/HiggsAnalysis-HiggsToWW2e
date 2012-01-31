#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "RecoCaloTools/Navigation/interface/EcalPreshowerNavigator.h"
#include "TH1F.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/ESHelper.h"

using namespace std;

ESHelper::ESHelper(std::map<DetId, EcalRecHit> rechits_map, const CaloSubdetectorGeometry* geometry_p, const CaloSubdetectorTopology* topology_p) {
  _rechits_map = rechits_map;
  _geometry = geometry_p;
  _topology = topology_p;
}

vector<float> ESHelper::getESHits(double X, double Y, double Z, int row) {
  //cout<<row<<endl;

  vector<float> esHits;

  //double X = bcPtr->x();
  //double Y = bcPtr->y();
  //double Z = bcPtr->z();
  const GlobalPoint point(X,Y,Z);

  DetId esId1 = (dynamic_cast<const EcalPreshowerGeometry*>(_geometry))->getClosestCellInPlane(point, 1);
  DetId esId2 = (dynamic_cast<const EcalPreshowerGeometry*>(_geometry))->getClosestCellInPlane(point, 2);
  ESDetId esDetId1 = (esId1 == DetId(0)) ? ESDetId(0) : ESDetId(esId1);
  ESDetId esDetId2 = (esId2 == DetId(0)) ? ESDetId(0) : ESDetId(esId2);  

  map<DetId, EcalRecHit>::iterator it;
  ESDetId next;
  ESDetId strip1;
  ESDetId strip2;

  strip1 = esDetId1;
  strip2 = esDetId2;

  EcalPreshowerNavigator theESNav1(strip1, _topology);
  theESNav1.setHome(strip1);
  
  EcalPreshowerNavigator theESNav2(strip2, _topology);
  theESNav2.setHome(strip2);

  if (row == 1) {
    if (strip1 != ESDetId(0)) strip1 = theESNav1.north();
    if (strip2 != ESDetId(0)) strip2 = theESNav2.east();
  } else if (row == -1) {
    if (strip1 != ESDetId(0)) strip1 = theESNav1.south();
    if (strip2 != ESDetId(0)) strip2 = theESNav2.west();
  }
  // Plane 1 
  if (strip1 == ESDetId(0)) {
    for (int i=0; i<31; ++i) esHits.push_back(0);
  } else {
    
    it = _rechits_map.find(strip1);
    if (it->second.energy() > 1.0e-10) esHits.push_back(it->second.energy());  
    else esHits.push_back(0);
    //cout<<"center : "<<strip1<<" "<<it->second.energy()<<endl;      

    // east road 
    for (int i=0; i<15; ++i) {
      next = theESNav1.east();
      if (next != ESDetId(0)) {
        it = _rechits_map.find(next);
        if (it->second.energy() > 1.0e-10) esHits.push_back(it->second.energy());  
        else esHits.push_back(0);
        //cout<<"east "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"east "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }

    // west road 
    theESNav1.setHome(strip1);
    theESNav1.home();
    for (int i=0; i<15; ++i) {
      next = theESNav1.west();
      if (next != ESDetId(0)) {
        it = _rechits_map.find(next);
        if (it->second.energy() > 1.0e-10) esHits.push_back(it->second.energy());  
        else esHits.push_back(0);
        //cout<<"west "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"west "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }
  }


  if (strip2 == ESDetId(0)) {
    for (int i=0; i<31; ++i) esHits.push_back(0);
  } else {

    it = _rechits_map.find(strip2);
    if (it->second.energy() > 1.0e-10) esHits.push_back(it->second.energy());  
    else esHits.push_back(0);
    //cout<<"center : "<<strip2<<" "<<it->second.energy()<<endl;      

    // north road 
    for (int i=0; i<15; ++i) {
      next = theESNav2.north();
      if (next != ESDetId(0)) {
        it = _rechits_map.find(next);
        if (it->second.energy() > 1.0e-10) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"north "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;  
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"north "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }

    // south road 
    theESNav2.setHome(strip2);
    theESNav2.home();
    for (int i=0; i<15; ++i) {
      next = theESNav2.south();
      if (next != ESDetId(0)) {
        it = _rechits_map.find(next);
        if (it->second.energy() > 1.0e-10) esHits.push_back(it->second.energy());  
        else esHits.push_back(0);
        //cout<<"south "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"south "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }
  }

  return esHits;
}

vector<float> ESHelper::getESShape(vector<float> ESHits0) {
  const int nBIN = 21;
  vector<float> esShape;

  TH1F *htmpF = new TH1F("htmpF","",nBIN,0,nBIN);
  TH1F *htmpR = new TH1F("htmpR","",nBIN,0,nBIN);
  htmpF->Reset(); htmpR->Reset();

  Float_t effsigmaRR=0.;

  for(int ibin=0; ibin<((nBIN+1)/2); ibin++) {
    if (ibin==0) {
      htmpF->SetBinContent((nBIN+1)/2,ESHits0[ibin]);
      htmpR->SetBinContent((nBIN+1)/2,ESHits0[ibin+31]);
    } else { // hits sourd the seed
      htmpF->SetBinContent((nBIN+1)/2+ibin,ESHits0[ibin]);
      htmpF->SetBinContent((nBIN+1)/2-ibin,ESHits0[ibin+15]);
      htmpR->SetBinContent((nBIN+1)/2+ibin,ESHits0[ibin+31]);
      htmpR->SetBinContent((nBIN+1)/2-ibin,ESHits0[ibin+31+15]);
    }
  }

  // ---- Effective Energy Deposit Width ---- //
  double EffWidthSigmaXX = 0.;
  double EffWidthSigmaYY = 0.;
  double totalEnergyXX   = 0.;
  double totalEnergyYY   = 0.;
  double EffStatsXX      = 0.;
  double EffStatsYY      = 0.;
  for (int id_X=1; id_X<=21; id_X++) {
    totalEnergyXX  += htmpF->GetBinContent(id_X);
    EffStatsXX     += htmpF->GetBinContent(id_X)*(id_X-11)*(id_X-11);
    totalEnergyYY  += htmpR->GetBinContent(id_X);
    EffStatsYY     += htmpR->GetBinContent(id_X)*(id_X-11)*(id_X-11);
  }
  // If denominator == 0, effsigmaRR = 0;
  EffWidthSigmaXX  = (totalEnergyXX>0.)  ? sqrt(fabs(EffStatsXX  / totalEnergyXX))   : 0.;
  EffWidthSigmaYY  = (totalEnergyYY>0.)  ? sqrt(fabs(EffStatsYY  / totalEnergyYY))   : 0.;
  effsigmaRR =  ((totalEnergyXX  + totalEnergyYY) >0.) ? sqrt(EffWidthSigmaXX  * EffWidthSigmaXX  + EffWidthSigmaYY  * EffWidthSigmaYY)  : 0.;

  esShape.push_back(effsigmaRR);
  esShape.push_back(EffWidthSigmaXX);
  esShape.push_back(EffWidthSigmaYY);

  delete htmpF;
  delete htmpR;

  return esShape;
}

vector<float> ESHelper::getESEn(vector<float> ESHits0) {
  const int nBIN = 21;
  vector<float> esEn;

  float esf_e1 = 0.;    float esr_e1 = 0.;    // Sum of energy deposited in p/m  0 strip
  float esf_e3 = 0.;    float esr_e3 = 0.;    // Sum of energy deposited in p/m  1 strip
  float esf_e5 = 0.;    float esr_e5 = 0.;    // Sum of energy deposited in p/m  2 strips
  float esf_e7 = 0.;    float esr_e7 = 0.;    // Sum of energy deposited in p/m  3 strips
  float esf_e11= 0.;    float esr_e11= 0.;    // Sum of energy deposited in p/m  5 strips
  float esf_e21= 0.;    float esr_e21= 0.;    // Sum of energy deposited in p/m 10 strips

  // Insert energy from central strip
  esf_e1  += ESHits0[0];     esr_e1  += ESHits0[31];
  esf_e3  += ESHits0[0];     esr_e3  += ESHits0[31];
  esf_e5  += ESHits0[0];     esr_e5  += ESHits0[31];
  esf_e7  += ESHits0[0];     esr_e7  += ESHits0[31];
  esf_e11 += ESHits0[0];     esr_e11 += ESHits0[31];
  esf_e21 += ESHits0[0];     esr_e21 += ESHits0[31];
  // hits sround the central strip
  for(int ibin=1; ibin<((nBIN+1)/2); ibin++) {
    if (ibin<=1) {
      esf_e3 += ESHits0[ibin];
      esf_e3 += ESHits0[ibin+15];
      esr_e3 += ESHits0[ibin+31];
      esr_e3 += ESHits0[ibin+31+15];
    }
    if (ibin<=2) {
      esf_e5 += ESHits0[ibin];
      esf_e5 += ESHits0[ibin+15];
      esr_e5 += ESHits0[ibin+31];
      esr_e5 += ESHits0[ibin+31+15];
    }
    if (ibin<=3) {
      esf_e7 += ESHits0[ibin];
      esf_e7 += ESHits0[ibin+15];
      esr_e7 += ESHits0[ibin+31];
      esr_e7 += ESHits0[ibin+31+15];
    }
    if (ibin<=5) {
      esf_e11 += ESHits0[ibin];
      esf_e11 += ESHits0[ibin+15];
      esr_e11 += ESHits0[ibin+31];
      esr_e11 += ESHits0[ibin+31+15];
    }
    esf_e21 += ESHits0[ibin];
    esf_e21 += ESHits0[ibin+15];
    esr_e21 += ESHits0[ibin+31];
    esr_e21 += ESHits0[ibin+31+15];
  }

  esEn.push_back(esf_e1);  esEn.push_back(esr_e1);
  esEn.push_back(esf_e3);  esEn.push_back(esr_e3);
  esEn.push_back(esf_e5);  esEn.push_back(esr_e5);
  esEn.push_back(esf_e7);  esEn.push_back(esr_e7);
  esEn.push_back(esf_e11); esEn.push_back(esr_e11);
  esEn.push_back(esf_e21); esEn.push_back(esr_e21);

  return esEn;
}

