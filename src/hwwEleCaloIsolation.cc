// my includes
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/hwwEleCaloIsolation.h"

//CLHEP
#include <CLHEP/Vector/LorentzVector.h>

// temporary version.
// we want to check a cone between the HCAL cluster and the electron track.
// but this seems not to be possible for the time being....
// or at least I do not understand how to do....

using namespace std;
hwwEleCaloIsolation::hwwEleCaloIsolation (){}

hwwEleCaloIsolation::hwwEleCaloIsolation (const PixelMatchGsfElectron *gsfEle, const HBHERecHitCollection hcalrhColl, edm::ESHandle<CaloGeometry> caloGeo) : 
  _myGsfEle(gsfEle),
  _hcalrh(hcalrhColl),   
  _caloGeo(caloGeo)
{
  _extRadius = 0.20;   
}

hwwEleCaloIsolation::~hwwEleCaloIsolation (){}


void hwwEleCaloIsolation::setExtRadius (float extRadius){_extRadius = extRadius; }

float hwwEleCaloIsolation::getEtHcal () const
{
  float dummyEt = 0 ;

  // electron
  Hep3Vector elePAtVtx(_myGsfEle->trackMomentumAtVtx().x(), _myGsfEle->trackMomentumAtVtx().y(), _myGsfEle->trackMomentumAtVtx().z() );
  float ele_et = (_myGsfEle->caloEnergy())*(sin(elePAtVtx.theta()));    
  double seedEta = _myGsfEle->superCluster()->seed()->position().eta();
  double seedPhi = _myGsfEle->superCluster()->seed()->position().phi();

  // hcal calo rechits
  for(HBHERecHitCollection::const_iterator this_rh = _hcalrh.begin(); this_rh != _hcalrh.end(); this_rh++ ){ 

    HcalDetId hcalDetID(this_rh->id());
    const CaloCellGeometry* hcalCell = 
      _caloGeo->getSubdetectorGeometry(hcalDetID)->getGeometry(hcalDetID);
    double hcalEta   = hcalCell->getPosition().eta();
    double hcalPhi   = hcalCell->getPosition().phi();
    double hcalTheta = hcalCell->getPosition().theta(); 

    double PI = 3.141592653589;
    float deltaEta = hcalEta - seedEta;
    float deltaPhi = fabs(hcalPhi - seedPhi);
    if (deltaPhi>PI){ deltaPhi = 2*PI-deltaPhi;}
    double dr = sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);

    double hcalEt  = (*this_rh).energy()*sin(hcalTheta);     
    if ( (fabs(dr) < _extRadius) && (hcalEt > 0.) ){ dummyEt += hcalEt; } 
    
  }//end loop over calo rechit
  
  // sum rechit et / ele et
  dummyEt = dummyEt/ele_et; 

  return dummyEt;
}


// minimum distance from rechits upper a given Et cut
float hwwEleCaloIsolation::minDeltaR (float minEt) const
{
  float minDR = 100000. ;

 // electron
  Hep3Vector elePAtVtx(_myGsfEle->trackMomentumAtVtx().x(), _myGsfEle->trackMomentumAtVtx().y(), _myGsfEle->trackMomentumAtVtx().z() );
  double seedEta = _myGsfEle->superCluster()->seed()->position().eta();
  double seedPhi = _myGsfEle->superCluster()->seed()->position().phi();

  // hcal calo rechits
  for(HBHERecHitCollection::const_iterator this_rh = _hcalrh.begin(); this_rh != _hcalrh.end(); this_rh++ ){ 
    HcalDetId hcalDetID(this_rh->id());
    const CaloCellGeometry* hcalCell = 
      _caloGeo->getSubdetectorGeometry(hcalDetID)->getGeometry(hcalDetID);
    double hcalEta   = hcalCell->getPosition().eta();
    double hcalPhi   = hcalCell->getPosition().phi();
    double hcalTheta = hcalCell->getPosition().theta(); 

    double hcalEt  = (*this_rh).energy()*sin(hcalTheta);     
    if (hcalEt<minEt){continue;}

    double PI = 3.141592653589;
    float deltaEta = hcalEta - seedEta;
    float deltaPhi = fabs(hcalPhi - seedPhi);
    if (deltaPhi>PI){ deltaPhi = 2*PI-deltaPhi;}
    double dr = sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);

    if ( fabs(dr) < minDR ){ minDR = dr; }
    
  }//end loop over calo rechit
  
  return minDR;  
}

bool hwwEleCaloIsolation::isIsolated (float etCut) const
{
  bool dummyIsolation = true ;
  
  if (hwwEleCaloIsolation::getEtHcal() > etCut )   // 0.05 
    dummyIsolation = false ;
  
  return dummyIsolation ;
}
