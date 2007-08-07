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

  // electron et
  float ele_et = _myGsfEle->et() ;

  // electron: sc position
  math::XYZPoint theCaloPosition = _myGsfEle->caloPosition () ;
  GlobalPoint pclu (theCaloPosition.x(), theCaloPosition.y(), theCaloPosition.z() );

  // is isolated?
  const CaloGeometry* caloGeomP = _caloGeo.product();
  for (int subdet=0; subdet<=7; subdet++) {
    const CaloSubdetectorGeometry* sdg=caloGeomP->getSubdetectorGeometry(DetId::Hcal,subdet);
    
    if (sdg!=0) {
      // get the list of detids within range (from geometry)
      CaloSubdetectorGeometry::DetIdSet dis=sdg->getCells(pclu, _extRadius);
      
      // loop over detids...
      HBHERecHitCollection::const_iterator j,je=_hcalrh.end();      
      
      for (CaloSubdetectorGeometry::DetIdSet::iterator i=dis.begin(); i!=dis.end(); i++) {
	if (i->subdetId()!=subdet) continue; 

	j=_hcalrh.find(*i);	
	
	if (j!=je){
	  double hcalHit_eta = caloGeomP->getPosition(j->detid()).eta();
	  double hcalHit_Et  = j->energy()*sin(2*atan(exp(-hcalHit_eta))); 
	  // dummyEt += j->energy();   // why ene and not et????
	  dummyEt += hcalHit_Et;
	}
      }
    }
  }

  // sum rechit et / ele et
  dummyEt = dummyEt/ele_et; 
  
  return dummyEt;
}

bool hwwEleCaloIsolation::isIsolated (float etCut) const
{
  bool dummyIsolation = true ;
  
  if (hwwEleCaloIsolation::getEtHcal() > etCut )   // 0.05 
    dummyIsolation = false ;
  
  return dummyIsolation ;
}
