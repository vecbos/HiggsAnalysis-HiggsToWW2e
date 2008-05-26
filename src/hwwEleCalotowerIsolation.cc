#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/hwwEleCalotowerIsolation.h"

hwwEleCalotowerIsolation::hwwEleCalotowerIsolation(const reco::PixelMatchGsfElectron *gsfEle,
						   const CaloTowerCollection *calotowers)
{

  m_gsfEle = gsfEle;
  m_calotowers = calotowers;
  m_EtHcal = 0.;
  m_EtEcal = 0.;
  

}

float hwwEleCalotowerIsolation::getEtHcal(bool relative) {

  getEtTowers();

  float xval=0;

  if( relative ) {
    math::XYZVector electronMomentum = m_gsfEle->momentum () ;
    double elePt = sqrt( electronMomentum.perp2() );
    xval = m_EtHcal / elePt;
  }
  else xval = m_EtHcal;

  return xval;

}

float hwwEleCalotowerIsolation::getEtEcal(bool relative) {

  getEtTowers();

  float xval=0;

  if( relative ) {
    math::XYZVector electronMomentum = m_gsfEle->momentum () ;
    double elePt = sqrt( electronMomentum.perp2() );
    xval = m_EtEcal / elePt;
  }
  else xval = m_EtEcal;

  return xval;

}

void hwwEleCalotowerIsolation::getEtTowers() {

  double tmpSumHadEt = 0.0;
  double tmpSumEmEt = 0.0;

  reco::SuperClusterRef sc = m_gsfEle->get<reco::SuperClusterRef>();
  
  double scEta=sc->eta();
  double scPhi=sc->phi();
      
  // loop on calotowers
  CaloTowerCollection::const_iterator tower;
  for(tower = m_calotowers->begin(); tower != m_calotowers->end(); ++tower){

    double hadEt = tower->hadEt();
    double emEt = tower->emEt();
    
    double towerEta=tower->eta();
    double towerPhi=tower->phi();
    double twoPi= 2*M_PI;
    if(towerPhi<0) towerPhi+=twoPi;
    if(scPhi<0) scPhi+=twoPi;
    double deltaPhi=fabs(towerPhi-scPhi);
    if(deltaPhi>twoPi) deltaPhi-=twoPi;
    if(deltaPhi>M_PI) deltaPhi=twoPi-deltaPhi;
    double deltaEta = towerEta - scEta;
    
    double dr = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
    
    if ( dr < m_extRadius &&
	 dr >= m_intRadius ) {
      
      tmpSumHadEt += hadEt;
      tmpSumEmEt += emEt;
      
    }

  }

  m_EtHcal = tmpSumHadEt;
  m_EtEcal = tmpSumEmEt;

}
