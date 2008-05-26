#ifndef hwwEleCalotowerIsolation_h
#define hwwEleCalotowerIsolation_h

#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"

class hwwEleCalotowerIsolation {
 public:

  //! ctor
  hwwEleCalotowerIsolation() {}
  hwwEleCalotowerIsolation(const reco::PixelMatchGsfElectron *gsfEle,
			   const CaloTowerCollection *calotowers);

  //! dtor
  ~hwwEleCalotowerIsolation() {}

  //! methods
  void setIntRadius (float intRadius) { m_intRadius = intRadius; }
  void setExtRadius (float extRadius) { m_extRadius = extRadius; }

  //! get the Et sums in a cone. If relative => give sumEt/pT electron
  float getEtHcal(bool relative = true);
  float getEtEcal(bool relative = true);

 private:
  
  void getEtTowers();

  const reco::PixelMatchGsfElectron *m_gsfEle;
  const CaloTowerCollection *m_calotowers;

  float m_intRadius, m_extRadius;
  float m_EtHcal, m_EtEcal;

};

#endif

