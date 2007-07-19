#ifndef HWWEleId_h
#define HWWEleId_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "FWCore/Framework/interface/Event.h"#include "FWCore/Framework/interface/EventSetup.h"#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <functional>
#include <vector>
#include <map>


class HWWEleId{

 public:

  typedef reco::PixelMatchGsfElectron electron ;
  typedef const electron * electronPtr ;
  typedef reco::PixelMatchGsfElectronCollection collection ;
  typedef std::vector<const reco::PixelMatchGsfElectron *> container ;
  typedef container::const_iterator const_iterator ;

  //! ctor
  HWWEleId (const edm::ParameterSet& conf) ;
  //!dtor
  ~HWWEleId () ;

  //! iterator to the begin of the selected collection
  const_iterator begin () const { return m_selected.begin () ; }
    //! iterator to the end of the selected collection
  const_iterator end () const { return m_selected.end () ; }
  //! the actual selector method 
  void select (edm::Handle<collection>, 
               const edm::Event&, const edm::EventSetup& ) ;
  
 private:
 
  //! comment
  double 
  hwwEleSeedEnergyCorrector (const electron & Electron) ;
 
 private:
  
  //! the selected collection
  container m_selected ;
  
  //! comment
  edm::InputTag m_barrelClusterShapeAss ;
  //! comment
  edm::InputTag m_endcapClusterShapeAss ;


  //! selection variables
  double hoeCut [8] ;  
  double e9o25Cut [8] ;  
  double detaCut [8] ;                               
  double dphiInCut [8] ;                        
  double dphiOutCut [8] ;  
  double seeInfCut [8] ;  
  double seeSupCut [8] ;  
  double eopOutInfCut [8] ;          
  double eopOutSupCut [8] ;  


 
};  

#endif

/* FIXME

nn mi e' chiara la funzione (o l'uso) del m_selected, perche'
c'e' il container di puntatori? come funziona il selector dell'object?
e' solo un trucco per non perdere tempo a copiare?
NON riesco a capire esattamente che cosa dovrebbe fare! 
se poi ritorna iteratori, a chi puntano questi?

*/
