// my includes
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackExtraBase.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWEleId.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h" 
#include "Geometry/CaloGeometry/interface/CaloGeometry.h" 
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//CLHEP
#include <CLHEP/Vector/LorentzVector.h>

#include <iostream>


HWWEleId::HWWEleId (const edm::ParameterSet& conf) :
  m_barrelClusterShapeAss (conf.getParameter<edm::InputTag>
    ("barrelClusterShapeAssociation")),
  m_endcapClusterShapeAss (conf.getParameter<edm::InputTag>
    ("endcapClusterShapeAssociation"))
{
  
  hoeCut[0] = conf.getUntrackedParameter<double> ("hoeCut_0", 0.05) ;    
  e9o25Cut[0] = conf.getUntrackedParameter<double> ("e9o25Cut_0", 0.80) ;    
  detaCut[0] = conf.getUntrackedParameter<double> ("detaCut_0", 0.004) ;   
  dphiInCut[0] = conf.getUntrackedParameter<double> ("dphiInCut_0", 0.02) ;     
  dphiOutCut[0] = conf.getUntrackedParameter<double> ("dphiOutCut_0", 0.02) ; 
  seeInfCut[0] = conf.getUntrackedParameter<double> ("seeInfCut_0", 0.005) ;   
  seeSupCut[0] = conf.getUntrackedParameter<double> ("seeSupCut_0", 0.011) ;   
  eopOutInfCut[0] = conf.getUntrackedParameter<double> ("eopOutInfCut_0", 0.6) ;     
  eopOutSupCut[0] = conf.getUntrackedParameter<double> ("eopOutSupCut_0", 10000.) ;  
  
  hoeCut[1] = conf.getUntrackedParameter<double> ("hoeCut_1", 0.05) ;    
  e9o25Cut[1] = conf.getUntrackedParameter<double> ("e9o25Cut_1", 0.65) ;    
  detaCut[1] = conf.getUntrackedParameter<double> ("detaCut_1", 0.004) ;    
  dphiInCut[1] = conf.getUntrackedParameter<double> ("dphiInCut_1", 0.03) ;     
  dphiOutCut[1] = conf.getUntrackedParameter<double> ("dphiOutCut_1", 1000.) ; 
  seeInfCut[1] = conf.getUntrackedParameter<double> ("seeInfCut_1", 0.005) ;   
  seeSupCut[1] = conf.getUntrackedParameter<double> ("seeSupCut_1", 0.011) ;   
  eopOutInfCut[1] = conf.getUntrackedParameter<double> ("eopOutInfCut_1", 0.75) ;     
  eopOutSupCut[1] = conf.getUntrackedParameter<double> ("NAME_1", 10000.) ; 
                                                       
  hoeCut[2] = conf.getUntrackedParameter<double> ("hoeCut_2", 0.05) ;    
  e9o25Cut[2] = conf.getUntrackedParameter<double> ("e9o25Cut_2", 0.75) ;    
  detaCut[2] = conf.getUntrackedParameter<double> ("detaCut_2", 0.004) ;    
  dphiInCut[2] = conf.getUntrackedParameter<double> ("dphiInCut_2", 0.02) ;     
  dphiOutCut[2] = conf.getUntrackedParameter<double> ("dphiOutCut_2", 0.02) ; 
  seeInfCut[2] = conf.getUntrackedParameter<double> ("seeInfCut_2", 0.005) ;   
  seeSupCut[2] = conf.getUntrackedParameter<double> ("seeSupCut_2", 0.011) ;   
  eopOutInfCut[2] = conf.getUntrackedParameter<double> ("eopOutInfCut_2", 0.75) ;     
  eopOutSupCut[2] = conf.getUntrackedParameter<double> ("eopOutSupCut_2", 10000.) ; 
                                                       
  hoeCut[3] = conf.getUntrackedParameter<double> ("hoeCut_3", 0.05) ;    
  e9o25Cut[3] = conf.getUntrackedParameter<double> ("e9o25Cut_3", 0.65) ;    
  detaCut[3] = conf.getUntrackedParameter<double> ("detaCut_3", 0.005) ;    
  dphiInCut[3] = conf.getUntrackedParameter<double> ("dphiInCut_3", 0.04) ;     
  dphiOutCut[3] = conf.getUntrackedParameter<double> ("dphiOutCut_3", 1000.) ; 
  seeInfCut[3] = conf.getUntrackedParameter<double> ("seeInfCut_3", 0.005) ;   
  seeSupCut[3] = conf.getUntrackedParameter<double> ("seeSupCut_3", 0.011) ;   
  eopOutInfCut[3] = conf.getUntrackedParameter<double> ("eopOutInfCut_3", 0.75) ;     
  eopOutSupCut[3] = conf.getUntrackedParameter<double> ("eopOutSupCut_3", 10000.) ; 
                                                       
  hoeCut[4] = conf.getUntrackedParameter<double> ("hoeCut_4", 0.07) ;    
  e9o25Cut[4] = conf.getUntrackedParameter<double> ("e9o25Cut_4", 0.80) ;    
  detaCut[4] = conf.getUntrackedParameter<double> ("detaCut_4", 0.005) ;    
  dphiInCut[4] = conf.getUntrackedParameter<double> ("dphiInCut_4", 0.04) ;     
  dphiOutCut[4] = conf.getUntrackedParameter<double> ("dphiOutCut_4", 0.02) ; 
  seeInfCut[4] = conf.getUntrackedParameter<double> ("seeInfCut_4", 0.008) ;   
  seeSupCut[4] = conf.getUntrackedParameter<double> ("seeSupCut_4", 0.03) ;    
  eopOutInfCut[4] = conf.getUntrackedParameter<double> ("eopOutInfCut_4", 0.5) ;      
  eopOutSupCut[4] = conf.getUntrackedParameter<double> ("eopOutSupCut_4", 10000.) ; 
                                                       
  hoeCut[5] = conf.getUntrackedParameter<double> ("hoeCut_5", 0.07) ;    
  e9o25Cut[5] = conf.getUntrackedParameter<double> ("e9o25Cut_5", 0.70) ;    
  detaCut[5] = conf.getUntrackedParameter<double> ("detaCut_5", 0.005) ;    
  dphiInCut[5] = conf.getUntrackedParameter<double> ("dphiInCut_5", 0.04) ;     
  dphiOutCut[5] = conf.getUntrackedParameter<double> ("dphiOutCut_5", 1000.) ; 
  seeInfCut[5] = conf.getUntrackedParameter<double> ("seeInfCut_5", 0.008) ;   
  seeSupCut[5] = conf.getUntrackedParameter<double> ("seeSupCut_5", 0.03) ;    
  eopOutInfCut[5] = conf.getUntrackedParameter<double> ("eopOutInfCut_5", 0.8) ;      
  eopOutSupCut[5] = conf.getUntrackedParameter<double> ("eopOutSupCut_5", 10000.) ; 
                                                       
  hoeCut[6] = conf.getUntrackedParameter<double> ("hoeCut_6", 0.07) ;    
  e9o25Cut[6] = conf.getUntrackedParameter<double> ("e9o25Cut_6", 0.70) ;    
  detaCut[6] = conf.getUntrackedParameter<double> ("detaCut_6", 0.005) ;    
  dphiInCut[6] = conf.getUntrackedParameter<double> ("dphiInCut_6", 0.04) ;     
  dphiOutCut[6] = conf.getUntrackedParameter<double> ("dphiOutCut_6", 0.02) ; 
  seeInfCut[6] = conf.getUntrackedParameter<double> ("seeInfCut_6", 0.008) ;   
  seeSupCut[6] = conf.getUntrackedParameter<double> ("seeSupCut_6", 0.03) ;    
  eopOutInfCut[6] = conf.getUntrackedParameter<double> ("eopOutInfCut_6", 0.5) ;      
  eopOutSupCut[6] = conf.getUntrackedParameter<double> ("eopOutSupCut_6", 10000.) ; 
                                                       
  hoeCut[7] = conf.getUntrackedParameter<double> ("hoeCut_7", 0.07) ;    
  e9o25Cut[7] = conf.getUntrackedParameter<double> ("e9o25Cut_7", 0.65) ;    
  detaCut[7] = conf.getUntrackedParameter<double> ("detaCut_7", 0.005) ;    
  dphiInCut[7] = conf.getUntrackedParameter<double> ("dphiInCut_7", 0.05) ;     
  dphiOutCut[7] = conf.getUntrackedParameter<double> ("dphiOutCut_7", 1000.) ; 
  seeInfCut[7] = conf.getUntrackedParameter<double> ("seeInfCut_7", 0.008) ;   
  seeSupCut[7] = conf.getUntrackedParameter<double> ("seeSupCut_7", 0.022) ;   
  eopOutInfCut[7] = conf.getUntrackedParameter<double> ("eopOutInfCut_7", 0.8) ;      
  eopOutSupCut[7] = conf.getUntrackedParameter<double> ("eopOutSupCut_7", 10000.) ;  
  
//PG FIXME si puo' leggere come vettore?
//PG FIXME si puo' creare delle sotto-sezioni del cfg module?
}


// ------------------------------------------------------------


void 
HWWEleId::select (edm::Handle<collection> inputHandle, 
                  const edm::Event& iEvent, const edm::EventSetup& evtStp ) 
{
  m_selected.clear () ;

  typedef reco::BasicClusterShapeAssociationCollection clusterCollection ;
  // get barrel cluster shape collections
  edm::Handle<clusterCollection> barrelClShapes ;
  try { iEvent.getByLabel (m_barrelClusterShapeAss, barrelClShapes) ; }
  catch ( cms::Exception& ex ) { 
      edm::LogWarning ("HWWEleId") << "Can't get collection: " 
                                   << m_barrelClusterShapeAss ;
//      std::cerr << "Can't get collection: " 
//                << m_barrelClusterShapeAss << std::endl ;
      return ;                                  
    }

  // get endcap cluster shape collections
  edm::Handle<clusterCollection> endcapClShapes ;
  try { iEvent.getByLabel (m_endcapClusterShapeAss, endcapClShapes) ; }
  catch ( cms::Exception& ex ) { 
      edm::LogWarning ("HWWEleId") << "Can't get collection: " 
                                   << m_endcapClusterShapeAss ;
//      std::cerr << "Can't get collection: " 
//                << m_endcapClusterShapeAss << std::endl ;
      return ;                                  
    }

  //PG loop over input collection
  for (collection::const_iterator eleIt = inputHandle->begin () ;
       eleIt != inputHandle->end () ;
       ++eleIt)
    {
      // the track
      reco::GsfTrackRef myTrack = eleIt->gsfTrack () ;
      // the supercluster
      reco::SuperClusterRef mySc = eleIt->superCluster () ;
      // electron class
      int eleClass = eleIt->classification () ;
  
      // E and P 
      float eleFullCorrEne = eleIt->energy () ;         // corrected for Nxtal + f(eta) + E-P combination
      //      float caloCorrEne    = eleIt->caloEnergy () ;     // corrected for Nxtal + f(eta)
      //      float nxtalCorrEne   = mySc->energy () ;            // corrected for Nxtal 
      //      float rawEne         = mySc->rawEnergy () ;         // not corrected at all
      math::XYZVector trackAtVtx  = eleIt->momentum () ;    
      //      float trP = trackAtVtx.r () ;                       // not corrected
  
      // seed cluster variables 
      float seedEne  = mySc->seed ()->energy () ;
      float corrSeedEne = 10000.;            
//      GsfElectron MyEle2 = *eleIt ; //PG come mai questo??
      if ((int (eleClass/10) == 3) || (int (eleClass/10) == 13))
        {corrSeedEne = hwwEleSeedEnergyCorrector (*eleIt) ; }
//        {corrSeedEne = hwwEleSeedEnergyCorrector (MyEle2) ; }
      else {corrSeedEne = eleFullCorrEne ; }
      
      // shape variables
      float e3x3     = 1000. ;
      float e5x5     = 1000. ;
      float sigmaEE  = 1000. ;  
      DetId seedId   = mySc->seed ()->getHitsByDetId ()[0] ;
      if (seedId.subdetId () == EcalBarrel) {
        const reco::ClusterShapeRef & seedShapeRef = 
           barrelClShapes->find (mySc->seed ())->val ;
        e3x3    = seedShapeRef->e3x3 () ;
        e5x5    = seedShapeRef->e5x5 () ;
        sigmaEE = seedShapeRef->covEtaEta () ;
      }
      else{ 
        const reco::ClusterShapeRef& seedShapeRef = 
           endcapClShapes->find (mySc->seed ())->val ;
        e3x3    = seedShapeRef->e3x3 () ;
        e5x5    = seedShapeRef->e5x5 () ;
        sigmaEE = seedShapeRef->covEtaEta () ;
      }

      // problems with classification
      if (eleClass > 150) { continue ; }              
      // remove electrons within cracks for the moment!!!!!!
      if (eleClass == 40 ) { continue ; }  
      int index = -1 ;
      if ( eleClass == 0 ){ index = 0 ; }// barrel, golden
      else if ( eleClass == 10 ) { index = 1 ; }   // barrel, big brem
      else if ( eleClass == 20 ) { index = 2 ; }   // barrel, narrow
      else if ( eleClass > 29  && eleClass < 40) { index = 3 ; } // barrel, showering
      else if ( eleClass == 100) { index = 4 ; }  // endcap, golden
      else if ( eleClass == 110) { index = 5 ; }  // endcap, big brem 
      else if ( eleClass == 120) { index = 6 ; } // endcap, narrow
      else if ( eleClass > 129 && eleClass < 140) { index = 7 ; } // endcap, showering
      //FIXME else {/*if none is selected, is it possible?*/}     

      // miscellaneous for eleID 
      float hOverE         = eleIt->hadronicOverEm () ;          
      //      float EoPcorr        = eleIt->eSuperClusterOverP () ;      
      //      float EoPnotcorr     = rawEne/trP ;
      //      float EoPout         = eleIt->eSeedClusterOverPout () ;
      float EoPoutCorr     = eleIt->eSeedClusterOverPout ()*corrSeedEne/seedEne ;
      float deltaEtaAtVtx  = eleIt->deltaEtaSuperClusterTrackAtVtx () ;
      float deltaPhiAtVtx  = eleIt->deltaPhiSuperClusterTrackAtVtx () ;
      float deltaPhiAtCalo = eleIt->deltaPhiSeedClusterTrackAtCalo () ;

      if (( fabs (hOverE) > hoeCut[index] )                       ||	    
          ( e3x3/e5x5 < e9o25Cut[index] )                         ||    
          ( fabs (deltaEtaAtVtx) > detaCut[index])                ||	 	
          ( fabs (deltaPhiAtCalo) > dphiOutCut[index])            ||	    
          ( fabs (deltaPhiAtVtx) > dphiInCut[index])              ||	   
          ( sigmaEE < 0 )                                         ||
          ( sigmaEE < seeInfCut[index] * seeInfCut[index] )       ||
          ( sigmaEE > seeSupCut[index] * seeSupCut[index] )       ||
          ( EoPoutCorr < eopOutInfCut[index] )                    ||
          ( EoPoutCorr > eopOutSupCut[index] )) 
        { continue ; }
    
      //PG fill the selected's collection with survived electrons
      m_selected.push_back (& (* eleIt)) ;
    } //PG loop over input collection   

  return ;
}        


// ------------------------------------------------------------


double 
HWWEleId::hwwEleSeedEnergyCorrector (const electron & Electron) 
{
  // Eseed should be corrected for 
  // Escale = Nxtal + the part of f(eta) not related to brem subclusters
  // Tricky: for the time being only Nxtal
 
  reco::SuperClusterRef supercluster = Electron.superCluster () ;
  
  double Ebrem = 0.;  
  reco::basicCluster_iterator bc;
  for (bc = supercluster->clustersBegin () ; 
       bc != supercluster->clustersEnd () ; 
       ++bc) {
    LogDebug ("HWWEleId") << (*bc)->energy () 
			  << " " << supercluster->seed ()->energy ();
    Ebrem = Ebrem + (*bc)->energy () ;
  }
  Ebrem = Ebrem - supercluster->seed ()->energy () ;  
  return supercluster->energy () - Ebrem ;
}


// ------------------------------------------------------------


HWWEleId::~HWWEleId()
{}
