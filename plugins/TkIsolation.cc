//C++ includes
#include <vector>
#include <functional>

//ROOT includes
#include <Math/VectorUtil.h>

//CMSSW includes
#include "HiggsAnalysis/HiggsToWW2e/plugins/TkIsolation.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TrajectoryCleaning/interface/TrajectoryCleanerBySharedHits.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "TrackingTools/PatternTools/interface/TrajectoryFitter.h"
#include "RecoTracker/CkfPattern/interface/TransientInitialStateEstimator.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

using namespace std ;
using namespace ROOT::Math::VectorUtil ;

TkIsolation::TkIsolation ()
{
}


TkIsolation::TkIsolation (const reco::PixelMatchGsfElectron* electron,
				      const reco::TrackCollection* trackCollection)   :
  electron_(electron) ,
  trackCollection_(trackCollection)  
  
{
  electronCollection_ = 0 ; 
  extRadius_ = 0.25 ;
  intRadius_ = 0.015 ;
  ptLow_ = 1.5 ; 
  lip_ = 0.1 ; 
}

TkIsolation::TkIsolation (const reco::PixelMatchGsfElectron* electron, 
				      const reco::TrackCollection* trackCollection,
				      const reco::PixelMatchGsfElectronCollection* electronCollection) : 
  electron_(electron) ,
  trackCollection_(trackCollection) ,
  electronCollection_(electronCollection)  
{
  extRadius_ = 0.25 ;
  intRadius_ = 0.015 ;
  ptLow_ = 1.5 ; 
  lip_ = 0.1 ; 
}  

TkIsolation::~TkIsolation ()
{
}

void TkIsolation::setExtRadius (double extRadius)
{
  extRadius_ = extRadius ;
}

void TkIsolation::setIntRadius (double intRadius)
{  
  intRadius_ = intRadius ;
}

void TkIsolation::setPtLow (double ptLow)
{  
  ptLow_ = ptLow ;
}

void TkIsolation::setLip (double lip)
{  
  lip_ = lip ;
}

int TkIsolation::getNumberTracks () const
{  
  //counter for the tracks in the isolation cone
  int dummyCounter = 0 ;    

  //Take the electron track
  reco::GsfTrackRef tmpTrack = electron_->gsfTrack() ;
  math::XYZVector tmpElectronMomentumAtVtx = (*tmpTrack).innerMomentum () ; 

  for ( reco::TrackCollection::const_iterator itrTr  = (*trackCollection_).begin() ; 
                                              itrTr != (*trackCollection_).end()   ; 
	   			              ++itrTr ) 
    {
	math::XYZVector tmpTrackMomentumAtVtx = (*itrTr).innerMomentum () ; 
	double this_pt  = sqrt( tmpTrackMomentumAtVtx.Perp2 () );
	if ( this_pt < ptLow_ ) 
	  continue ;  
	if (fabs( (*itrTr).dz() - (*tmpTrack).dz() ) > lip_ )
          continue ;
	double dr = DeltaR(tmpTrackMomentumAtVtx,tmpElectronMomentumAtVtx) ;
	if ( fabs(dr) < extRadius_ && 
	     fabs(dr) > intRadius_ )
	  ++dummyCounter ;       
    }//end loop over tracks                 

  return dummyCounter ;
}

double TkIsolation::getPtTracks () const
{
  //dummy counter for the pT of tracks inside the cone
  double dummypT = 0 ;
  
  //Take the electron track
  reco::GsfTrackRef tmpTrack = electron_->gsfTrack() ;
  math::XYZVector tmpElectronMomentumAtVtx = (*tmpTrack).innerMomentum () ; 

  std::cout << "# tracks " << (*trackCollection_).size() << std::endl ;

  for ( reco::TrackCollection::const_iterator itrTr  = (*trackCollection_).begin() ; 
                                              itrTr != (*trackCollection_).end()   ; 
	   			              ++itrTr) 
    {
	math::XYZVector tmpTrackMomentumAtVtx = (*itrTr).innerMomentum () ; 
	double this_pt  = sqrt( tmpTrackMomentumAtVtx.Perp2 () );
	double dr = DeltaR(tmpTrackMomentumAtVtx,tmpElectronMomentumAtVtx) ;
	if ( this_pt < ptLow_ ) 
	  continue ;  
	if (fabs( (*itrTr).dz() - (*tmpTrack).dz() ) > lip_ )
          continue ;
	if ( fabs(dr) < extRadius_ && 
	     fabs(dr) > intRadius_ )
             dummypT += this_pt ;
    }

   return dummypT ;
}

