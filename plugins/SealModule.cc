#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "RecoMuon/TrackerSeedGenerator/plugins/CollectionCombiner.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWEleAmbiguityResolve.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWTreeDumper.h"


typedef ObjectSelector< HWWEleAmbiguityResolve, 
			reco::GsfElectronRefVector
			> AmbResolver ;

typedef CollectionCombiner< reco::SuperClusterCollection > SuperClusterCombiner;
typedef CollectionCombiner< reco::BasicClusterCollection > BasicClusterCombiner;

DEFINE_FWK_MODULE (AmbResolver) ;
DEFINE_FWK_MODULE(SuperClusterCombiner);
DEFINE_FWK_MODULE(BasicClusterCombiner);
DEFINE_FWK_MODULE (HWWTreeDumper);  
