#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "RecoMuon/TrackerSeedGenerator/plugins/CollectionCombiner.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWEleAmbiguityResolve.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/LeptonTrackFilter.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/LeptonCandidateFilter.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/SuperClusterSeedsProducer.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/ElectronAndPhotonSuperClusterProducer.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/ChargedPFMetProducer.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWTreeDumper.h"


typedef ObjectSelector< HWWEleAmbiguityResolve, 
			reco::GsfElectronRefVector
			> AmbResolver ;

typedef ObjectSelector< LeptonTrackFilter, 
			reco::TrackRefVector
			> leptonTrackFilter ;

typedef ObjectSelector< LeptonCandidateFilter, 
			reco::CandidateCollection
			> leptonCandidateFilter ;

typedef CollectionCombiner< reco::SuperClusterCollection > SuperClusterCombiner;
typedef CollectionCombiner< reco::BasicClusterCollection > BasicClusterCombiner;
typedef CollectionCombiner< reco::PFMETCollection > PfMetCombiner;

DEFINE_FWK_MODULE (AmbResolver) ;
DEFINE_FWK_MODULE (leptonTrackFilter);
DEFINE_FWK_MODULE (leptonCandidateFilter);
DEFINE_FWK_MODULE (SuperClusterSeedsProducer);
DEFINE_FWK_MODULE (ElectronAndPhotonSuperClusterProducer);
DEFINE_FWK_MODULE (ChargedPFMetProducer);
DEFINE_FWK_MODULE(BasicClusterCombiner);
DEFINE_FWK_MODULE(SuperClusterCombiner);
DEFINE_FWK_MODULE(PfMetCombiner);
DEFINE_FWK_MODULE (HWWTreeDumper);  
