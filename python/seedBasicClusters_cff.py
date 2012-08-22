import FWCore.ParameterSet.Config as cms

# first merge the EB/EE basic clusters in one collection
mergedBasicClusters = cms.EDProducer("BasicClusterCombiner",
                                     labels = cms.VInputTag(cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters"),
                                                            cms.InputTag("multi5x5SuperClusters","multi5x5EndcapBasicClusters"))
                                     )

# this is for ECAL clusters
seedBasicClusters = cms.EDProducer("SuperClusterSeedsProducer",
                                   BasicClusterLabel = cms.InputTag("mergedBasicClusters"),
                                   SuperClusterLabel = cms.InputTag("electronAndPhotonSuperClusters")
                                   )


seedBasicClustersSequence = cms.Sequence( mergedBasicClusters
                                          * seedBasicClusters )
