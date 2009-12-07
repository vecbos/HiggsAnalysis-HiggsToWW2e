import FWCore.ParameterSet.Config as cms

mergedBasicClusters = cms.EDProducer("BasicClusterCombiner",
                                     labels = cms.VInputTag(cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters"),
                                                            cms.InputTag("multi5x5BasicClusters","multi5x5EndcapBasicClusters"))
                                     )
