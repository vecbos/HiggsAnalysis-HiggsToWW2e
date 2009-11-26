import FWCore.ParameterSet.Config as cms

process = cms.Process("HiggsToWW2e")

process.extend(cms.include("RecoEcal/EgammaClusterProducers/data/geometryForClustering.cff"))
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration/StandardSequences/RawToDigi_Data_cff")
process.load("Configuration.StandardSequences.ReconstructionCosmics_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR09_P_V6::All'

# run a dedicated ECAL reco sequence
process.load("RecoLocalCalo.Configuration.ecalLocalRecoSequenceCosmics_cff")
process.load("RecoEcal.EgammaClusterProducers.islandClusteringSequence_cff")

# not cosmics tracker reco
process.load("RecoTracker.Configuration.RecoTracker_cff")
process.load("TrackingTools.Configuration.TrackingTools_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
process.load("RecoVertex.Configuration.RecoVertex_cff")
process.load("RecoPixelVertexing.Configuration.RecoPixelVertexing_cff")

# --- apply the muon correction from the common code
process.load("HiggsAnalysis.HiggsToWW2Leptons.HWWMetCorrector_cfi")

# --- jet met sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.metProducerSequence_cff")

# --- electron sequences ---
process.load("RecoEgamma.EgammaIsolationAlgos.eleIsolationSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.ambiguityResolvedElectrons_cfi")
process.load("HiggsAnalysis.HiggsToWW2e.electronIdSequence_cff")

# --- track sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.trackCandidates_cfi")

# --- tree dumper ---
process.load("HiggsAnalysis.HiggsToWW2e.treeDumper_cfi")
process.treeDumper.ecalBarrelBCCollection = 'islandBasicClusters:islandBarrelBasicClusters'
process.treeDumper.ecalEndcapBCCollection = 'islandBasicClusters:islandEndcapBasicClusters'
process.treeDumper.metCollection = 'met'
process.treeDumper.nameFile = 'default.root'
process.treeDumper.dumpTriggerResults = True
process.treeDumper.dumpBCs = True
process.treeDumper.removeBadCrystalsInBCs = False
process.treeDumper.dumpTracks = True
process.treeDumper.dumpElectrons = False
process.treeDumper.dumpMuons = False
process.treeDumper.dumpJets = False
process.treeDumper.dumpMet = False
process.treeDumper.dumpVertices = False
process.treeDumper.dumpParticleFlowObjects = False
process.treeDumper.dumpGenInfo = False
process.treeDumper.dumpMCTruth = False
process.treeDumper.dumpGenMet = False
process.treeDumper.dumpGenJets = False
process.treeDumper.dumpSignalKfactor = False
process.treeDumper.saveFatTrk = True
process.treeDumper.saveJet1BTag = False
process.treeDumper.saveJet2BTag = False
process.treeDumper.dumpTree = True

process.options = cms.untracked.PSet(
      fileMode =  cms.untracked.string('NOMERGE')
      )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            debugFlag = cms.untracked.bool(True),
                            debugVebosity = cms.untracked.uint32(10),
                            fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/m/meridian/firstCollSkims/firstColl_bscFilter_122294.root',
                                                              #'rfio:/castor/cern.ch/user/m/meridian/firstCollSkims/firstColl_bscFilter_122314.root'
                                                              #'rfio:/castor/cern.ch/user/m/meridian/firstCollSkims/firstColl_bscFilter_122318.root'
                                                              )
                            )

process.digisSequence = cms.Sequence ( process.RawToDigi ) 

process.localrecoSequence = cms.Sequence ( process.trackerlocalreco * process.ecalLocalRecoSequence )

process.globalrecoSequence = cms.Sequence ( process.offlineBeamSpot+process.recopixelvertexing*process.ckftracks+process.islandClusteringSequence*process.islandSuperClusters+process.vertexreco)

process.p = cms.Path ( process.digisSequence *
                       process.localrecoSequence *
                       process.globalrecoSequence *
                       process.trackCandidates)
                       
process.q = cms.EndPath ( process.treeDumper )
