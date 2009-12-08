import FWCore.ParameterSet.Config as cms

process = cms.Process("KsPi0Pi0Analysis")

process.extend(cms.include("RecoEcal/EgammaClusterProducers/data/geometryForClustering.cff"))
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR09_P_V6::All'

# --- electron sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.basicClusterMerger_cfi")

# --- tree dumper ---
process.load("HiggsAnalysis.HiggsToWW2e.treeDumper_cfi")
process.treeDumper.nameFile = 'default_data.root'
process.treeDumper.dumpTriggerResults = True
process.treeDumper.dumpGenInfo = False
process.treeDumper.dumpMCTruth = False
process.treeDumper.dumpSignalKfactor = False
process.treeDumper.dumpGenMet = False
process.treeDumper.dumpGenJets = False
process.treeDumper.dumpSCs = False
process.treeDumper.dumpBCs = True
process.treeDumper.dumpTracks = True
process.treeDumper.dumpGsfTracks = False
process.treeDumper.dumpVertices = True
process.treeDumper.dumpParticleFlowObjects = False
process.treeDumper.saveFatTrk = False
process.treeDumper.saveTrackDeDx = True
process.treeDumper.saveJet1BTag = False
process.treeDumper.saveJet2BTag = False
process.treeDumper.dumpElectrons = False
process.treeDumper.dumpPFlowElectrons = False
process.treeDumper.dumpMuons = False
process.treeDumper.dumpJets = False
process.treeDumper.dumpMet = False
process.treeDumper.dumpTree = True

process.options = cms.untracked.PSet(
      fileMode =  cms.untracked.string('NOMERGE')
      )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(300) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            debugFlag = cms.untracked.bool(True),
                            debugVebosity = cms.untracked.uint32(10),
                            fileNames = cms.untracked.vstring(
'file:/cmsrm/pc21_2/meridian/data/BeamCommissioning09_BSCFilter_v3/bscFilter_123596_7.root'    
                            )
                            )

process.p = cms.Path ( process.mergedBasicClusters )                       
process.q = cms.EndPath ( process.treeDumper )
