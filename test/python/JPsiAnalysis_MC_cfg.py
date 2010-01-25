import FWCore.ParameterSet.Config as cms

process = cms.Process("JPsiAnalysis")

process.extend(cms.include("RecoEcal/EgammaClusterProducers/data/geometryForClustering.cff"))
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'STARTUP3X_V8H::All'

# --- electron sequences ---
process.load("RecoEgamma.EgammaIsolationAlgos.eleIsolationSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.superClusterMerger_cfi")
process.load("HiggsAnalysis.HiggsToWW2e.basicClusterMerger_cfi")
process.load("HiggsAnalysis.HiggsToWW2e.ambiguityResolvedElectrons_cfi")
process.load("HiggsAnalysis.HiggsToWW2e.electronIdSequence_cff")

# --- tree dumper ---
process.load("HiggsAnalysis.HiggsToWW2e.treeDumper_cfi")
process.treeDumper.nameFile = 'default_MC.root'
process.treeDumper.dumpTriggerResults = False
process.treeDumper.dumpGenInfo = True
process.treeDumper.dumpMCTruth = True
process.treeDumper.dumpSignalKfactor = False
process.treeDumper.dumpSCs = True
process.treeDumper.dumpBCs = False
process.treeDumper.dumpTracks = True
process.treeDumper.dumpGsfTracks = True
process.treeDumper.dumpVertices = True
process.treeDumper.dumpParticleFlowObjects = True
process.treeDumper.saveFatTrk = True
process.treeDumper.saveTrackDeDx = True
process.treeDumper.saveJet1BTag = False
process.treeDumper.saveJet2BTag = False
process.treeDumper.dumpElectrons = True
process.treeDumper.dumpPFlowElectrons = True
process.treeDumper.dumpPFpreId = False
process.treeDumper.dumpMuons = False
process.treeDumper.dumpJets = False
process.treeDumper.dumpMet = False
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
                            fileNames = cms.untracked.vstring('file:/cmsrm/pc21/emanuele/data/Pool/jpsiEE_CMSSW_3_3_6_v6_97.root')
                            )

process.p = cms.Path ( process.mergedSuperClusters * 
                       process.eIdSequence *
                       process.eleIsolationSequence *
                       process.ambiguityResolvedElectrons )
                       
process.q = cms.EndPath ( process.treeDumper )
