import FWCore.ParameterSet.Config as cms

process = cms.Process("VecBosAnalysis")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR10_P_V11::All'

# --- jet met sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.metProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagPFJetsProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagJPTJetsProducerSequence_cff")

# to correct calo met ---
#process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5CaloJet
#process.metMuonJESCorAK5 = metJESCorAK5CaloJet.clone()
#process.metMuonJESCorAK5.inputUncorMetLabel = "met"
#process.metMuonJESCorAK5.useTypeII = True
#process.metMuonJESCorAK5.hasMuonsCorr = False

# --- electron sequences ---
process.load("RecoEgamma.EgammaIsolationAlgos.eleIsolationSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.superClusterMerger_cfi")
process.load("HiggsAnalysis.HiggsToWW2e.basicClusterMerger_cfi")
process.load("HiggsAnalysis.HiggsToWW2e.ambiguityResolvedElectrons_cfi")
process.load("HiggsAnalysis.HiggsToWW2e.electronIdSequence_cff")

# --- calotowers sequence ---
process.load("HiggsAnalysis.HiggsToWW2e.lowThrCaloTowers_cfi")

# --- ECAL clusters merging in a unique collection ---
process.load("HiggsAnalysis.HiggsToWW2e.superClusterMerger_cfi")

# --- tree dumper ---
process.load("HiggsAnalysis.HiggsToWW2e.treeDumper_cfi")
process.treeDumper.nameFile = 'default_data.root'
process.treeDumper.dumpTriggerResults = True
process.treeDumper.dumpHLTObjects = True
process.treeDumper.dumpGenInfo = False
process.treeDumper.dumpMCTruth = False
process.treeDumper.dumpSignalKfactor = False
process.treeDumper.dumpGenMet = False
process.treeDumper.dumpGenJets = False
process.treeDumper.dumpTracks = True
process.treeDumper.dumpGsfTracks = True
process.treeDumper.dumpSCs = True
process.treeDumper.dumpBCs = True
process.treeDumper.dumpVertices = True
process.treeDumper.dumpCaloTowers = False
process.treeDumper.dumpParticleFlowObjects = True
process.treeDumper.saveFatTrk = True
process.treeDumper.saveTrackDeDx = True
process.treeDumper.dumpTree = True

process.options = cms.untracked.PSet(
      fileMode =  cms.untracked.string('NOMERGE')
      )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(300) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
#                            skipEvents = cms.untracked.uint32(6764),
                            fileNames = cms.untracked.vstring('/store/data/Run2010A/EG/RECO/Jul16thReReco-v2/0063/10753EB3-2591-DF11-822A-002618943800.root')
                            )

process.p = cms.Path ( process.mergedSuperClusters * process.mergedBasicClusters *
                       process.ourJetSequence *
                       process.newBtaggingSequence * process.newPFJetBtaggingSequence * process.newJPTJetBtaggingSequence *
                       process.eIdSequence *
                       process.eleIsolationSequence *
                       process.ambiguityResolvedElectrons )
                       
process.q = cms.EndPath ( process.treeDumper )
