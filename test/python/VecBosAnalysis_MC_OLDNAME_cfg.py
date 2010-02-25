import FWCore.ParameterSet.Config as cms

process = cms.Process("VecBosAnalysis")

process.extend(cms.include("RecoEcal/EgammaClusterProducers/data/geometryForClustering.cff"))
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_31X_V8::All'

# --- jet met sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.metProducerSequence_cff")

# --- electron sequences ---
process.load("RecoEgamma.EgammaIsolationAlgos.eleIsolationSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.ambiguityResolvedElectrons_cfi")
process.load("HiggsAnalysis.HiggsToWW2e.electronIdSequence_cff")

# --- calotowers sequence ---
process.load("HiggsAnalysis.HiggsToWW2e.lowThrCaloTowers_cfi")

# --- ECAL clusters merging in a unique collection ---
process.load("HiggsAnalysis.HiggsToWW2e.superClusterMerger_cfi")

# --- tree dumper ---
process.load("HiggsAnalysis.HiggsToWW2e.treeDumper_cfi")
process.treeDumper.nameFile = 'default_MC.root'
process.treeDumper.dumpTriggerResults = True
process.treeDumper.dumpGenInfo = True
process.treeDumper.dumpSignalKfactor = False
process.treeDumper.dumpTracks = True
process.treeDumper.dumpGsfTracks = True
process.treeDumper.dumpSCs = True
process.treeDumper.dumpVertices = True
process.treeDumper.dumpCaloTowers = True
process.treeDumper.dumpGenJets = True
process.treeDumper.dumpParticleFlowObjects = True
process.treeDumper.saveFatTrk = True
process.treeDumper.saveTrackDeDx = True
process.treeDumper.saveJet1BTag = True
process.treeDumper.saveJet2BTag = False
process.treeDumper.dumpTree = True

# old ak5 naming
######################################################
process.treeDumper.genJetCollection = "antikt5GenJets"
process.treeDumper.jetCollection2 = "antikt5CaloJets"
process.treeDumper.PFjetCollection2 = "antikt5PFJets"
process.jetVertexAlpha2.JET_ALGO = "antikt5CaloJets"
process.L2L3CorJetAK5Calo.src = "antikt5CaloJets"
process.L2L3CorJetAK5PF.src = "antikt5PFJets"
######################################################

process.options = cms.untracked.PSet(
      fileMode =  cms.untracked.string('NOMERGE')
      )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            debugFlag = cms.untracked.bool(True),
                            debugVebosity = cms.untracked.uint32(10),
                            fileNames = cms.untracked.vstring('file:/cmsrm/pc23_2/emanuele/Pool/9AC9CBF4-DB97-DE11-AAD4-001E0B5A645A.root')
                            )

process.p = cms.Path ( process.lowThrCaloTowers * process.mergedSuperClusters *
                       process.ourJetSequence * process.newBtaggingSequence *
                       process.eIdSequence *
                       process.eleIsolationSequence *
                       process.ambiguityResolvedElectrons )
                       
process.q = cms.EndPath ( process.treeDumper )
