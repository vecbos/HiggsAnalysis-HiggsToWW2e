import FWCore.ParameterSet.Config as cms

process = cms.Process("VecBosAnalysis")

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START3X_V26A::All'

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
process.load("HiggsAnalysis.HiggsToWW2e.basicClusterMerger_cfi")
process.load("HiggsAnalysis.HiggsToWW2e.superClusterMerger_cfi")

# --- to recover the ak5 GenJets that are not re-recoed in 33X samples ---
process.load('RecoJets.Configuration.GenJetParticles_cff')
process.load("RecoJets.Configuration.RecoGenJets_cff")

# --- tree dumper ---
process.load("HiggsAnalysis.HiggsToWW2e.treeDumper_cfi")
process.treeDumper.nameFile = 'default_MC.root'
process.treeDumper.dumpTriggerResults = True
process.treeDumper.dumpGenInfo = True
process.treeDumper.dumpSignalKfactor = False
process.treeDumper.dumpTracks = True
process.treeDumper.dumpElectrons = True
process.treeDumper.dumpGsfTracks = True
process.treeDumper.dumpSCs = True
process.treeDumper.dumpBCs = True
process.treeDumper.dumpVertices = True
process.treeDumper.dumpCaloTowers = False
process.treeDumper.dumpGenJets = True
process.treeDumper.dumpParticleFlowObjects = True
process.treeDumper.saveFatTrk = True
process.treeDumper.saveTrackDeDx = True
process.treeDumper.saveJet1BTag = True
process.treeDumper.dumpTree = True

# in order to take the correct HLT informations in Spring10 samples
process.treeDumper.TriggerResultsTag = "TriggerResults::REDIGI"

process.options = cms.untracked.PSet(
      fileMode =  cms.untracked.string('NOMERGE')
      )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = cms.untracked.vstring('/store/mc/Spring10/QCD_BCtoE_Pt80to170/GEN-SIM-RECO/MC_3XY_V25_S09_preproduction-v2/0120/08901658-0F3B-DF11-8E12-001D0967DA17.root')
                            )

process.p = cms.Path ( process.lowThrCaloTowers * process.mergedBasicClusters * process.mergedSuperClusters *
                       process.genParticlesForJets * process.ak5GenJets * # added for re-recoed V9 Summer09 samples where the ak5GenJet collection was dropped
                       process.ourJetSequence * process.newBtaggingSequence *
                       process.gsfElectrons *
                       process.eIdSequence *
                       process.eleIsolationSequence *
                       process.ambiguityResolvedElectrons )
                       
process.q = cms.EndPath ( process.treeDumper )
