import FWCore.ParameterSet.Config as cms

process = cms.Process("HiggsToWW2e")

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_36Y_V4::All'

# --- common code to comput gg fusion signal k factor
process.load("HiggsAnalysis.HiggsToWW2Leptons.HWWKFactorProducer_cfi")

# --- apply the muon correction from the common code
process.load("HiggsAnalysis.HiggsToWW2Leptons.HWWMetCorrector_cfi")

# --- jet met sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.metProducerSequence_cff")

# --- electron sequences ---
process.load("RecoEgamma.EgammaIsolationAlgos.eleIsolationSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.ambiguityResolvedElectrons_cfi")
process.load("HiggsAnalysis.HiggsToWW2e.electronIdSequence_cff")

# --- ECAL clusters merging in a unique collection ---
process.load("HiggsAnalysis.HiggsToWW2e.superClusterMerger_cfi")

# --- tree dumper ---
process.load("HiggsAnalysis.HiggsToWW2e.treeDumper_cfi")
process.treeDumper.nameFile = 'default.root'
process.treeDumper.dumpTriggerResults = True
process.treeDumper.dumpGenInfo = True
process.treeDumper.dumpSignalKfactor = True
process.treeDumper.dumpTracks = True
process.treeDumper.dumpGsfTracks = True
process.treeDumper.dumpVertices = True
process.treeDumper.dumpParticleFlowObjects = True
process.treeDumper.dumpSCs = True
process.treeDumper.saveFatTrk = True
process.treeDumper.saveJet1BTag = True
process.treeDumper.saveJet2BTag = False
process.treeDumper.dumpTree = True

process.options = cms.untracked.PSet(
      fileMode =  cms.untracked.string('NOMERGE')
      )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = cms.untracked.vstring('/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V26A_356ReReco-v1/0009/FEFC70B6-F53D-DF11-B57E-003048679150.root')
                            )

process.p = cms.Path ( process.KFactorProducer * process.muonCorrectedMET *
                       process.mergedSuperClusters *
                       process.ourJetSequence * process.newBtaggingSequence *
                       process.eIdSequence *
                       process.eleIsolationSequence *
                       process.ambiguityResolvedElectrons)
                       
process.q = cms.EndPath ( process.treeDumper )
