import FWCore.ParameterSet.Config as cms

process = cms.Process("HiggsToWW2e")

process.extend(cms.include("RecoEcal/EgammaClusterProducers/data/geometryForClustering.cff"))
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR09_P_V6::All'

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

# --- track sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.trackCandidates_cfi")

# --- tree dumper ---
process.load("HiggsAnalysis.HiggsToWW2e.treeDumper_cfi")
process.treeDumper.nameFile = 'default.root'
process.treeDumper.dumpTriggerResults = True
process.treeDumper.dumpGenInfo = False
process.treeDumper.dumpMCTruth = False
process.treeDumper.dumpGenMet = False
process.treeDumper.dumpGenJets = False
process.treeDumper.dumpSignalKfactor = False
process.treeDumper.dumpTracks = True
process.treeDumper.dumpVertices = False
process.treeDumper.dumpParticleFlowObjects = False
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
                            fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/e/emanuele/Beams09/First11Events23Nov09.root')
                            )

process.p = cms.Path ( process.muonCorrectedMET *
                       process.jetSequence * process.pfjetSCSequence *
                       process.eIdSequence *
                       process.eleIsolationSequence *
                       process.ambiguityResolvedElectrons *
                       process.trackCandidates)
                       
process.q = cms.EndPath ( process.treeDumper )
