import FWCore.ParameterSet.Config as cms

process = cms.Process("HiggsToWW2e")

process.extend(cms.include("RecoEcal/EgammaClusterProducers/data/geometryForClustering.cff"))

# --- common preselection code ---
process.load("HiggsAnalysis.HiggsToWW2Leptons.HWWPreselectionSequence_cff")

# --- common code to comput gg fusion signal k factor
process.load("HiggsAnalysis.HiggsToWW2Leptons.HWWKFactorProducer_cfi")

# --- jet met sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.metProducerSequence_cff")

# --- electron sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.ambiguityResolvedElectrons_cfi")
process.load("HiggsAnalysis.HiggsToWW2e.alternativeIsolationSequence_cff")

# --- tree dumper ---
process.load("HiggsAnalysis.HiggsToWW2e.treeDumper_cfi")
process.treeDumper.nameFile = 'default.root'
process.treeDumper.dumpTriggerResults = True
process.treeDumper.dumpPreselInfo = True
process.treeDumper.dumpGenInfo = False
process.treeDumper.dumpSignalKfactor = True
process.treeDumper.dumpGenInfoMcAtNlo = False
process.treeDumper.dumpParticleFlowObjects = True
process.treeDumper.saveFatTrk = True
process.treeDumper.dumpTree = True

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                            debugFlag = cms.untracked.bool(True),
                            debugVebosity = cms.untracked.uint32(10),
                            fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/e/emanuele/RECO/relval_Zee_CMSSW_2_1_7.root')
                            )

process.p = cms.Path ( process.KFactorProducer *
                       process.higgsToWW2LeptonsPreselectionSequence *
                       process.jetSequence *
                       process.metSequence *
                       # process.alternativeIsolationSequence * # needed only when running on AODs
                       process.ambiguityResolvedElectrons )

process.q = cms.EndPath ( process.treeDumper )