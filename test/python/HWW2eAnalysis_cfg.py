import FWCore.ParameterSet.Config as cms

process = cms.Process("HiggsToWW2e")

process.extend(cms.include("RecoEcal/EgammaClusterProducers/data/geometryForClustering.cff"))

# --- common preselection code ---
# TODO????

# --- common code to comput gg fusion signal k factor
# TODO????

# --- jet met sequence ---
process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequence_cff")

# --- electron sequence ---
process.load("HiggsAnalysis.HiggsToWW2e.ambiguityResolvedElectrons_cfi")
process.ambiguityResolvedElectrons.reducedElectronsRefCollectionLabel = "pixelMatchGsfElectronsRef" # temporary until the common code is not migrated
process.ambiguityResolvedElectrons.doRefCheck = False

process.load("HiggsAnalysis.HiggsToWW2e.electronIdSequence_cff")

# --- tree dumper ---
process.load("HiggsAnalysis.HiggsToWW2e.treeDumper_cfi")
process.treeDumper.nameFile = 'default.root'
process.treeDumper.dumpTriggerResults = True
process.treeDumper.dumpPreselInfo = False
process.treeDumper.dumpGenInfo = False
process.treeDumper.dumpSignalKfactor = False
process.treeDumper.dumpGenInfoMcAtNlo = False
process.treeDumper.saveFatTrk = False
process.treeDumper.dumpTree = True

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                            debugFlag = cms.untracked.bool(True),
                            debugVebosity = cms.untracked.uint32(10),
                            fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/e/emanuele/RECO/relval_Zee_CMSSW_2_1_7.root')
                            )

process.p = cms.Path ( process.jetSequence *
#                       higgsToWW2LeptonsPreselectionSequence *
                       process.ambiguityResolvedElectrons * process.eIdSequence )

process.q = cms.EndPath ( process.treeDumper )
