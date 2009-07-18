import FWCore.ParameterSet.Config as cms

process = cms.Process("HiggsToWW2e")

process.extend(cms.include("RecoEcal/EgammaClusterProducers/data/geometryForClustering.cff"))
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'IDEAL_V9::All'

# --- common preselection code ---
process.load("HiggsAnalysis.HiggsToWW2Leptons.HWWPreselectionSequence_cff")

# --- common code to comput gg fusion signal k factor
process.load("HiggsAnalysis.HiggsToWW2Leptons.HWWKFactorProducer_cfi")

# --- jet met sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.metProducerSequence_cff")

# --- electron sequences ---
process.load("RecoEgamma.EgammaIsolationAlgos.eleIsolationSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.ambiguityResolvedElectrons_cfi")
process.ambiguityResolvedElectrons.reducedElectronsRefCollectionLabel = "isolatedElectronsRef"
process.ambiguityResolvedElectrons.doRefCheck = True # i.e. take only the isolated electrons and resolve them
process.load("HiggsAnalysis.HiggsToWW2e.electronIdSequence_cff")
# the FastSim stores only the reduced rechits collection (the rh belonging to supercluster)
process.eleIsoDepositEcalFromHits.ExtractorPSet.barrelEcalHits = "reducedEcalRecHitsEB"
process.eleIsoDepositEcalFromHits.ExtractorPSet.endcapEcalHits = "reducedEcalRecHitsEE"
process.electronEcalRecHitIsolationLcone.ecalBarrelRecHitProducer = "reducedEcalRecHitsEB"
process.electronEcalRecHitIsolationLcone.ecalBarrelRecHitCollection = ""
process.electronEcalRecHitIsolationLcone.ecalEndcapRecHitProducer = "reducedEcalRecHitsEE"
process.electronEcalRecHitIsolationLcone.ecalEndcapRecHitCollection = ""
process.electronEcalRecHitIsolationScone.ecalBarrelRecHitProducer = "reducedEcalRecHitsEB"
process.electronEcalRecHitIsolationScone.ecalBarrelRecHitCollection = ""
process.electronEcalRecHitIsolationScone.ecalEndcapRecHitProducer = "reducedEcalRecHitsEE"
process.electronEcalRecHitIsolationScone.ecalEndcapRecHitCollection = ""


# --- track sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.trackCandidates_cfi")

# --- tree dumper ---
process.load("HiggsAnalysis.HiggsToWW2e.treeDumper_cfi")
process.treeDumper.nameFile = 'default.root'
process.treeDumper.dumpTriggerResults = True
process.treeDumper.dumpPreselInfo = True
process.treeDumper.dumpGenInfo = True
process.treeDumper.dumpSignalKfactor = True
process.treeDumper.dumpGenInfoMcAtNlo = False
process.treeDumper.dumpTracks = True
process.treeDumper.dumpVertices = True
process.treeDumper.dumpParticleFlowObjects = True
process.treeDumper.saveFatTrk = True
process.treeDumper.dumpTree = True
# the FastSim stores only the reduced rechits collection (the rh belonging to supercluster)
process.treeDumper.ecalBarrelRecHits = "reducedEcalRecHitsEB"
process.treeDumper.ecalEndcapRecHits = "reducedEcalRecHitsEE"
# the fast sim does not contain the track extra infos
process.treeDumper.saveFatEcal = False
process.treeDumper.saveFatTrk = False
# the btag is only run on IC jets. FastSim does not contain the rechits to re-run
process.treeDumper.saveJet1BTag = False
process.treeDumper.saveJet2BTag = False
# only the IC PF jets are in the RECO
process.treeDumper.PFjetCollection1 = 'L2L3CorJetIC5PF'
process.treeDumper.PFjetCollection2 = 'iterativeCone5PFJets'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            debugFlag = cms.untracked.bool(True),
                            debugVebosity = cms.untracked.uint32(10),
                            fileNames = cms.untracked.vstring('file:/cmsrm/pc21/emanuele/data/Higgs2.2.X/Wlnu_FastSim_1kEvents_1.root')
                            )

process.p = cms.Path ( process.KFactorProducer *
                       process.higgsToWW2LeptonsPreselectionSequence *
                       process.jetSequence *
                       process.pfjetICSequence *
                       process.eIdSequence *
                       process.eleIsoFromDepsHcalFromHits * process.egammaIsolationSequence * process.eleIsolationSequence *
                       process.ambiguityResolvedElectrons *
                       process.trackCandidates )

process.q = cms.EndPath ( process.treeDumper )
