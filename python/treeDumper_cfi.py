import FWCore.ParameterSet.Config as cms

treeDumper = cms.EDFilter("HWWTreeDumper",
                         TriggerResultsTag = cms.InputTag("TriggerResults::HLT"),
                          electronCollection = cms.InputTag("ambiguityResolvedElectrons"),
                          muonCollection = cms.InputTag("isolatedMuons"), # preselection
                          ecalBarrelSCCollection = cms.InputTag("correctedHybridSuperClusters"),
                          ecalEndcapSCCollection = cms.InputTag("multi5x5SuperClusters","multi5x5EndcapSuperClusters"),
                          ecalBarrelRecHits = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
                          ecalEndcapRecHits = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
                          electronIdCutsLabel = cms.untracked.string('egammaIDCutsLoose'),
                          electronIdLikelihoodLabel = cms.untracked.string('egammaIDLikelihood'),
                          tracksForIsolationProducer = cms.InputTag("generalTracks"),
                          calotowersForIsolationProducer = cms.InputTag("towerMaker"),
                          jetCollection1 = cms.InputTag("iterativeCone5CaloJets"),
                          genJetCollection1 = cms.InputTag("iterativeCone5GenJets"),
                          jetVertexAlphaCollection1 = cms.InputTag("jetVertexAlpha1","Var"),
                          jetCollection2 = cms.InputTag("sisCone5CaloJets"),
                          genJetCollection2 = cms.InputTag("sisCone5GenJets"),
                          jetVertexAlphaCollection2 = cms.InputTag("jetVertexAlpha2","Var"),
                          PFjetCollection1 = cms.InputTag("iterativeCone5PFJets"),
                          PFjetCollection2 = cms.InputTag("sisCone5PFJets"),
                          metCollection = cms.InputTag("muonCorrectedMET"), # preselection
                          genMetCollection = cms.InputTag("genMet"),
                          PFmetCollection = cms.InputTag("pfMET"),
                          hepMcCollection = cms.InputTag("source"),
                          genInfoCollection = cms.InputTag("source"),
                          genWeightCollection = cms.untracked.string('CSA07WeightProducer'),
                          nameFile = cms.untracked.string('analysisTree.root'),
                          nameTree = cms.untracked.string('ntp1'),
                          # switch ON/OFF the candidate collections to dump
                          dumpElectrons = cms.untracked.bool(True),
                          dumpMuons = cms.untracked.bool(True),
                          dumpSCs = cms.untracked.bool(False),
                          dumpJets = cms.untracked.bool(True),
                          dumpMet = cms.untracked.bool(True),
                          # switch ON/OFF the particle flow objects to dump
                          dumpParticleFlowObjects = cms.untracked.bool(False),
                          # switch ON/OFF the additional informations to dump
                          saveEcal = cms.untracked.bool(True),
                          saveFatEcal = cms.untracked.bool(True),
                          saveTrk = cms.untracked.bool(True),
                          saveFatTrk = cms.untracked.bool(True),
                          saveEleID = cms.untracked.bool(True),
                          saveJetAlpha = cms.untracked.bool(True),
                          # MC truth
                          mcTruthCollection = cms.InputTag("genParticles"),
                          electronMatchMap = cms.InputTag("electronMatchMap"),
                          muonMatchMap = cms.InputTag("muonMatchMap"),
                          dumpGenMet = cms.untracked.bool(True),
                          dumpGenJets = cms.untracked.bool(True),
                          dumpMCTruth = cms.untracked.bool(True),
                          dumpGenInfo = cms.untracked.bool(True),
                          dumpPreselInfo = cms.untracked.bool(True),
                          dumpSignalKfactor = cms.untracked.bool(True),
                          dumpGenInfoMcAtNlo = cms.untracked.bool(False),
                          doMCEleMatch = cms.untracked.bool(False),
                          doMCMuonMatch = cms.untracked.bool(False),
                          # trigger results
                          dumpTriggerResults = cms.untracked.bool(False),
                          # effectively dump the data into the tree
                          dumpTree = cms.untracked.bool(False)
                          )
