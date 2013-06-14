import FWCore.ParameterSet.Config as cms

runOnAOD = 1
prescale = 1
runReReco = True
if (runReReco == True):
    print "FIRST RERECO IS DONE"
applyLaserCalibration = True
print "USING PRESCALE FACTOR "+str(prescale)
if (runReReco == True and applyLaserCalibration == False):
    print "BEWARE!! RUNNING WITH NO ECAL LASER CALIBRATIONS"

process = cms.Process('reRECO')
#process.GlobalTag.globaltag = 'GR_R_42_V19::All'
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'GR_P_V20::All'
process.GlobalTag.globaltag = 'GR_R_42_V19::All'

print "USING GT "+str(process.GlobalTag.globaltag)
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("EcalIntercalibConstantsRcd"),
             tag = cms.string("EcalIntercalibConstants_v10_offline"),
             connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_ECAL")
             )
    ,cms.PSet(record = cms.string("EcalADCToGeVConstantRcd"),
              tag = cms.string("EcalADCToGeVConstant_v10_offline"),
              connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_ECAL")
              )
    ,cms.PSet(record = cms.string("EcalLaserAPDPNRatiosRcd"),
              #           tag = cms.string("EcalLaserAPDPNRatios_2011fit_noVPT_nolim_online"),
#              tag = cms.string("EcalLaserAPDPNRatios_test_20110625"),
#              tag = cms.string("EcalLaserAPDPNRatios_2011V3_online"),
              tag = cms.string("EcalLaserAPDPNRatios_prompt"),
              connect = cms.untracked.string("frontier://PromptProd/CMS_COND_311X_ECAL_LAS")
          )
 #beam spot to arrive to very last runs after 167151
    ,cms.PSet(record = cms.string("BeamSpotObjectsRcd"),
              tag = cms.string("BeamSpotObjects_PCL_byLumi_v0_prompt"),
              connect = cms.untracked.string("frontier://PromptProd/CMS_COND_31X_BEAMSPOT")
              )

)

process.prescaler = cms.EDFilter("Prescaler",
                                    prescaleFactor = cms.int32(prescale),
                                    prescaleOffset = cms.int32(0)
                                    )
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')



# --- jet met sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.metProducerSequence_cff")
#process.load("HiggsAnalysis.HiggsToWW2e.btagProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagPFJetsProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagPFPUcorrJetsProducerSequence_cff")
process.load("HiggsAnalysis.HiggsToWW2e.btagPFNoPUJetsProducerSequence_cff")
#process.load("HiggsAnalysis.HiggsToWW2e.btagJPTJetsProducerSequence_cff")
process.load("WWAnalysis.Tools.chargedMetProducer_cfi")
process.chargedMetProducer.collectionTag = "particleFlow"
process.chargedMetProducer.vertexTag = "offlinePrimaryVertices"

process.load("HiggsAnalysis.HiggsToWW2e.jetProducerSequenceFastJet_cff")
#process.newJetTracksAssociatorAtVertex.jets = 'ak5CaloJetsL2L3Residual'
#process.newSoftElectronTagInfos.jets = 'ak5CaloJetsL2L3Residual'
#process.newSoftMuonTagInfos.jets = 'ak5CaloJetsL2L3Residual'
process.newPFPUcorrJetTracksAssociatorAtVertex.jets = 'ak5PFJetsL1FastL2L3Residual'
process.newPFPUcorrJetsSoftElectronTagInfos.jets = 'ak5PFJetsL1FastL2L3Residual'
process.newPFPUcorrJetsSoftMuonTagInfos.jets = 'ak5PFJetsL1FastL2L3Residual'
#process.newPFJetNoPUTracksAssociatorAtVertex.jets = 'ak5PFJetsNoPUL1FastL2L3Residual'
#process.newPFJetsNoPUSoftElectronTagInfos.jets = 'ak5PFJetsNoPUL1FastL2L3Residual'
#process.newPFJetsNoPUSoftMuonTagInfos.jets = 'ak5PFJetsNoPUL1FastL2L3Residual'
#process.newJetTracksAssociatorAtVertex.jets = 'ak5CaloJetsL1FastL2L3'
#process.newSoftElectronTagInfos.jets = 'ak5CaloJetsL1FastL2L3'
#process.newSoftMuonTagInfos.jets = 'ak5CaloJetsL1FastL2L3'

# to correct calo met ---
#process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5CaloJet
#process.metMuonJESCorAK5 = metJESCorAK5CaloJet.clone()
#process.metMuonJESCorAK5.inputUncorMetLabel = "met"
#process.metMuonJESCorAK5.useTypeII = True
#process.metMuonJESCorAK5.hasMuonsCorr = False

# --- track sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.leptonLinkedTracks_cfi")

# --- electron sequences ---
process.load("HiggsAnalysis.HiggsToWW2e.electronIdSequence_cff")

# --- pf isolation sequence ---
process.load("HiggsAnalysis.HiggsToWW2e.leptonPFIsoSequence_cff")

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
process.treeDumper.dumpSCs = True
process.treeDumper.dumpBCs = False
process.treeDumper.dumpVertices = True
process.treeDumper.dumpCaloTowers = False
process.treeDumper.dumpParticleFlowObjects = False
process.treeDumper.dumpPFCandidates = False

process.treeDumper.dumpElectrons = cms.untracked.bool(True)
process.treeDumper.dumpCalibratedElectrons = cms.untracked.bool(False)
process.treeDumper.dumpPFlowElectrons = cms.untracked.bool(True)
process.treeDumper.dumpPFpreId = cms.untracked.bool(False)
process.treeDumper.dumpPhotons = cms.untracked.bool(True)
process.treeDumper.dumpMuons = cms.untracked.bool(False)
process.treeDumper.dumpPFTaus = cms.untracked.bool(False)
process.treeDumper.dumphpsPFTaus = cms.untracked.bool(False)
process.treeDumper.dumphpsTancTaus = cms.untracked.bool(False)
process.treeDumper.dumpPFCandidates = cms.untracked.bool(False)
process.treeDumper.dumpTracks = cms.untracked.bool(False)
process.treeDumper.dumpGsfTracks = cms.untracked.bool(True)
process.treeDumper.dumpMuonTracks = cms.untracked.bool(False)

process.treeDumper.dumpK0s = cms.untracked.bool(False)
process.treeDumper.dumpCaloTowers = cms.untracked.bool(False)
process.treeDumper.dumpJets = cms.untracked.bool(True)
process.treeDumper.dumpPUcorrPFJet = cms.untracked.bool(True)
process.treeDumper.dumpMet = cms.untracked.bool(True)

process.treeDumper.saveEcal = cms.untracked.bool(True)
process.treeDumper.saveFatEcal = cms.untracked.bool(True)
process.treeDumper.saveTrk = cms.untracked.bool(True)
process.treeDumper.saveFatTrk = cms.untracked.bool(False)
process.treeDumper.saveTrackDeDx = cms.untracked.bool(False)
process.treeDumper.saveEleID = cms.untracked.bool(True)
process.treeDumper.savePFEleGsfTrk = cms.untracked.bool(True)
process.treeDumper.savePFEleBasic = cms.untracked.bool(True)
process.treeDumper.saveJetBTag = cms.untracked.bool(True)
process.treeDumper.savePFTauBasic = cms.untracked.bool(False)
process.treeDumper.saveLeadPFCand = cms.untracked.bool(False)
process.treeDumper.savePFTauDiscriminators = cms.untracked.bool(False)



process.treeDumper.dumpTree = True
if (runOnAOD == 1) :
    process.treeDumper.saveFatTrk = False
    process.treeDumper.saveTrackDeDx = False
    process.treeDumper.dumpPFlowElectrons = False
    process.treeDumper.dumpHcalNoiseFlags = False
    process.treeDumper.AODHcalNoiseFlags = False
else :
    process.treeDumper.saveFatTrk = True
    process.treeDumper.saveTrackDeDx = True
    process.treeDumper.dumpPFlowElectrons = True
    process.treeDumper.dumpHcalNoiseFlags = True
    process.treeDumper.AODHcalNoiseFlags = False

process.options = cms.untracked.PSet(
      fileMode =  cms.untracked.string('NOMERGE')
      )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(300) )

process.source = cms.Source("PoolSource",
#                            noEventSort = cms.untracked.bool(True),
#                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
#                            skipEvents = cms.untracked.uint32(6764),
                            fileNames = cms.untracked.vstring('/store/data/Run2011A/DoubleElectron/RAW/v1/000/170/854/AE4FD824-4CB3-E011-A3DD-BCAEC53296F2.root')
                            )
if (runReReco == True):
    process.source.inputCommands = cms.untracked.vstring("drop *", "keep *_*_*_HLT")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

########################################
#no Laser Correction applied
########################################
if (applyLaserCalibration == False):
    process.ecalRecHit.laserCorrection = cms.bool(False)
# Path and EndPath definitions
if (runReReco == True):
    process.raw2digi_step = cms.Path(process.prescaler * process.RawToDigi)
    process.L1Reco_step = cms.Path(process.prescaler * process.L1Reco)
    process.reconstruction_step = cms.Path(process.prescaler * process.reconstruction)

process.vecbosReco = cms.Path ( process.prescaler *
                                process.leptonLinkedTracks
                                * process.mergedSuperClusters 
                                * process.chargedMetProducer
                                * process.metSequence
                                * process.pfIsolationAllSequence
                                * process.ourJetSequenceDataReduced
                                * process.newPFPUcorrJetBtaggingSequence
                                * process.eIdSequence
                                * process.FastjetForIsolation
                       )
#process.endjob_step = cms.Path(process.prescaler *process.endOfProcess)
#process.RECOoutput_step = cms.EndPath(process.RECOoutput)
process.RECOoutput_step = cms.Path(process.prescaler *process.treeDumper)
#process.q = cms.EndPath (  )
# Schedule definition
if (runReReco == True):
    process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.vecbosReco,process.RECOoutput_step)
else :
    process.schedule = cms.Schedule(process.vecbosReco,process.RECOoutput_step)

