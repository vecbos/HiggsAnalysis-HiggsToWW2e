import FWCore.ParameterSet.Config as cms

# Select good muons to remove from cone (skim definition. To choose electrons in the overlap removal. Not perfect)
goodMuons = cms.EDFilter("MuonRefSelector",
                         src = cms.InputTag("muons"),
                         cut = cms.string("pt > 10 && " +
                                          "isGlobalMuon && isTrackerMuon"
                                          )
                         )

# Select good electrons to remove from cone (skim definition. To choose electrons in the overlap removal. Not perfect)
goodElectrons = cms.EDFilter("GsfElectronRefSelector",
                             src = cms.InputTag("gsfElectrons"),
                             cut = cms.string( "pt > 10 &&" +
                                               " abs(deltaEtaSuperClusterTrackAtVtx) < 0.010 &&" +
                                               " (( isEB && sigmaIetaIeta < 0.011) ||" +
                                               "  (!isEB && sigmaIetaIeta < 0.031))")
                             )

# create isolation 'deposits'
pfIsoChargedHadronsDZ = cms.EDFilter("PdgIdPFCandidateSelector",
                                     pdgId = cms.vint32(211, -211, 321, -321, 999211, 2212, -2212),
                                     src = cms.InputTag("particleFlow")
                                     )
pfIsoNeutralHadrons = cms.EDFilter("PdgIdPFCandidateSelector",
                                   pdgId = cms.vint32(111, 130, 310, 2112),
                                   src = cms.InputTag("pfNoPileUp")
                                   )
pfIsoChargedHadrons = pfIsoNeutralHadrons.clone ( pdgId = [211, -211, 321, -321, 999211, 2212, -2212] )
pfIsoPhotons        = pfIsoNeutralHadrons.clone ( pdgId = [22] )

# make the actual IsoDeposits
isoDepMuonWithChargedIsoDZ = cms.EDProducer("CandIsoDepositProducer",
                                            src = cms.InputTag("muons"),
                                            MultipleDepositsFlag = cms.bool(False),
                                            trackType = cms.string('candidate'),
                                            ExtractorPSet = cms.PSet(
    Diff_z = cms.double(0.1),
    ComponentName = cms.string('CandViewExtractor'),
    DR_Max = cms.double(1.0),
    Diff_r = cms.double(99999.99),
    inputCandView = cms.InputTag("pfIsoChargedHadronsDZ"),
    DR_Veto = cms.double(1e-05),
    DepositLabel = cms.untracked.string('')
    )
)

isoDepMuonWithChargedIso = cms.EDProducer("CandIsoDepositProducer",
                                          src = cms.InputTag("muons"),
                                          MultipleDepositsFlag = cms.bool(False),
                                          trackType = cms.string('candidate'),
                                          ExtractorPSet = cms.PSet(
    Diff_z = cms.double(99999.99),
    ComponentName = cms.string('CandViewExtractor'),
    DR_Max = cms.double(1.0),
    Diff_r = cms.double(99999.99),
    inputCandView = cms.InputTag("pfIsoChargedHadrons"),
    DR_Veto = cms.double(1e-05),
    DepositLabel = cms.untracked.string('')
    )
)

isoDepMuonWithNeutralIso = isoDepMuonWithChargedIso.clone()
isoDepMuonWithPhotonIso = isoDepMuonWithChargedIso.clone()
isoDepMuonWithNeutralIso.ExtractorPSet.inputCandView = "pfIsoNeutralHadrons"
isoDepMuonWithPhotonIso.ExtractorPSet.inputCandView = "pfIsoPhotons"

isoDepElectronWithChargedIsoDZ = isoDepMuonWithChargedIsoDZ.clone( src = "gsfElectrons" )
isoDepElectronWithChargedIso = isoDepMuonWithChargedIso.clone( src = "gsfElectrons" )
isoDepElectronWithNeutralIso = isoDepMuonWithNeutralIso.clone( src = "gsfElectrons" )
isoDepElectronWithPhotonIso  = isoDepMuonWithPhotonIso.clone( src = "gsfElectrons" )

# Mu vetos
muVetos       = cms.vstring('0.01')
elChargedVeto = cms.vstring('0.01')
elNeutralVeto = cms.vstring('0.07')
elPhotonVeto  = cms.vstring('RectangularEtaPhiVeto(-0.025,0.025,-0.5,0.5)')

# insert them into the pat leptons
# ha, made you look, they are actually down below in the electron and muon sections

# make the crazy sequence
pfIsoSequence = cms.Sequence(
    pfIsoNeutralHadrons +
    pfIsoChargedHadrons +
    pfIsoChargedHadronsDZ +
    pfIsoPhotons * (
    isoDepMuonWithChargedIsoDZ +
    isoDepMuonWithChargedIso +
    isoDepMuonWithNeutralIso +
    isoDepMuonWithPhotonIso +
    isoDepElectronWithChargedIsoDZ +
    isoDepElectronWithChargedIso +
    isoDepElectronWithNeutralIso +
    isoDepElectronWithPhotonIso
    )
)


# add the pf isolation values
# muons
# ucsd guess
from RecoEgamma.EgammaIsolationAlgos.eleIsoFromDepsModules_cff import *
muonPfChargedDeps = eleIsoFromDepsHcalFromTowers.clone()
muonPfChargedDeps.deposits[-1].src = "isoDepMuonWithChargedIso"
muonPfChargedDeps.deposits[-1].label = cms.string("pfCharged")
muonPfChargedDeps.deposits[-1].deltaR = 0.4
muonPfChargedDeps.deposits[-1].vetos  = [ 'Threshold(0.5)' ]
muonPfChargedDeps.deposits[-1].vetos += [ veto for veto in muVetos ]
muonPfChargedDeps.deposits[-1].vetos += [ 'goodMuons:'+veto for veto in muVetos ]
muonPfChargedDeps.deposits[-1].vetos += [ 'goodElectrons:'+veto for veto in elChargedVeto ]

muonPfNeutralDeps = eleIsoFromDepsHcalFromTowers.clone()
muonPfNeutralDeps.deposits[-1].src = "isoDepMuonWithNeutralIso"
muonPfNeutralDeps.deposits[-1].label = cms.string("pfNeutral")
muonPfNeutralDeps.deposits[-1].deltaR = 0.4
muonPfNeutralDeps.deposits[-1].vetos  = [ 'Threshold(0.5)' ]
muonPfNeutralDeps.deposits[-1].vetos += [ veto for veto in muVetos ]
muonPfNeutralDeps.deposits[-1].vetos += [ 'goodMuons:'+veto for veto in muVetos ]
muonPfNeutralDeps.deposits[-1].vetos += [ 'goodElectrons:'+veto for veto in elNeutralVeto ]

muonPfPhotonDeps = eleIsoFromDepsHcalFromTowers.clone()
muonPfPhotonDeps.deposits[-1].src = "isoDepMuonWithPhotonIso"
muonPfPhotonDeps.deposits[-1].label = cms.string("pfPhoton")
muonPfPhotonDeps.deposits[-1].deltaR = 0.4
muonPfPhotonDeps.deposits[-1].vetos  = [ 'Threshold(0.5)' ]
muonPfPhotonDeps.deposits[-1].vetos += [ veto for veto in muVetos ]
muonPfPhotonDeps.deposits[-1].vetos += [ 'goodMuons:'+veto for veto in muVetos ]
muonPfPhotonDeps.deposits[-1].vetos += [ 'goodElectrons:'+veto for veto in elPhotonVeto ]



electronPfChargedDeps = eleIsoFromDepsHcalFromTowers.clone()
electronPfChargedDeps.deposits[-1].src = "isoDepElectronWithChargedIso"
electronPfChargedDeps.deposits[-1].label = cms.string("pfCharged")
electronPfChargedDeps.deposits[-1].deltaR = 0.4
electronPfChargedDeps.deposits[-1].vetos  = [ 'Threshold(0.5)' ]
electronPfChargedDeps.deposits[-1].vetos += [ veto for veto in muVetos ]
electronPfChargedDeps.deposits[-1].vetos += [ 'goodMuons:'+veto for veto in muVetos ]
electronPfChargedDeps.deposits[-1].vetos += [ 'goodElectrons:'+veto for veto in elChargedVeto ]

electronPfNeutralDeps = eleIsoFromDepsHcalFromTowers.clone()
electronPfNeutralDeps.deposits[-1].src = "isoDepElectronWithNeutralIso"
electronPfNeutralDeps.deposits[-1].label = cms.string("pfNeutral")
electronPfNeutralDeps.deposits[-1].deltaR = 0.4
electronPfNeutralDeps.deposits[-1].vetos  = [ 'Threshold(0.5)' ]
electronPfNeutralDeps.deposits[-1].vetos += [ veto for veto in muVetos ]
electronPfNeutralDeps.deposits[-1].vetos += [ 'goodMuons:'+veto for veto in muVetos ]
electronPfNeutralDeps.deposits[-1].vetos += [ 'goodElectrons:'+veto for veto in elNeutralVeto ]

electronPfPhotonDeps = eleIsoFromDepsHcalFromTowers.clone()
electronPfPhotonDeps.deposits[-1].src = "isoDepElectronWithPhotonIso"
electronPfPhotonDeps.deposits[-1].label = cms.string("pfPhoton")
electronPfPhotonDeps.deposits[-1].deltaR = 0.4
electronPfPhotonDeps.deposits[-1].vetos  = [ 'Threshold(0.5)' ]
electronPfPhotonDeps.deposits[-1].vetos += [ veto for veto in muVetos ]
electronPfPhotonDeps.deposits[-1].vetos += [ 'goodMuons:'+veto for veto in muVetos ]
electronPfPhotonDeps.deposits[-1].vetos += [ 'goodElectrons:'+veto for veto in elPhotonVeto ]




# muons
# smurf default
muonPfGenericChargedDeps = eleIsoFromDepsHcalFromTowers.clone()
muonPfGenericChargedDeps.deposits[-1].src = "isoDepMuonWithChargedIsoDZ"
muonPfGenericChargedDeps.deposits[-1].label = cms.string("pfGenericCharged")
muonPfGenericChargedDeps.deposits[-1].deltaR = 0.3
muonPfGenericChargedDeps.deposits[-1].vetos = [ '0.01' ]

muonPfGenericNeutralDeps = eleIsoFromDepsHcalFromTowers.clone()
muonPfGenericNeutralDeps.deposits[-1].src = "isoDepMuonWithNeutralIso"
muonPfGenericNeutralDeps.deposits[-1].label = cms.string("pfGenericNeutral")
muonPfGenericNeutralDeps.deposits[-1].deltaR = 0.3
muonPfGenericNeutralDeps.deposits[-1].vetos  = [ 'Threshold(1.0)' ]

muonPfGenericPhotonDeps = eleIsoFromDepsHcalFromTowers.clone()
muonPfGenericPhotonDeps.deposits[-1].src = "isoDepMuonWithPhotonIso"
muonPfGenericPhotonDeps.deposits[-1].label = cms.string("pfGenericPhoton")
muonPfGenericPhotonDeps.deposits[-1].deltaR = 0.3
muonPfGenericPhotonDeps.deposits[-1].vetos = [ ]

# electrons
# smurf default
electronPfGenericChargedDeps = eleIsoFromDepsHcalFromTowers.clone()
electronPfGenericChargedDeps.deposits[-1].src = "isoDepElectronWithChargedIsoDZ"
electronPfGenericChargedDeps.deposits[-1].label = cms.string("pfGenericCharged")
electronPfGenericChargedDeps.deposits[-1].deltaR = 0.4
electronPfGenericChargedDeps.deposits[-1].vetos  = [ '0.01' ]

electronPfGenericNeutralDeps = eleIsoFromDepsHcalFromTowers.clone()
electronPfGenericNeutralDeps.deposits[-1].src = "isoDepElectronWithNeutralIso"
electronPfGenericNeutralDeps.deposits[-1].label = cms.string("pfGenericNeutral")
electronPfGenericNeutralDeps.deposits[-1].deltaR = 0.4
electronPfGenericNeutralDeps.deposits[-1].vetos  = [ 'Threshold(1.0)','0.07' ]

electronPfGenericPhotonDeps = eleIsoFromDepsHcalFromTowers.clone()
electronPfGenericPhotonDeps.deposits[-1].src = "isoDepElectronWithPhotonIso"
electronPfGenericPhotonDeps.deposits[-1].label = cms.string("pfGenericPhoton")
electronPfGenericPhotonDeps.deposits[-1].deltaR = 0.4
electronPfGenericPhotonDeps.deposits[-1].vetos  = [ 'RectangularEtaPhiVeto(-0.025,0.025,-0.5,0.5)' ]







# muons
# smurf default no overlap
muonPfGenericNoOverChargedDeps = eleIsoFromDepsHcalFromTowers.clone()
muonPfGenericNoOverChargedDeps.deposits[-1].src = "isoDepMuonWithChargedIsoDZ"
muonPfGenericNoOverChargedDeps.deposits[-1].label = cms.string("pfGenericNoOverCharged")
muonPfGenericNoOverChargedDeps.deposits[-1].deltaR = 0.3
muonPfGenericNoOverChargedDeps.deposits[-1].vetos = [ '0.01' ]
muonPfGenericNoOverChargedDeps.deposits[-1].vetos += [ 'goodMuons:0.01','goodElectrons:0.01' ]

muonPfGenericNoOverNeutralDeps = eleIsoFromDepsHcalFromTowers.clone()
muonPfGenericNoOverNeutralDeps.deposits[-1].src = "isoDepMuonWithNeutralIso"
muonPfGenericNoOverNeutralDeps.deposits[-1].label = cms.string("pfGenericNoOverNeutral")
muonPfGenericNoOverNeutralDeps.deposits[-1].deltaR = 0.3
muonPfGenericNoOverNeutralDeps.deposits[-1].vetos  = [ 'Threshold(1.0)' ]
muonPfGenericNoOverNeutralDeps.deposits[-1].vetos += [ 'goodMuons:0.01','goodElectrons:0.01' ]

muonPfGenericNoOverPhotonDeps = eleIsoFromDepsHcalFromTowers.clone()
muonPfGenericNoOverPhotonDeps.deposits[-1].src = "isoDepMuonWithPhotonIso"
muonPfGenericNoOverPhotonDeps.deposits[-1].label = cms.string("pfGenericNoOverPhoton")
muonPfGenericNoOverPhotonDeps.deposits[-1].deltaR = 0.3
muonPfGenericNoOverPhotonDeps.deposits[-1].vetos = [ ]
muonPfGenericNoOverPhotonDeps.deposits[-1].vetos += [ 'goodMuons:0.01','goodElectrons:0.01' ]


# electrons
# smurf default no overlap
electronPfGenericNoOverChargedDeps = eleIsoFromDepsHcalFromTowers.clone()
electronPfGenericNoOverChargedDeps.deposits[-1].src = "isoDepElectronWithChargedIsoDZ"
electronPfGenericNoOverChargedDeps.deposits[-1].label = cms.string("pfGenericNoOverCharged")
electronPfGenericNoOverChargedDeps.deposits[-1].deltaR = 0.4
electronPfGenericNoOverChargedDeps.deposits[-1].vetos  = [ '0.01' ]
electronPfGenericNoOverChargedDeps.deposits[-1].vetos += [ 'goodMuons:0.01','goodElectrons:0.01' ]

electronPfGenericNoOverNeutralDeps = eleIsoFromDepsHcalFromTowers.clone()
electronPfGenericNoOverNeutralDeps.deposits[-1].src = "isoDepElectronWithNeutralIso"
electronPfGenericNoOverNeutralDeps.deposits[-1].label = cms.string("pfGenericNoOverNeutral")
electronPfGenericNoOverNeutralDeps.deposits[-1].deltaR = 0.4
electronPfGenericNoOverNeutralDeps.deposits[-1].vetos  = [ 'Threshold(1.0)','0.07' ]
electronPfGenericNoOverNeutralDeps.deposits[-1].vetos += [ 'goodMuons:0.01','goodElectrons:0.01' ]

electronPfGenericNoOverPhotonDeps = eleIsoFromDepsHcalFromTowers.clone()
electronPfGenericNoOverPhotonDeps.deposits[-1].src = "isoDepElectronWithPhotonIso"
electronPfGenericNoOverPhotonDeps.deposits[-1].label = cms.string("pfGenericNoOverPhoton")
electronPfGenericNoOverPhotonDeps.deposits[-1].deltaR = 0.4
electronPfGenericNoOverPhotonDeps.deposits[-1].vetos  = [ 'RectangularEtaPhiVeto(-0.025,0.025,-0.5,0.5)' ]
electronPfGenericNoOverPhotonDeps.deposits[-1].vetos += [ 'goodMuons:0.01','goodElectrons:0.01' ]

# this is the simplest PF isolation (combined)
import WWAnalysis.Tools.electronPFIsoMapProd_cfi
electronCombinedPFIsoMapProducer = WWAnalysis.Tools.electronPFIsoMapProd_cfi.electronPFIsoMapProd.clone()
electronCombinedPFIsoMapProducer.vtxLabel = 'offlinePrimaryVertices' # if the event has the first vertex bad, will be discarded offline.

import WWAnalysis.Tools.muonPFIsoMapProd_cfi
muonCombinedPFIsoMapProducer = WWAnalysis.Tools.muonPFIsoMapProd_cfi.muonPFIsoMapProd.clone()
muonCombinedPFIsoMapProducer.vtxLabel = 'offlinePrimaryVertices' # if the event has the first vertex bad, will be discarded offline.

# this is the PF candidate isolation with pfnopu input and custom vetoes
from MyAnalysis.IsolationTools.electronPFIsolations_cff import *
from MyAnalysis.IsolationTools.muonPFIsolations_cff import *

from CommonTools.ParticleFlow.pfNoPileUp_cff import *
pfPileUp.PFCandidates = "particleFlow"
pfNoPileUp.bottomCollection = "particleFlow"

pfPUSequence = cms.Sequence( pfPileUp * pfNoPileUp )

pfIsoStdSequence = cms.Sequence( goodMuons * goodElectrons * pfIsoSequence )

pfIsolationsDefault = cms.Sequence( muonPfChargedDeps * muonPfNeutralDeps * muonPfPhotonDeps *
                                    electronPfChargedDeps * electronPfNeutralDeps * electronPfPhotonDeps )

pfIsolationGenericDefault = cms.Sequence( muonPfGenericChargedDeps * muonPfGenericNeutralDeps * muonPfGenericPhotonDeps *
                                          electronPfGenericChargedDeps * electronPfGenericNeutralDeps * electronPfGenericPhotonDeps )

pfIsolationGenericNoOverlapDefault = cms.Sequence( muonPfGenericNoOverChargedDeps * muonPfGenericNoOverNeutralDeps * muonPfGenericNoOverPhotonDeps *
                                                   electronPfGenericNoOverChargedDeps * electronPfGenericNoOverNeutralDeps * electronPfGenericNoOverPhotonDeps )

pfIsolationCombined = cms.Sequence( electronCombinedPFIsoMapProducer * muonCombinedPFIsoMapProducer )

pfIsolationSingleType = cms.Sequence( electronPFIsoChHad * electronPFIsoNHad * electronPFIsoPhoton * muonPFIsoChHad * muonPFIsoNHad * muonPFIsoPhoton )

pfIsolationAllSequence = cms.Sequence( pfPUSequence * pfIsoStdSequence * pfIsolationsDefault * pfIsolationGenericDefault * pfIsolationGenericNoOverlapDefault * pfIsolationCombined * pfIsolationSingleType )

