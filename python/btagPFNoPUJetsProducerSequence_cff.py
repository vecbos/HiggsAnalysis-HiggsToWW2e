import FWCore.ParameterSet.Config as cms

# b-tagging general configuration
from RecoJets.JetAssociationProducers.ak5JTA_cff import *
from RecoBTag.Configuration.RecoBTag_cff import *

# create a new jets and tracks association
import RecoJets.JetAssociationProducers.ak5JTA_cff
newPFJetNoPUTracksAssociatorAtVertex = RecoJets.JetAssociationProducers.ak5JTA_cff.ak5JetTracksAssociatorAtVertex.clone()
newPFJetNoPUTracksAssociatorAtVertex.jets = "ak5PFJetsNoPUL1FastL2L3Residual"
newPFJetNoPUTracksAssociatorAtVertex.tracks = "generalTracks"

# impact parameter b-tag
import RecoBTag.Configuration.RecoBTag_cff
newPFJetsNoPUImpactParameterTagInfos = RecoBTag.Configuration.RecoBTag_cff.impactParameterTagInfos.clone()
newPFJetsNoPUImpactParameterTagInfos.jetTracks = "newPFJetNoPUTracksAssociatorAtVertex"
newTrackCountingHighEffBPFJetNoPUTags = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighEffBJetTags.clone()
newTrackCountingHighEffBPFJetNoPUTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsNoPUImpactParameterTagInfos") )
newTrackCountingHighPurBPFJetNoPUTags = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighPurBJetTags.clone()
newTrackCountingHighPurBPFJetNoPUTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsNoPUImpactParameterTagInfos") )
newJetProbabilityBPFJetNoPUTags = RecoBTag.Configuration.RecoBTag_cff.jetProbabilityBJetTags.clone()
newJetProbabilityBPFJetNoPUTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsNoPUImpactParameterTagInfos") )
newJetBProbabilityBPFJetNoPUTags = RecoBTag.Configuration.RecoBTag_cff.jetBProbabilityBJetTags.clone()
newJetBProbabilityBPFJetNoPUTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsNoPUImpactParameterTagInfos") )
newImpactParameterMVABPFJetNoPUTags = RecoBTag.Configuration.RecoBTag_cff.impactParameterMVABJetTags.clone()
newImpactParameterMVABPFJetNoPUTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsNoPUImpactParameterTagInfos") )

# secondary vertex b-tag
newPFJetsNoPUSecondaryVertexTagInfos = RecoBTag.Configuration.RecoBTag_cff.secondaryVertexTagInfos.clone()
newPFJetsNoPUSecondaryVertexTagInfos.trackIPTagInfos = "newPFJetsNoPUImpactParameterTagInfos"
newSimpleSecondaryVertexHighEffBPFJetNoPUTags = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexHighEffBJetTags.clone()
newSimpleSecondaryVertexHighEffBPFJetNoPUTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsNoPUSecondaryVertexTagInfos") )
newSimpleSecondaryVertexHighPurBPFJetNoPUTags = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexHighPurBJetTags.clone()
newSimpleSecondaryVertexHighPurBPFJetNoPUTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsNoPUSecondaryVertexTagInfos") )
newCombinedSecondaryVertexBPFJetNoPUTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexBJetTags.clone()
newCombinedSecondaryVertexBPFJetNoPUTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsNoPUImpactParameterTagInfos"), cms.InputTag("newPFJetsNoPUSecondaryVertexTagInfos") )
newCombinedSecondaryVertexMVABPFJetNoPUTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexMVABJetTags.clone()
newCombinedSecondaryVertexMVABPFJetNoPUTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsNoPUImpactParameterTagInfos"), cms.InputTag("newPFJetsNoPUSecondaryVertexTagInfos") )

# soft electron b-tag
newPFJetsNoPUSoftElectronTagInfos = RecoBTag.Configuration.RecoBTag_cff.softElectronTagInfos.clone()
newPFJetsNoPUSoftElectronTagInfos.jets = "ak5PFJetsL2L3Residual"
newSoftElectronBPFJetNoPUTags = RecoBTag.Configuration.RecoBTag_cff.softElectronBJetTags.clone()
newSoftElectronBPFJetNoPUTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsNoPUSoftElectronTagInfos") )
newSoftElectronByIP3dBPFJetNoPUTags = RecoBTag.Configuration.RecoBTag_cff.softElectronByIP3dBJetTags.clone()
newSoftElectronByIP3dBPFJetNoPUTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsNoPUSoftElectronTagInfos") )
newSoftElectronByPtBPFJetNoPUTags = RecoBTag.Configuration.RecoBTag_cff.softElectronByPtBJetTags.clone()
newSoftElectronByPtBPFJetNoPUTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsNoPUSoftElectronTagInfos") )

# soft muon b-tag
newPFJetsNoPUSoftMuonTagInfos = RecoBTag.Configuration.RecoBTag_cff.softMuonTagInfos.clone()
newPFJetsNoPUSoftMuonTagInfos.jets = "ak5PFJetsL2L3Residual"
newSoftMuonBPFJetNoPUTags = RecoBTag.Configuration.RecoBTag_cff.softMuonBJetTags.clone()
newSoftMuonBPFJetNoPUTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsNoPUSoftMuonTagInfos") )
newSoftMuonByIP3dBPFJetNoPUTags = RecoBTag.Configuration.RecoBTag_cff.softMuonByIP3dBJetTags.clone()
newSoftMuonByIP3dBPFJetNoPUTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsNoPUSoftMuonTagInfos") )
newSoftMuonByPtBPFJetNoPUTags = RecoBTag.Configuration.RecoBTag_cff.softMuonByPtBJetTags.clone()
newSoftMuonByPtBPFJetNoPUTags.tagInfos = cms.VInputTag( cms.InputTag("newPFJetsNoPUSoftMuonTagInfos") )

# prepare a path running the new modules
newPFJetNoPUTracksAssociator = cms.Sequence(
    newPFJetNoPUTracksAssociatorAtVertex
    )

newPFJetNoPUBtaggingIP = cms.Sequence(
    newPFJetsNoPUImpactParameterTagInfos * (
       newTrackCountingHighEffBPFJetNoPUTags +
       newTrackCountingHighPurBPFJetNoPUTags +
       newJetProbabilityBPFJetNoPUTags +
       newJetBProbabilityBPFJetNoPUTags )
    )

newPFJetNoPUBtaggingSV = cms.Sequence(
    newPFJetsNoPUImpactParameterTagInfos *
    newPFJetsNoPUSecondaryVertexTagInfos * (
       newSimpleSecondaryVertexHighEffBPFJetNoPUTags +
       newSimpleSecondaryVertexHighPurBPFJetNoPUTags +
       newCombinedSecondaryVertexBPFJetNoPUTags +
       newCombinedSecondaryVertexMVABPFJetNoPUTags )
    )

newPFJetNoPUBtaggingEle = cms.Sequence(
    softElectronCands * # already run for calojets
    newPFJetsNoPUSoftElectronTagInfos *
       newSoftElectronBPFJetNoPUTags +
       newSoftElectronByIP3dBPFJetNoPUTags +
       newSoftElectronByPtBPFJetNoPUTags
    )

newPFJetNoPUBtaggingMu = cms.Sequence(
    newPFJetsNoPUSoftMuonTagInfos * (
       newSoftMuonBPFJetNoPUTags +
       newSoftMuonByIP3dBPFJetNoPUTags +
       newSoftMuonByPtBPFJetNoPUTags )
    )

newPFJetNoPUBtagging = cms.Sequence(
    newPFJetNoPUBtaggingIP +
    newPFJetNoPUBtaggingSV +
    newPFJetNoPUBtaggingEle +
    newPFJetNoPUBtaggingMu )

newPFJetNoPUBtaggingSequence = cms.Sequence(
    newPFJetNoPUTracksAssociator *
       newPFJetNoPUBtagging )
