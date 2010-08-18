import FWCore.ParameterSet.Config as cms

# b-tagging general configuration
from RecoJets.JetAssociationProducers.ak5JTA_cff import *
from RecoBTag.Configuration.RecoBTag_cff import *

# create a new jets and tracks association
import RecoJets.JetAssociationProducers.ak5JTA_cff
newJPTJetTracksAssociatorAtVertex = RecoJets.JetAssociationProducers.ak5JTA_cff.ak5JetTracksAssociatorAtVertex.clone()
newJPTJetTracksAssociatorAtVertex.jets = "ak5JPTJetsL2L3"
newJPTJetTracksAssociatorAtVertex.tracks = "generalTracks"

# impact parameter b-tag
import RecoBTag.Configuration.RecoBTag_cff
newJPTJetsImpactParameterTagInfos = RecoBTag.Configuration.RecoBTag_cff.impactParameterTagInfos.clone()
newJPTJetsImpactParameterTagInfos.jetTracks = "newJPTJetTracksAssociatorAtVertex"
newTrackCountingHighEffBJPTJetTags = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighEffBJetTags.clone()
newTrackCountingHighEffBJPTJetTags.tagInfos = cms.VInputTag( cms.InputTag("newJPTJetsImpactParameterTagInfos") )
newTrackCountingHighPurBJPTJetTags = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighPurBJetTags.clone()
newTrackCountingHighPurBJPTJetTags.tagInfos = cms.VInputTag( cms.InputTag("newJPTJetsImpactParameterTagInfos") )
newJetProbabilityBJPTJetTags = RecoBTag.Configuration.RecoBTag_cff.jetProbabilityBJetTags.clone()
newJetProbabilityBJPTJetTags.tagInfos = cms.VInputTag( cms.InputTag("newJPTJetsImpactParameterTagInfos") )
newJetBProbabilityBJPTJetTags = RecoBTag.Configuration.RecoBTag_cff.jetBProbabilityBJetTags.clone()
newJetBProbabilityBJPTJetTags.tagInfos = cms.VInputTag( cms.InputTag("newJPTJetsImpactParameterTagInfos") )
newImpactParameterMVABJPTJetTags = RecoBTag.Configuration.RecoBTag_cff.impactParameterMVABJetTags.clone()
newImpactParameterMVABJPTJetTags.tagInfos = cms.VInputTag( cms.InputTag("newJPTJetsImpactParameterTagInfos") )

# secondary vertex b-tag
newJPTJetsSecondaryVertexTagInfos = RecoBTag.Configuration.RecoBTag_cff.secondaryVertexTagInfos.clone()
newJPTJetsSecondaryVertexTagInfos.trackIPTagInfos = "newJPTJetsImpactParameterTagInfos"
newSimpleSecondaryVertexBJPTJetTags = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexBJetTags.clone()
newSimpleSecondaryVertexBJPTJetTags.tagInfos = cms.VInputTag( cms.InputTag("newJPTJetsSecondaryVertexTagInfos") )
newCombinedSecondaryVertexBJPTJetTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexBJetTags.clone()
newCombinedSecondaryVertexBJPTJetTags.tagInfos = cms.VInputTag( cms.InputTag("newJPTJetsImpactParameterTagInfos"), cms.InputTag("newJPTJetsSecondaryVertexTagInfos") )
newCombinedSecondaryVertexMVABJPTJetTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexMVABJetTags.clone()
newCombinedSecondaryVertexMVABJPTJetTags.tagInfos = cms.VInputTag( cms.InputTag("newJPTJetsImpactParameterTagInfos"), cms.InputTag("newJPTJetsSecondaryVertexTagInfos") )

# soft electron b-tag
newJPTJetsSoftElectronTagInfos = RecoBTag.Configuration.RecoBTag_cff.softElectronTagInfos.clone()
newJPTJetsSoftElectronTagInfos.jets = "ak5JPTJetsL2L3"
newSoftElectronBJPTJetTags = RecoBTag.Configuration.RecoBTag_cff.softElectronBJetTags.clone()
newSoftElectronBJPTJetTags.tagInfos = cms.VInputTag( cms.InputTag("newJPTJetsSoftElectronTagInfos") )
newSoftElectronByIP3dBJPTJetTags = RecoBTag.Configuration.RecoBTag_cff.softElectronByIP3dBJetTags.clone()
newSoftElectronByIP3dBJPTJetTags.tagInfos = cms.VInputTag( cms.InputTag("newJPTJetsSoftElectronTagInfos") )
newSoftElectronByPtBJPTJetTags = RecoBTag.Configuration.RecoBTag_cff.softElectronByPtBJetTags.clone()
newSoftElectronByPtBJPTJetTags.tagInfos = cms.VInputTag( cms.InputTag("newJPTJetsSoftElectronTagInfos") )

# soft muon b-tag
newJPTJetsSoftMuonTagInfos = RecoBTag.Configuration.RecoBTag_cff.softMuonTagInfos.clone()
newJPTJetsSoftMuonTagInfos.jets = "ak5JPTJetsL2L3"
newSoftMuonBJPTJetTags = RecoBTag.Configuration.RecoBTag_cff.softMuonBJetTags.clone()
newSoftMuonBJPTJetTags.tagInfos = cms.VInputTag( cms.InputTag("newJPTJetsSoftMuonTagInfos") )
newSoftMuonByIP3dBJPTJetTags = RecoBTag.Configuration.RecoBTag_cff.softMuonByIP3dBJetTags.clone()
newSoftMuonByIP3dBJPTJetTags.tagInfos = cms.VInputTag( cms.InputTag("newJPTJetsSoftMuonTagInfos") )
newSoftMuonByPtBJPTJetTags = RecoBTag.Configuration.RecoBTag_cff.softMuonByPtBJetTags.clone()
newSoftMuonByPtBJPTJetTags.tagInfos = cms.VInputTag( cms.InputTag("newJPTJetsSoftMuonTagInfos") )

# prepare a path running the new modules
newJPTJetTracksAssociator = cms.Sequence(
    newJPTJetTracksAssociatorAtVertex
    )

newJPTJetBtaggingIP = cms.Sequence(
    newJPTJetsImpactParameterTagInfos * (
       newTrackCountingHighEffBJPTJetTags +
       newTrackCountingHighPurBJPTJetTags +
       newJetProbabilityBJPTJetTags +
       newJetBProbabilityBJPTJetTags )
    )

newJPTJetBtaggingSV = cms.Sequence(
    newJPTJetsImpactParameterTagInfos *
    newJPTJetsSecondaryVertexTagInfos * (
       newSimpleSecondaryVertexBJPTJetTags +
       newCombinedSecondaryVertexBJPTJetTags +
       newCombinedSecondaryVertexMVABJPTJetTags )
    )

newJPTJetBtaggingEle = cms.Sequence(
    softElectronCands * # already run for calojets
    newJPTJetsSoftElectronTagInfos *
       newSoftElectronBJPTJetTags +
       newSoftElectronByIP3dBJPTJetTags +
       newSoftElectronByPtBJPTJetTags
    )

newJPTJetBtaggingMu = cms.Sequence(
    newJPTJetsSoftMuonTagInfos * (
       newSoftMuonBJPTJetTags +
       newSoftMuonByIP3dBJPTJetTags +
       newSoftMuonByPtBJPTJetTags )
    )

newJPTJetBtagging = cms.Sequence(
    newJPTJetBtaggingIP +
    newJPTJetBtaggingSV +
    newJPTJetBtaggingEle +
    newJPTJetBtaggingMu )

newJPTJetBtaggingSequence = cms.Sequence(
    newJPTJetTracksAssociator *
       newJPTJetBtagging )
