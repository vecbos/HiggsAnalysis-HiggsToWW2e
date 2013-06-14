import FWCore.ParameterSet.Config as cms

from EgammaAnalysis.ElectronTools.calibratedElectrons_cfi import *
#calibratedElectrons.isMC = cms.bool(False)                          
#calibratedElectrons.inputDataset = cms.string("22Jan2013ReReco")    
calibratedElectrons.isMC = cms.bool(True)
calibratedElectrons.inputDataset = cms.string("Summer12_LegacyPaper")  
calibratedElectrons.updateEnergyError = cms.bool(True)
calibratedElectrons.correctionsType = cms.int32(2)   
calibratedElectrons.combinationType = cms.int32(3)   
calibratedElectrons.lumiRatio = cms.double(1.0)      
calibratedElectrons.verbose = cms.bool(False)
calibratedElectrons.synchronization = cms.bool(False) 

from EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi import *
eleRegressionEnergy.inputElectronsTag = cms.InputTag('gsfElectrons')
eleRegressionEnergy.inputCollectionType = cms.uint32(0)
eleRegressionEnergy.energyRegressionType = cms.uint32(2)
eleRegressionEnergy.useRecHitCollections = cms.bool(True)
eleRegressionEnergy.produceValueMaps = cms.bool(True)

eCalibSequence = cms.Sequence( eleRegressionEnergy * calibratedElectrons)

