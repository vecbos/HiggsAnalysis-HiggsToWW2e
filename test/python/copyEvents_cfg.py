from FWCore.ParameterSet.Config import *

process = Process("copyEvents")

process.include( "FWCore/MessageLogger/data/MessageLogger.cfi" )

process.maxEvents = untracked.PSet( input = untracked.int32(300) )
process.source = Source( "PoolSource",
                         fileNames = untracked.vstring( "/store/relval/CMSSW_2_1_7/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/0C3B40D7-F87D-DD11-A9FB-000423D998BA.root" )
                         )

process.out = OutputModule( "PoolOutputModule",
                            fileName = untracked.string( "relvalZee217.root" ),
                            outputCommands = untracked.vstring("keep * ")
                            )

process.o = EndPath( process.out )
