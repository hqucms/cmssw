import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("TESTDQM")

options = VarParsing('analysis')
options.register('moduleType',
                 'LD',
                 VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "type of the module, e.g., LD, HD, LD3, etc.")
options.parseArguments()


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 50000
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(options.inputFiles))

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

# Logical mapping
process.load('Geometry.HGCalMapping.hgCalModuleInfoESSource_cfi')
process.load('Geometry.HGCalMapping.hgCalSiModuleInfoESSource_cfi')

# HGCal DQM
process.load('DQM.HGCal.hgCalDigisClient_cfi')
process.load('DQM.HGCal.hgCalDigisClientHarvester_cfi')
process.hgCalDigisClient.Prescale = 1000

process.DQMStore = cms.Service("DQMStore")
process.load("DQMServices.FileIO.DQMFileSaverOnline_cfi")
process.dqmSaver.tag = 'HGCAL'

# path
process.p = cms.Path(process.hgCalDigisClient * process.hgCalDigisClientHarvester * process.dqmSaver)

# configure test beam conditions
from DPGAnalysis.HGCalTools.tb2023_cfi import configTBConditions, addPerformanceReports
configTBConditions(process, moduleType=options.moduleType)
addPerformanceReports(process)
