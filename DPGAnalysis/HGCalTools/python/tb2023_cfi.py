import FWCore.ParameterSet.Config as cms


def addPerformanceReports(process, addMemCheck=False):

    # add timing and mem (too slow) for FWK jobs report
    process.Timing = cms.Service("Timing",
                                 summaryOnly=cms.untracked.bool(True),
                                 useJobReport=cms.untracked.bool(True))

    if addMemCheck:
        process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
                                                ignoreTotal=cms.untracked.int32(1),
                                                jobReportOutputOnly=cms.untracked.bool(True))

    return process


def configTBConditions(process, key='default', moduleType='LD'):
    """ maybe this should be done with eras/modifiers? """

    module_locators = {
        'LD': 'Geometry/HGCalMapping/data/modulelocator_tb.txt',
        'HD': 'Geometry/HGCalMapping/data/modulelocator_tb_HD.txt',
        'LD3': 'Geometry/HGCalMapping/data/modulelocator_tb_LD3.txt',
    }
    process.hgCalModuleInfoESSource.filename = module_locators[moduleType]
    process.hgCalSiModuleInfoESSource.filename = 'Geometry/HGCalMapping/data/WaferCellMapTraces.txt'

    hex_templates = {
        'LD': '/eos/cms/store/group/dpg_hgcal/comm_hgcal/ykao/hexagons_20230801.root',
        'HD': '/eos/cms/store/group/dpg_hgcal/comm_hgcal/ykao/geometry_HD_full_wafer.root',
        'LD3': '/eos/cms/store/group/dpg_hgcal/comm_hgcal/ykao/geometry_LD3_partial_wafer_20230913.root',
    }
    if hasattr(process, 'hgCalDigisClientHarvester'):
        process.hgCalDigisClientHarvester.HexTemplateFile = hex_templates[moduleType]

    pedestals = {
        'default': '/eos/cms/store/group/dpg_hgcal/comm_hgcal/ykao/calibration_parameters_v2.txt',
    }

    if hasattr(process, 'hgCalPedestalsESSource'):
        process.hgCalPedestalsESSource.filename = pedestals[key]
    if hasattr(process, 'hgcalCalibESProducer'):
        process.hgcalCalibESProducer.filename = pedestals[key]
    # if hasattr(process,'hgcalConfigurationESProducer'):
    ###    process.hgcalConfigurationESProducer.filename = yamls[key]

    return process
