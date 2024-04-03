#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ModuleFactory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "CondFormats/DataRecord/interface/HGCalMappingModuleIndexerRcd.h"
#include "CondFormats/DataRecord/interface/HGCalMappingModuleRcd.h"
#include "CondFormats/HGCalObjects/interface/HGCalMappingModuleIndexer.h"
#include "CondFormats/HGCalObjects/interface/HGCalMappingParameterHost.h"
#include "CondFormats/HGCalObjects/interface/alpaka/HGCalMappingParameterDevice.h"
#include "DataFormats/HGCalDigi/interface/HGCalElectronicsId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace hgcal {

    class HGCalMappingModuleESProducer : public ESProducer {
    public:
      //
      HGCalMappingModuleESProducer(const edm::ParameterSet& iConfig)
          : ESProducer(iConfig), filename_(iConfig.getParameter<edm::FileInPath>("filename")) {
        auto cc = setWhatProduced(this);
        moduleIndexTkn_ = cc.consumes(iConfig.getParameter<edm::ESInputTag>("moduleindexer"));
      }

      //
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
        edm::ParameterSetDescription desc;
        desc.add<edm::FileInPath>("filename")->setComment("module locator file");
        desc.add<edm::ESInputTag>("moduleindexer", edm::ESInputTag(""))->setComment("Dense module index tool");
        descriptions.addWithDefaultLabel(desc);
      }

      //
      std::optional<HGCalMappingModuleParamHost> produce(const HGCalMappingModuleRcd& iRecord) {
        //get cell and module indexer
        auto modIndexer = iRecord.get(moduleIndexTkn_);

        // load dense indexing
        const uint32_t size = modIndexer.maxModulesIdx_;
        HGCalMappingModuleParamHost moduleParams(size, cms::alpakatools::host());
        for (size_t i = 0; i < size; i++)
          moduleParams.view()[i].valid() = false;

        // load module mapping parameters
        std::ifstream file(filename_.fullPath());
        std::string line, typecode;
        size_t iline(0);
        int plane, i1, i2, zside, celltype;
        uint16_t fedid, slinkidx, captureblock, econdidx, captureblockidx, typeidx;
        bool isSiPM;
        uint32_t idx, eleid, detid(0);

        while (std::getline(file, line)) {
          iline++;
          if (iline == 1)
            continue;

          std::istringstream stream(line);
          stream >> plane >> i1 >> i2 >> typecode >> econdidx >> captureblock >> captureblockidx >> slinkidx >> fedid >>
              zside;

          idx = modIndexer.getIndexForModule(fedid, captureblockidx, econdidx);

          typeidx = modIndexer.getTypeForModule(fedid, captureblockidx, econdidx);

          auto celltypes = modIndexer.convertTypeCode(typecode);
          isSiPM = celltypes.first;
          celltype = celltypes.second;

          eleid = HGCalElectronicsId((zside > 0), fedid, captureblockidx, econdidx, 0, 0).raw();

          if (!isSiPM) {
            int zp(zside > 0 ? 1 : -1);
            DetId::Detector det = plane <= 26 ? DetId::Detector::HGCalEE : DetId::Detector::HGCalHSi;
            detid = HGCSiliconDetId(det, zp, celltype, plane, i1, i2, 0, 0).rawId();
          }

          auto module = moduleParams.view()[idx];
          module.valid() = true;
          module.zside() = (zside > 0);
          module.isSiPM() = isSiPM;
          module.celltype() = celltype;
          module.plane() = plane;
          module.i1() = i1;
          module.i2() = i2;
          module.typeidx() = typeidx;
          module.fedid() = fedid;
          module.slinkidx() = slinkidx;
          module.captureblock() = captureblock;
          module.econdidx() = econdidx;
          module.captureblockidx() = captureblockidx;
          module.eleid() = eleid;
          module.detid() = detid;
        }

        return moduleParams;

      }  // end of produce()

    private:
      edm::ESGetToken<HGCalMappingModuleIndexer, HGCalMappingModuleIndexerRcd> moduleIndexTkn_;
      const edm::FileInPath filename_;
    };

  }  // namespace hgcal

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_EVENTSETUP_ALPAKA_MODULE(hgcal::HGCalMappingModuleESProducer);
