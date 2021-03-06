#include "Geometry/MuonNumbering/interface/GEMNumberingScheme.h"
#include "Geometry/MuonNumbering/interface/MuonBaseNumber.h"
#include "Geometry/MuonNumbering/interface/MuonDDDConstants.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//#define LOCAL_DEBUG

GEMNumberingScheme::GEMNumberingScheme(const MuonDDDConstants& muonConstants) { initMe(muonConstants); }

GEMNumberingScheme::GEMNumberingScheme(const DDCompactView& cpv) {
  MuonDDDConstants muonConstants(cpv);
  initMe(muonConstants);
}

void GEMNumberingScheme::initMe(const MuonDDDConstants& muonConstants) {
  int theLevelPart = muonConstants.getValue("level");
  theRegionLevel = muonConstants.getValue("mg_region") / theLevelPart;
  theStationLevel = muonConstants.getValue("mg_station") / theLevelPart;
  theRingLevel = muonConstants.getValue("mg_ring") / theLevelPart;
  theSectorLevel = muonConstants.getValue("mg_sector") / theLevelPart;
  theRollLevel = muonConstants.getValue("mg_roll") / theLevelPart;
#ifdef LOCAL_DEBUG
  edm::LogVerbatim("GEMNumberingScheme") << "Initialize GEMNumberingScheme"
                                         << "\ntheRegionLevel " << theRegionLevel << "\ntheStationLevel "
                                         << theStationLevel << "\ntheRingLevel " << theRingLevel << "\ntheSectorLevel "
                                         << theSectorLevel << "\ntheRollLevel " << theRollLevel;
#endif
}

int GEMNumberingScheme::baseNumberToUnitNumber(const MuonBaseNumber& num) {
#ifdef LOCAL_DEBUG
  edm::LogVerbatim("GEMNumberingScheme") << "GEMNumbering " << num.getLevels();
  for (int level = 1; level <= num.getLevels(); level++) {
    edm::LogVerbatim("GEMNumberingScheme") << level << " " << num.getSuperNo(level) << " " << num.getBaseNo(level);
  }
#endif

  int maxLevel = theRollLevel;
  if (num.getLevels() != maxLevel) {
    edm::LogWarning("GEMNumberingScheme")
        << "MuonGEMNumberingScheme::BNToUN: BaseNumber has " << num.getLevels() << " levels, need " << maxLevel;
    return 0;
  }

  int region(0), ring(0), station(0), layer(0), chamber(0), roll(0);

  //decode significant GEM levels

  if (num.getBaseNo(theRegionLevel) == 0)
    region = 1;
  else
    region = -1;

  // All GEM super chambers in stations 1 and 2 are on ring 1.
  // The long super chambers in station 2 are assigned *station 3* due
  // to the current limitation in the definition of the GEMDetId,
  // i.e. only 2 layers available per station.
  //  ring    = num.getSuperNo(theRingLevel);
  // GEM are only on the first ring
  ring = 1;
  station = num.getSuperNo(theStationLevel);
#ifdef LOCAL_DEBUG
  edm::LogVerbatim("GEMNumberingScheme") << "GEMNumbering: Ring " << ring << " Station "
                                         << num.getSuperNo(theStationLevel) << ":" << station;
#endif

  roll = num.getBaseNo(theRollLevel) + 1;
  const int copyno = num.getBaseNo(theSectorLevel) + 1;
  if (copyno < 50) {
    if (copyno % 2 == 0) {
      layer = 2;
      chamber = copyno - 1;
    } else {
      layer = 1;
      chamber = copyno;
    }
  } else {
    int copynp = copyno - 50;
    if (copynp % 2 != 0) {
      layer = 2;
      chamber = copynp - 1;
    } else {
      layer = 1;
      chamber = copynp;
    }
  }

  // collect all info

#ifdef LOCAL_DEBUG
  edm::LogVerbatim("GEMNumberingScheme") << "GEMNumberingScheme: Region " << region << " Ring " << ring << " Station "
                                         << station << " Layer " << layer << " Chamber " << chamber << " Roll " << roll;
#endif

  // Build the actual numbering
  GEMDetId id(region, ring, station, layer, chamber, roll);

#ifdef LOCAL_DEBUG
  edm::LogVerbatim("GEMNumberingScheme") << id.rawId() << " DetId " << id;
#endif

  return id.rawId();
}
