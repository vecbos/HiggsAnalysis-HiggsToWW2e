// -*- C++ -*-
#ifndef ESHELPER_H
#define ESHELPER_H

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"

#include <vector>
#include <map>

class ESHelper {
public:
  ESHelper(std::map<DetId, EcalRecHit> rechits_map, const CaloSubdetectorGeometry* geometry, const CaloSubdetectorTopology* topology);
  ~ESHelper() {};

  std::vector<float> getESHits(double x, double y, double z, int row);
  std::vector<float> getESShape(std::vector<float> ESHits0);
  std::vector<float> getESEn(std::vector<float> ESHits0);

private:
  std::map<DetId, EcalRecHit> _rechits_map;
  const CaloSubdetectorGeometry *_geometry;
  const CaloSubdetectorTopology *_topology;
};

#endif

