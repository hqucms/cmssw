#include <iostream>
#include <memory>
#include <cmath>
#include <numeric>

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" //for pat
#include "DataFormats/Candidate/interface/Candidate.h" // for reco

#include "PhysicsTools/NearestNeighbors/interface/NearestNeighbors.h"

#include "ABCNetMakeInputs.h"

using namespace abcnet;
using PointCloud = cms::nanoflann::PointCloud<float, 2>; //default value for the DIM parameter is 3, but want to gather using 2 dims only (eta, phi). See template definition L11 of PhysicsTools/NearestNeighbors/interface/NearestNeighbors.h

std::tuple< std::unordered_map<std::string, std::vector<float>>, std::vector<float> > ABCNetMakeInputs::makeFeatureMap (const reco::CandidateView * pfCol, const std::vector<size_t> & indices, bool debug) {
  
  //store PF candidates and their pts into vectors
  std::vector<const pat::PackedCandidate*> PFCands;
  
  for (const auto &i : indices) {
    const pat::PackedCandidate *lPack = dynamic_cast<const pat::PackedCandidate*>(&pfCol->at(i));
    if(lPack == nullptr) { 
      // throw error
      throw edm::Exception(edm::errors::LogicError,"ABCNetProducer: cannot get weights since inputs are not PackedCandidates");
    }
    PFCands.push_back(lPack);
  } // end loop over PF candidates

  //fill feature map
  std::unordered_map<std::string, std::vector<float>> fts;
  PointCloud::Points points(4000, {{0.f, 0.f}});

  for (unsigned i = 0; i < PFCands.size() && i < 4000; ++i) {
    const auto &aPF = *PFCands[i];

    points[i] = {{float(aPF.eta()), float(aPF.phi())}};
    
    fts["PFCandEta"].push_back(aPF.eta()); //f0
    fts["PFCandPhi"].push_back(aPF.phi()); //f1
    fts["PFCandLogPt"].push_back(log(aPF.pt())); //f2
    fts["PFCandLogE"].push_back(log(aPF.energy())); //f3c
    fts["PFCandCharge"].push_back(aPF.charge()); //f4
    if (abs(aPF.pdgId()) == 11) fts["PFCandIsEle"].push_back(1.0); else fts["PFCandIsEle"].push_back(0.0); //f5
    if (abs(aPF.pdgId()) == 13) fts["PFCandIsMu"].push_back(1.0); else fts["PFCandIsMu"].push_back(0.0); //f6
    if (abs(aPF.pdgId()) == 211 || abs(aPF.pdgId()) == 2 || abs(aPF.pdgId()) == 1) fts["PFCandIsHad"].push_back(1.0); else fts["PFCandIsHad"].push_back(0.0); //f7
    if (abs(aPF.pdgId()) == 130) fts["PFCandIsNeutralHad"].push_back(1.0); else fts["PFCandIsNeutralHad"].push_back(0.0); //f8
    if (abs(aPF.pdgId()) == 22) fts["PFCandIsPhoton"].push_back(1.0); else fts["PFCandIsPhoton"].push_back(0.0); //f9
    fts["PFCandNumHits"].push_back(aPF.numberOfHits()); //f14
    fts["PFCandNumLayersHit"].push_back(aPF.trackerLayersWithMeasurement()); //f15
    fts["PFCandFromPV"].push_back(aPF.fromPV()); //f16
    fts["PFCandTrackHighPurity"].push_back(aPF.trackHighPurity()); //f17
    const auto *trk = aPF.bestTrack();
    if (trk) {
      if ( isinf(aPF.dz()/aPF.dzError()) ) fts["PFCandDZSig"].push_back(999.0); else fts["PFCandDZSig"].push_back(aPF.dz()/aPF.dzError()); //f10
      if ( isinf(aPF.dxy()/aPF.dxyError()) ) fts["PFCandDXYSig"].push_back(999.0); else fts["PFCandDXYSig"].push_back(aPF.dxy()/aPF.dxyError()); //f11
      fts["PFCandNormChi2"].push_back(trk->normalizedChi2()); //f13
      fts["PFCandQuality"].push_back(trk->qualityMask()); //f18
    } else {
      fts["PFCandDZSig"].push_back(0.0);
      fts["PFCandDXYSig"].push_back(0.0);
      fts["PFCandNormChi2"].push_back(999.0);
      fts["PFCandQuality"].push_back(0.0);
    }
    if (aPF.pdgId() == 130 || aPF.pdgId() == 1) fts["PFCandHCalFrac"].push_back(aPF.hcalFraction()); else if (aPF.isIsolatedChargedHadron()) fts["PFCandHCalFrac"].push_back(aPF.rawHcalFraction()); else fts["PFCandHCalFrac"].push_back(0.0); //f12
  }
  if (debug) {
    std::cout << "*** NOW CHECKING THE SORTING OF PF CANDIDATES ***" << std::endl;
    if ( !std::is_sorted(fts["PFCandLogPt"].begin(), fts["PFCandLogPt"].end(), std::greater<float>() ) ) std::cout << "*** WARNING *** PFCandidates are not sorted in Pt!!!" << std::endl;
    else std::cout << "No issues with the sorting of PF candidates detected" << std::endl;
  }

  auto knn_indices = PointCloud::knn<float>(points, points, 21);  // queries.size() * num_neighbors
  std::vector<float> kNNs;
  kNNs.reserve(4000 * 20);
  for (unsigned i = 0; i < 4000; ++i){
    for (unsigned j = 1; j < 21; ++j){
      kNNs.push_back(knn_indices.at(i * 21 + j));
    }
  }

  return {fts, kNNs};
};
