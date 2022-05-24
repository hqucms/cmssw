#include <iostream>
#include <memory>
#include <cmath>
#include <numeric>

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" //for pat
#include "DataFormats/Candidate/interface/Candidate.h" // for reco
#include "ABCNetMakeInputs.h"

using namespace abcnet;

std::unordered_map<std::string, std::vector<float>> ABCNetMakeInputs::makeFeatureMap (const reco::CandidateView * pfCol, bool debug) {
  std::unordered_map<std::string, std::vector<float>> fts;
  for (auto const & aPF : *pfCol ) {

    const pat::PackedCandidate *lPack = dynamic_cast<const pat::PackedCandidate*>(&aPF);
    fts["PFCandEta"].push_back(lPack->eta()); //f0
    fts["PFCandPhi"].push_back(lPack->phi()); //f1
    fts["PFCandLogPt"].push_back(log(lPack->pt())); //f2
    fts["PFCandLogE"].push_back(log(lPack->energy())); //f3
    fts["PFCandCharge"].push_back(lPack->charge()); //f4
    if (abs(lPack->pdgId()) == 11) fts["PFCandIsEle"].push_back(1.0); else fts["PFCandIsEle"].push_back(0.0); //f5
    if (abs(lPack->pdgId()) == 13) fts["PFCandIsMu"].push_back(1.0); else fts["PFCandIsMu"].push_back(0.0); //f6
    if (abs(lPack->pdgId()) == 211 || abs(lPack->pdgId()) == 2 || abs(lPack->pdgId()) == 1) fts["PFCandIsHad"].push_back(1.0); else fts["PFCandIsHad"].push_back(0.0); //f7
    if (abs(lPack->pdgId()) == 130) fts["PFCandIsK"].push_back(1.0); else fts["PFCandIsK"].push_back(0.0); //f8
    if (abs(lPack->pdgId()) == 22) fts["PFCandIsPhoton"].push_back(1.0); else fts["PFCandIsPhoton"].push_back(0.0); //f9
    fts["PFCandNumHits"].push_back(lPack->numberOfHits()); //f14
    fts["PFCandNumLayersHit"].push_back(lPack->trackerLayersWithMeasurement()); //f15
    fts["PFCandFromPV"].push_back(lPack->fromPV()); //f16
    fts["PFCandTrackHighPurity"].push_back(lPack->trackHighPurity()); //f17
    if (lPack->bestTrack()) {
      if ( isinf(lPack->dz()/lPack->dzError()) ) fts["PFCandDZSig"].push_back(999.0); else fts["PFCandDZSig"].push_back(lPack->dz()/lPack->dzError()); //f10
      if ( isinf(lPack->dxy()/lPack->dxyError()) ) fts["PFCandDXYSig"].push_back(999.0); else fts["PFCandDXYSig"].push_back(lPack->dxy()/lPack->dxyError()); //f11
      const auto *trk = lPack->bestTrack();
      fts["PFCandNormChi2"].push_back(trk->normalizedChi2()); //f13
      fts["PFCandQuality"].push_back(trk->qualityMask()); //f18
    }
    else {
      fts["PFCandDZSig"].push_back(0.0);
      fts["PFCandDXYSig"].push_back(0.0);
      fts["PFCandNormChi2"].push_back(999.0);
      fts["PFCandQuality"].push_back(0.0);
    }
    if (lPack->pdgId() == 130 || lPack->pdgId() == 1) fts["PFCandHCalFrac"].push_back(lPack->hcalFraction()); else if (lPack->isIsolatedChargedHadron()) fts["PFCandHCalFrac"].push_back(lPack->rawHcalFraction()); else fts["PFCandHCalFrac"].push_back(0.0); //f12
  }
  
  return fts;

};
