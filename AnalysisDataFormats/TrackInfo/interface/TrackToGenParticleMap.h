#ifndef TrackToGenParticleMap_h
#define TrackToGenParticleMap_h

#include <vector>
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

namespace reco {
    // map
    typedef std::map<edm::Ptr<reco::Candidate>, edm::Ptr<pat::PackedGenParticle>> TrackToGenParticleMap;
}

#endif