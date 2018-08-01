////////////////////////////////////////////////////////////////////////
// Class:       CosmicTrackHitTagger
// Plugin Type: producer (art v2_05_00)
// File:        CosmicTrackHitTagger_module.cc
//
// Generated at Fri Aug 18 16:54:54 2017 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"


#include "uboone/UBXSec/Algorithms/CosmicTagToolInterface.h"
#include "uboone/UBXSec/Algorithms/CosmicTagByHitIntegral.h"

#include <memory>

/*
namespace ubana {
  struct SimpleHit {
    double time;
    double wire;
    double integral;
  };
}

using SimpleHitVector = std::vector<ubana::SimpleHit>;
*/

class CosmicTrackHitTagger;


class CosmicTrackHitTagger : public art::EDProducer {
public:
  explicit CosmicTrackHitTagger(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CosmicTrackHitTagger(CosmicTrackHitTagger const &) = delete;
  CosmicTrackHitTagger(CosmicTrackHitTagger &&) = delete;
  CosmicTrackHitTagger & operator = (CosmicTrackHitTagger const &) = delete;
  CosmicTrackHitTagger & operator = (CosmicTrackHitTagger &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  /// Takes a vector of hits and returns a vector of hits only in the collection plane
  std::vector<art::Ptr<recob::Hit>> FilterHits(std::vector<art::Ptr<recob::Hit>>);

  std::string _pfp_producer;
  std::string _track_producer;
  std::string _hit_producer;

  bool _debug;

  // Insitantiate algorithm
  ::ubana::CosmicTagByHitIntegral _algo;
};


CosmicTrackHitTagger::CosmicTrackHitTagger(fhicl::ParameterSet const & p) {

  _pfp_producer       = p.get<std::string>("PFPartProducer", "pandoraCosmic");
  _track_producer     = p.get<std::string>("TrackProducer", "pandoraCosmic");
  _hit_producer       = p.get<std::string>("HitProducer", "swtrigger");
  _debug              = p.get<bool>("DebugMode", "false");
}

void CosmicTrackHitTagger::produce(art::Event & e) {

   // load PFParticles
  if (_debug) { std::cout << "Loading PFParticles from producer " << _pfp_producer << std::endl; }
  art::Handle<std::vector<recob::PFParticle> > pfp_h;
  e.getByLabel(_pfp_producer,pfp_h);

  // make sure pfparticles look good
  if(!pfp_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate PFParticles!"<<std::endl;
 /*   e.put(std::move(cosmicTagTrackVector));
    e.put(std::move(assnOutCosmicTagTrack));
    e.put(std::move(assnOutCosmicTagPFParticle));
*/
    return; //throw std::exception();
  }

  // load Tracks
  if (_debug) { std::cout << "Loading Tracks from producer " << _track_producer << std::endl; }
  art::Handle<std::vector<recob::Track> > track_h;
  e.getByLabel(_track_producer, track_h);

  // make sure pfparticles look good
  if(!track_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Tracks!"<<std::endl;
  /*  e.put(std::move(cosmicTagTrackVector));
    e.put(std::move(assnOutCosmicTagTrack));
    e.put(std::move(assnOutCosmicTagPFParticle));
*/
    return; //throw std::exception();
  }

  // grab tracks associated with PFParticles
  art::FindManyP<recob::Track> pfp_track_assn_v(pfp_h, e, _track_producer);

  // grab hits associated to Tracks
  art::FindManyP<recob::Hit> hits_from_track(track_h, e, _track_producer);

  std::vector<art::Ptr<recob::PFParticle> > PFPVec;
  art::fill_ptr_vector(PFPVec, pfp_h);



  // PFParticle loop
  for (size_t i = 0; i < PFPVec.size(); i++) {

    bool isCosmic = false;
    auto pfp = PFPVec.at(i);

    // grab associated tracks
    std::vector<art::Ptr<recob::Track>> track_v = pfp_track_assn_v.at(i);

    if (_debug) {
      std::cout << "[CosmicTrackHitTagger] Looping through pfpart number " << i << std::endl;
      //std::cout << "PFPart has " << track_v.size() << " tracks associated" << std::endl;
    }

    // Track loop
    for (auto track : track_v) {

      if (_debug) {
        std::cout << "[CosmicTrackHitTagger] Track start " << track->Vertex().X() << ", " << track->Vertex().Y() << ", " << track->Vertex().Z() << std::endl;
        std::cout << "[CosmicTrackHitTagger] Track end   " << track->End().X() << ", " << track->End().Y() << ", " << track->End().Z() << std::endl;
      }

      std::vector<art::Ptr<recob::Hit>> hits = hits_from_track.at(track.key());
      hits = this->FilterHits(hits);
      
      ubana::SimpleHitVector simple_hit_v;

      for (auto hit : hits) {

        ubana::SimpleHit_t simple_hit;
   
        simple_hit.time     = hit->PeakTime();
        simple_hit.wire     = hit->WireID().Wire;
        simple_hit.integral = hit->Integral();

        simple_hit_v.emplace_back(simple_hit);
      }

      isCosmic = _algo.IsCosmic(simple_hit_v);
    }

    std::cout << "[CosmicTrackHitTagger] PFP " << pfp->Self() << " " << (isCosmic ? "is" : "is not") << " a cosmic" << std::endl;
  }
}

std::vector<art::Ptr<recob::Hit>> CosmicTrackHitTagger::FilterHits(std::vector<art::Ptr<recob::Hit>> hits){

  std::vector<art::Ptr<recob::Hit>> output;

  for (auto hit : hits) {

    if (hit->View() != 2) continue;

    output.emplace_back(hit);
  }

  return output;

}

DEFINE_ART_MODULE(CosmicTrackHitTagger)
