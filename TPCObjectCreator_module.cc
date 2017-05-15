////////////////////////////////////////////////////////////////////////
// Class:       TPCObjectCreator
// Plugin Type: producer (art v2_05_00)
// File:        TPCObjectCreator_module.cc
//
// Generated at Mon May 15 10:56:05 2017 by Marco Del Tutto using cetskelgen
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

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <memory>

namespace ubana {
  class TPCObjectCreator;
}


class ubana::TPCObjectCreator : public art::EDProducer {
public:
  explicit TPCObjectCreator(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TPCObjectCreator(TPCObjectCreator const &) = delete;
  TPCObjectCreator(TPCObjectCreator &&) = delete;
  TPCObjectCreator & operator = (TPCObjectCreator const &) = delete;
  TPCObjectCreator & operator = (TPCObjectCreator &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  /**
   *  @brief Gets all the tracks and PFP for a single Pandora slice
   *
   *  @param pfParticleList the list of PFP
   *  @param pfParticleToTrackMap map from PFP to tracks
   *  @param particle the PFP
   *  @param pfp_v output, a vector of PFP (the TPC object)
   *  @param track_v output, a of tracks (the TPC object)   */
  void CollectTracksAndPFP(lar_pandora::PFParticlesToTracks pfParticleToTrackMap, lar_pandora::PFParticleVector pfParticleList, art::Ptr<recob::PFParticle> particle, lar_pandora::PFParticleVector &pfp_v, lar_pandora::TrackVector &track_v);

  /**
   *  @brief Constructs TPC objects using Pandora PFP slices
   *
   *  @param pfParticleList the list of PFP
   *  @param pfParticleToTrackMap map from PFP to tracks
   *  @param pfParticleToVertexMap map from PFP to vertices
   *  @param _particleLabel the PFP producer module
   *  @param pfp_v_v output, a vector of vector of PFP (a vector of TPC objects)
   *  @param track_v_v output, a vector of vector of tracks (a vector of TPC objects)   */
  void GetTPCObjects(lar_pandora::PFParticleVector pfParticleList, lar_pandora::PFParticlesToTracks pfParticleToTrackMap, lar_pandora::PFParticlesToVertices  pfParticleToVertexMap, std::vector<lar_pandora::PFParticleVector> & pfp_v_v, std::vector<lar_pandora::TrackVector> & track_v_v);

private:

  std::string _particleLabel;
  std::string _vertexLabel;
  std::string _trackLabel;
};


ubana::TPCObjectCreator::TPCObjectCreator(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  _particleLabel = p.get<std::string>("PFParticleProducer"); 
  _vertexLabel   = p.get<std::string>("VertexProducer");
  _trackLabel    = p.get<std::string>("TrackProducer");
  // Call appropriate produces<>() functions here.
}

void ubana::TPCObjectCreator::produce(art::Event & e){

  //Vectors and maps we will use to store Pandora information
  lar_pandora::PFParticleVector pfParticleList;              //vector of PFParticles
  lar_pandora::PFParticlesToClusters pfParticleToClusterMap; //PFParticle-to-cluster map

  //Use LArPandoraHelper functions to collect Pandora information
  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _particleLabel, pfParticleList, pfParticleToClusterMap); //collect PFParticles and build map PFParticles to Clusters

  // Collect vertices and tracks
  lar_pandora::VertexVector           allPfParticleVertices;
  lar_pandora::PFParticlesToVertices  pfParticleToVertexMap;
  lar_pandora::LArPandoraHelper::CollectVertices(e, _vertexLabel, allPfParticleVertices, pfParticleToVertexMap);
  lar_pandora::TrackVector            allPfParticleTracks;
  lar_pandora::PFParticlesToTracks    pfParticleToTrackMap;
  lar_pandora::LArPandoraHelper::CollectTracks(e, _trackLabel, allPfParticleTracks, pfParticleToTrackMap);

  std::vector<lar_pandora::TrackVector     > track_v_v;
  std::vector<lar_pandora::PFParticleVector> pfp_v_v;

  this->GetTPCObjects(pfParticleList, pfParticleToTrackMap, pfParticleToVertexMap, pfp_v_v, track_v_v);


}

void ubana::TPCObjectCreator::GetTPCObjects(lar_pandora::PFParticleVector pfParticleList,
                                            lar_pandora::PFParticlesToTracks pfParticleToTrackMap,
                                            lar_pandora::PFParticlesToVertices  pfParticleToVertexMap,
                                            std::vector<lar_pandora::PFParticleVector> & pfp_v_v,
                                            std::vector<lar_pandora::TrackVector> & track_v_v) {

  track_v_v.clear();
  pfp_v_v.clear();

  std::cout << "[UBXSecHelper] Getting TPC Objects..." << std::endl;

  for (unsigned int n = 0; n < pfParticleList.size(); ++n) {
    const art::Ptr<recob::PFParticle> particle = pfParticleList.at(n);

    if(lar_pandora::LArPandoraHelper::IsNeutrino(particle)) {
      std::cout << "[UBXSecHelper] \t Creating TPC Object " << track_v_v.size() << std::endl;
      //std::cout << "IS NEUTRINO, pfp id " << particle->Self() << std::endl;
      lar_pandora::VertexVector nu_vertex_v;
      auto search = pfParticleToVertexMap.find(particle);
      if(search != pfParticleToVertexMap.end()) {
        nu_vertex_v = search->second;
      }

      double nu_vertex_xyz[3]={0.,0.,0.};
      nu_vertex_v[0]->XYZ(nu_vertex_xyz);
      //if (!this->InFV(nu_vertex_xyz)) continue;

      lar_pandora::TrackVector track_v;
      lar_pandora::PFParticleVector pfp_v;

      this->CollectTracksAndPFP(pfParticleToTrackMap, pfParticleList, particle, pfp_v, track_v);

      pfp_v_v.emplace_back(pfp_v);
      track_v_v.emplace_back(track_v);

      std::cout << "[UBXSecHelper] \t Number of pfp for this TPC object: "    << pfp_v.size()   << std::endl;
      for (auto pfp : pfp_v) {
        std::cout << "[UBXSecHelper] \t \t PFP " << pfp->Self() << " with pdg " << pfp->PdgCode();
        auto it = pfParticleToVertexMap.find(pfp);
        if (it == pfParticleToVertexMap.end()) {
           std::cout << " and vertex [vertex not available for this PFP]" << std::endl;
        } else {
          double xyz[3];
          (it->second)[0]->XYZ(xyz);
          std::cout << " and vertex " << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
        }
      }
      std::cout << std::endl;
      std::cout << "[UBXSecHelper] \t Number of tracks for this TPC object: " << track_v.size() << std::endl;

      //for (unsigned int i = 0; i < pfp_v.size(); i++) std::cout << "   pfp with ID " << pfp_v[i]->Self() << std::endl;
    } // end if neutrino
  } // end pfp loop
}

//______________________________________________________________________________________________________________________________________
void ubana::TPCObjectCreator::CollectTracksAndPFP(lar_pandora::PFParticlesToTracks pfParticleToTrackMap,
                                                  lar_pandora::PFParticleVector pfParticleList,
                                                  art::Ptr<recob::PFParticle> particle,
                                                  lar_pandora::PFParticleVector &pfp_v,
                                                  lar_pandora::TrackVector &track_v) {

  pfp_v.emplace_back(particle);

  lar_pandora::PFParticlesToTracks::const_iterator trackMapIter = pfParticleToTrackMap.find(particle);
  if (trackMapIter != pfParticleToTrackMap.end()) {
    lar_pandora::TrackVector tracks = trackMapIter->second;
    std::cout << "[UBXSecHelper] \t PFP " << particle->Self() << " has " << tracks.size() << " tracks ass." << std::endl;
    for (unsigned int trk = 0; trk < tracks.size(); trk++) {
      track_v.emplace_back(tracks[trk]);
    }
  }
  const std::vector<size_t> &daughterIDs = particle->Daughters();
  if(daughterIDs.size() == 0) return;
  else {
    for (unsigned int m = 0; m < daughterIDs.size(); ++m) {
      const art::Ptr<recob::PFParticle> daughter = pfParticleList.at(daughterIDs.at(m));
      CollectTracksAndPFP(pfParticleToTrackMap, pfParticleList, daughter, pfp_v, track_v);
    }
  }

}

DEFINE_ART_MODULE(ubana::TPCObjectCreator)
