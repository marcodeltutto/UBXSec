////////////////////////////////////////////////////////////////////////
// Class:       TPCObjectMaker
// Plugin Type: producer (art v2_05_00)
// File:        TPCObjectMaker_module.cc
//
// Generated at Mon May 15 10:56:05 2017 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 * \file TPCObjectMaker_module.cc
 *
 * \ingroup UBXSec
 * 
 * \brief LArSoft plugin to generate TPC Objects
 *
 * @author Marco Del Tutto
 */

/** \addtogroup UBXSec

   @{*/
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

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larsim/MCCheater/BackTracker.h"

#include "uboone/UBXSec/DataTypes/TPCObject.h"
#include "uboone/UBXSec/Algorithms/McPfpMatch.h"
#include "uboone/UBXSec/Algorithms/UBXSecHelper.h"

#include <memory>

const simb::Origin_t NEUTRINO_ORIGIN = simb::kBeamNeutrino;
const simb::Origin_t COSMIC_ORIGIN   = simb::kCosmicRay;

namespace ubana {
  class TPCObjectMaker;
}


class ubana::TPCObjectMaker : public art::EDProducer {
public:
  explicit TPCObjectMaker(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TPCObjectMaker(TPCObjectMaker const &) = delete;
  TPCObjectMaker(TPCObjectMaker &&) = delete;
  TPCObjectMaker & operator = (TPCObjectMaker const &) = delete;
  TPCObjectMaker & operator = (TPCObjectMaker &&) = delete;

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
   *  @param _pfp_producer the PFP producer module
   *  @param pfp_v_v output, a vector of vector of PFP (a vector of TPC objects)
   *  @param track_v_v output, a vector of vector of tracks (a vector of TPC objects)   */
  void GetTPCObjects(lar_pandora::PFParticleVector pfParticleList, lar_pandora::PFParticlesToTracks pfParticleToTrackMap, lar_pandora::PFParticlesToVertices  pfParticleToVertexMap, std::vector<lar_pandora::PFParticleVector> & pfp_v_v, std::vector<lar_pandora::TrackVector> & track_v_v);


  /**
   *  @brief Returns the nu PFP from a TPC object
   *
   *  @param pfp_v the TPC object (vector of PFP) */
  art::Ptr<recob::PFParticle> GetNuPFP(lar_pandora::PFParticleVector pfp_v);

private:

  ubxsec::McPfpMatch mcpfpMatcher;
  bool _is_mc;

  std::string _pfp_producer;
  std::string _vertexLabel;
  std::string _trackLabel;
  std::string _hitfinderLabel;
  std::string _geantModuleLabel;
  std::string _spacepointLabel;
  bool _debug;
};


ubana::TPCObjectMaker::TPCObjectMaker(fhicl::ParameterSet const & p)
{
  _pfp_producer       = p.get<std::string>("PFParticleProducer"); 
  _vertexLabel        = p.get<std::string>("VertexProducer");
  _trackLabel         = p.get<std::string>("TrackProducer");
  _hitfinderLabel     = p.get<std::string>("HitProducer");
  _geantModuleLabel   = p.get<std::string>("GeantModule");
  _spacepointLabel    = p.get<std::string>("SpacePointProducer");
  _debug              = p.get<bool>       ("Debug", false);

  produces< std::vector<ubana::TPCObject>>();
  produces< art::Assns<ubana::TPCObject,   recob::Track>>();
  produces< art::Assns<ubana::TPCObject,   recob::PFParticle>>();
}

void ubana::TPCObjectMaker::produce(art::Event & e){

  if (_debug) std::cout << "[TPCObjectMaker] Starts" << std::endl;
 
  _is_mc = !e.isRealData();

  if (_is_mc) mcpfpMatcher.Configure(e, _pfp_producer, _spacepointLabel, _hitfinderLabel, _geantModuleLabel);

  // Instantiate the output
  std::unique_ptr< std::vector< ubana::TPCObject > >                tpcObjectVector        (new std::vector<ubana::TPCObject>);
  std::unique_ptr< art::Assns<ubana::TPCObject, recob::Track>>      assnOutTPCObjectTrack  (new art::Assns<ubana::TPCObject, recob::Track>      );
  std::unique_ptr< art::Assns<ubana::TPCObject, recob::PFParticle>> assnOutTPCObjectPFP    (new art::Assns<ubana::TPCObject, recob::PFParticle> );

  art::ServiceHandle<cheat::BackTracker> bt;

  // Vectors and maps we will use to store Pandora information
  lar_pandora::PFParticleVector pfParticleList;              //vector of PFParticles
  lar_pandora::PFParticlesToClusters pfParticleToClusterMap; //PFParticle-to-cluster map

  // Use LArPandoraHelper functions to collect Pandora information
  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, pfParticleList, pfParticleToClusterMap); //collect PFParticles and build map PFParticles to Clusters

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


  // Do the MCParticle to PFParticle matching

  lar_pandora::MCParticlesToPFParticles matchedParticles;    // This is a map: MCParticle to matched PFParticle
  lar_pandora::MCParticlesToHits        matchedParticleHits;
  if (_is_mc) {
    mcpfpMatcher.GetRecoToTrueMatches(matchedParticles, matchedParticleHits);
  }



  // Loop over true particle and find the pfp with cosmic and neutrino origin 

  std::vector<art::Ptr<recob::PFParticle>> neutrinoOriginPFP;
  std::vector<art::Ptr<recob::PFParticle>> cosmicOriginPFP;
  neutrinoOriginPFP.clear();
  cosmicOriginPFP.clear();

  for (lar_pandora::MCParticlesToPFParticles::const_iterator iter1 = matchedParticles.begin(), iterEnd1 = matchedParticles.end();
      iter1 != iterEnd1; ++iter1) {

    art::Ptr<simb::MCParticle>  mc_par = iter1->first;   // The MCParticle 
    art::Ptr<recob::PFParticle> pf_par = iter1->second;  // The matched PFParticle

    const art::Ptr<simb::MCTruth> mc_truth = bt->TrackIDToMCTruth(mc_par->TrackId());

    if (!mc_truth) {
      std::cerr << "[TPCObjectMaker] Problem with MCTruth pointer." << std::endl;
      continue;
    }

    if (mc_truth->Origin() == COSMIC_ORIGIN) {
      cosmicOriginPFP.emplace_back(pf_par);
    }

    if (mc_truth->Origin() == NEUTRINO_ORIGIN) {
      neutrinoOriginPFP.emplace_back(pf_par);
    }
  }


  // Construct TPCObjects

  for (size_t i = 0; i < pfp_v_v.size(); i++){

    ::ubana::TPCObject obj;

    // Set tracks
    std::vector<recob::Track> trk_v;
    trk_v.clear();
    for (auto t : track_v_v[i]) trk_v.emplace_back((*t));
    obj.SetTracks(trk_v);

    // Set PFPs
    std::vector<recob::PFParticle> pfp_v;
    pfp_v.clear();
    for (auto p : pfp_v_v[i]) pfp_v.emplace_back((*p));
    obj.SetPFPs(pfp_v);

    // Set vertex
    art::Ptr<recob::PFParticle> pfp = this->GetNuPFP(pfp_v_v[i]);
    auto iter = pfParticleToVertexMap.find(pfp);
    if (iter != pfParticleToVertexMap.end()) {
      obj.SetVertex(*(iter->second[0]));
    }

    // Set origin
    ::ubana::TPCObjectOrigin origin = ubana::kUnknown;
    if (_is_mc)
      origin = UBXSecHelper::GetSliceOrigin(neutrinoOriginPFP, cosmicOriginPFP, pfp_v_v[i]); 
    obj.SetOrigin(origin);

    tpcObjectVector->emplace_back(obj);
    util::CreateAssn(*this, e, *tpcObjectVector, track_v_v[i],  *assnOutTPCObjectTrack);
    util::CreateAssn(*this, e, *tpcObjectVector, pfp_v_v[i],    *assnOutTPCObjectPFP);
  }



  // Put TPCObjects into the Event

  e.put(std::move(tpcObjectVector)); 
  e.put(std::move(assnOutTPCObjectTrack));
  e.put(std::move(assnOutTPCObjectPFP));


  if (_debug) std::cout << "[TPCObjectMaker] Ends" << std::endl;
}



//_____________________________________________________________________________________
art::Ptr<recob::PFParticle> ubana::TPCObjectMaker::GetNuPFP(lar_pandora::PFParticleVector pfp_v){

  for (unsigned int pfp = 0; pfp < pfp_v.size(); pfp++) {

    if(lar_pandora::LArPandoraHelper::IsNeutrino(pfp_v.at(pfp))) {
      return pfp_v.at(pfp);
    }
  }
  std::cout << "[TPCObjectMaker] No neutrino PFP found." << std::endl;

  art::Ptr<recob::PFParticle> temp;
  return temp;

}



//___________________________________________________________________________________________________
void ubana::TPCObjectMaker::GetTPCObjects(lar_pandora::PFParticleVector pfParticleList,
                                            lar_pandora::PFParticlesToTracks pfParticleToTrackMap,
                                            lar_pandora::PFParticlesToVertices  pfParticleToVertexMap,
                                            std::vector<lar_pandora::PFParticleVector> & pfp_v_v,
                                            std::vector<lar_pandora::TrackVector> & track_v_v) {

  track_v_v.clear();
  pfp_v_v.clear();

  if (_debug) std::cout << "[TPCObjectMaker] Getting TPC Objects..." << std::endl;

  for (unsigned int n = 0; n < pfParticleList.size(); ++n) {
    const art::Ptr<recob::PFParticle> particle = pfParticleList.at(n);

    if(lar_pandora::LArPandoraHelper::IsNeutrino(particle)) {
      if (_debug) std::cout << "[TPCObjectMaker] \t Creating TPC Object " << track_v_v.size() << std::endl;
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

      if (_debug) std::cout << "[TPCObjectMaker] \t Number of pfp for this TPC object: "    << pfp_v.size()   << std::endl;
      for (auto pfp : pfp_v) {
        if (_debug) std::cout << "[TPCObjectMaker] \t \t PFP " << pfp->Self() << " with pdg " << pfp->PdgCode();
        auto it = pfParticleToVertexMap.find(pfp);
        if (it == pfParticleToVertexMap.end()) {
           if (_debug) std::cout << " and vertex [vertex not available for this PFP]" << std::endl;
        } else {
          double xyz[3];
          (it->second)[0]->XYZ(xyz);
          if (_debug) std::cout << " and vertex " << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
        }
      }
      std::cout << std::endl;
      if (_debug) std::cout << "[TPCObjectMaker] \t Number of tracks for this TPC object: " << track_v.size() << std::endl;

      //for (unsigned int i = 0; i < pfp_v.size(); i++) std::cout << "   pfp with ID " << pfp_v[i]->Self() << std::endl;
    } // end if neutrino
  } // end pfp loop
}

//______________________________________________________________________________________________________________________________________
void ubana::TPCObjectMaker::CollectTracksAndPFP(lar_pandora::PFParticlesToTracks pfParticleToTrackMap,
                                                  lar_pandora::PFParticleVector pfParticleList,
                                                  art::Ptr<recob::PFParticle> particle,
                                                  lar_pandora::PFParticleVector &pfp_v,
                                                  lar_pandora::TrackVector &track_v) {

  pfp_v.emplace_back(particle);

  lar_pandora::PFParticlesToTracks::const_iterator trackMapIter = pfParticleToTrackMap.find(particle);
  if (trackMapIter != pfParticleToTrackMap.end()) {
    lar_pandora::TrackVector tracks = trackMapIter->second;
    if (_debug) std::cout << "[TPCObjectMaker] \t PFP " << particle->Self() << " has " << tracks.size() << " tracks ass." << std::endl;
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


DEFINE_ART_MODULE(ubana::TPCObjectMaker)

  /** @} */ // end of doxygen group
