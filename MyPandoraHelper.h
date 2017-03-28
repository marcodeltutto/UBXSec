#ifndef MY_PANDORA_HELPER_H
#define MY_PANDORA_HELPER_H

#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

typedef std::map< art::Ptr<recob::PFParticle>, unsigned int > RecoParticleToNMatchedHits;
typedef std::map< art::Ptr<simb::MCParticle>,  RecoParticleToNMatchedHits > ParticleMatchingMap;
typedef std::set< art::Ptr<recob::PFParticle> > PFParticleSet;
typedef std::set< art::Ptr<simb::MCParticle> >  MCParticleSet;

class MyPandoraHelper {

  public:


  /**
   *  @brief Perform matching between true and reconstructed particles
   *
   *  @param recoParticlesToHits the mapping from reconstructed particles to hits
   *  @param trueHitsToParticles the mapping from hits to true particles
   *  @param matchedParticles the output matches between reconstructed and true particles
   *  @param matchedHits the output matches between reconstructed particles and hits
   */
  static void GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, const lar_pandora::HitsToMCParticles &trueHitsToParticles,
       lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits);
   /**
   *  @brief Perform matching between true and reconstructed particles
   *
   *  @param recoParticlesToHits the mapping from reconstructed particles to hits
   *  @param trueHitsToParticles the mapping from hits to true particles
   *  @param matchedParticles the output matches between reconstructed and true particles
   *  @param matchedHits the output matches between reconstructed particles and hits
   *  @param recoVeto the veto list for reconstructed particles
   *  @param trueVeto the veto list for true particles
   */
  static void GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, const lar_pandora::HitsToMCParticles &trueHitsToParticles,
       lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits, PFParticleSet &recoVeto, MCParticleSet &trueVeto, bool _recursiveMatching);

  void GetRecoToTrueMatches(art::Event const & e, std::string _pfp_producer, std::string _spacepointLabel, std::string _hitfinderLabel, std::string _geantModuleLabel, lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits);

  static void GetTPCObjects(art::Event const & e, std::string _particleLabel, std::vector<lar_pandora::PFParticleVector> & pfp_v_v, std::vector<lar_pandora::TrackVector> & track_v_v);

  static void GetTPCObjects(lar_pandora::PFParticleVector pfParticleList, lar_pandora::PFParticlesToTracks pfParticleToTrackMap, lar_pandora::PFParticlesToVertices  pfParticleToVertexMap, std::vector<lar_pandora::PFParticleVector> & pfp_v_v, std::vector<lar_pandora::TrackVector> & track_v_v);

  static void CollectTracksAndPFP(lar_pandora::PFParticlesToTracks pfParticleToTrackMap, lar_pandora::PFParticleVector pfParticleList, art::Ptr<recob::PFParticle> particle, lar_pandora::PFParticleVector &pfp_v, lar_pandora::TrackVector &track_v);

  static void GetNuVertexFromTPCObject(art::Event const & e, std::string _particleLabel, lar_pandora::PFParticleVector pfp_v, double *reco_nu_vtx);

  static art::Ptr<recob::PFParticle> GetNuPFP(lar_pandora::PFParticleVector pfp_v);

  static bool InFV(double * nu_vertex_xyz);

  static int GetSliceOrigin(std::vector<art::Ptr<recob::PFParticle>> neutrinoOriginPFP, lar_pandora::PFParticleVector pfp_v);

  static void GetNumberOfHitsPerPlane(art::Event const & e, std::string _particleLabel, lar_pandora::TrackVector track_v, int & nhits_u, int & nhits_v, int & nhits_w );

  static bool IsCrossingBoundary(recob::Track track, int & vtx_ok);

  static bool IsCrossingTopBoundary(recob::Track track, int & vtx_ok);

  static bool GetLongestTrackFromTPCObj(lar_pandora::TrackVector track_v, recob::Track & out_track);

  static bool PointIsCloseToDeadRegion(double *reco_nu_vtx, int plane_no);

};

#endif //  MY_PANDORA_HELPER_H
