#ifndef MY_PANDORA_HELPER_H
#define MY_PANDORA_HELPER_H

#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

typedef std::map< art::Ptr<recob::PFParticle>, unsigned int > RecoParticleToNMatchedHits;
typedef std::map< art::Ptr<simb::MCParticle>,  RecoParticleToNMatchedHits > ParticleMatchingMap;
typedef std::set< art::Ptr<recob::PFParticle> > PFParticleSet;
typedef std::set< art::Ptr<simb::MCParticle> >  MCParticleSet;

class UBXSecHelper {

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

  /**
   *  @brief Constructs TPC objects using Pandora PFP slices
   *
   *  @param e the ART event
   *  @param _particleLabel the PFP producer module
   *  @param pfp_v_v output, a vector of vector of PFP (a vector of TPC objects)
   *  @param track_v_v output, a vector of vector of tracks (a vector of TPC objects)
   */
  static void GetTPCObjects(art::Event const & e, std::string _particleLabel, std::vector<lar_pandora::PFParticleVector> & pfp_v_v, std::vector<lar_pandora::TrackVector> & track_v_v);

  /**
   *  @brief Constructs TPC objects using Pandora PFP slices
   *
   *  @param pfParticleList the list of PFP
   *  @param pfParticleToTrackMap map from PFP to tracks
   *  @param pfParticleToVertexMap map from PFP to vertices
   *  @param _particleLabel the PFP producer module
   *  @param pfp_v_v output, a vector of vector of PFP (a vector of TPC objects)
   *  @param track_v_v output, a vector of vector of tracks (a vector of TPC objects)   */
  static void GetTPCObjects(lar_pandora::PFParticleVector pfParticleList, lar_pandora::PFParticlesToTracks pfParticleToTrackMap, lar_pandora::PFParticlesToVertices  pfParticleToVertexMap, std::vector<lar_pandora::PFParticleVector> & pfp_v_v, std::vector<lar_pandora::TrackVector> & track_v_v);

  /**
   *  @brief Gets all the tracks and PFP for a single Pandora slice
   *
   *  @param pfParticleList the list of PFP
   *  @param pfParticleToTrackMap map from PFP to tracks
   *  @param particle the PFP
   *  @param pfp_v output, a vector of PFP (the TPC object)
   *  @param track_v output, a of tracks (the TPC object)   */
  static void CollectTracksAndPFP(lar_pandora::PFParticlesToTracks pfParticleToTrackMap, lar_pandora::PFParticleVector pfParticleList, art::Ptr<recob::PFParticle> particle, lar_pandora::PFParticleVector &pfp_v, lar_pandora::TrackVector &track_v);

  /**
   *  @brief Returns the nu recon vertex from a TPC object
   *
   *  @param e the ART event
   *  @param _particleLabel the PFP procuder module
   *  @param pfp_v the TPC object (vector of PFP)
   *  @param reco_nu output, the nu vertex  */
  static void GetNuVertexFromTPCObject(art::Event const & e, std::string _particleLabel, lar_pandora::PFParticleVector pfp_v, double *reco_nu_vtx);

  /**
   *  @brief Returns the nu PFP from a TPC object
   *
   *  @param pfp_v the TPC object (vector of PFP) */
  static art::Ptr<recob::PFParticle> GetNuPFP(lar_pandora::PFParticleVector pfp_v);

  /**
   *  @brief Returns true if the point passed is in the fiducial volume
   *
   *  @param nu_vertex_xyz a 3 dimensional array */
  static bool InFV(double * nu_vertex_xyz);

  /**
   *  @brief Returns the origin of the TPC object (neutrino/cosmic)
   *
   *  @param neutrinoOriginPFP list of PFP with neutrino origin
   *  @param pfp_v the TPC object (vector of PFP)  */
  static int GetSliceOrigin(std::vector<art::Ptr<recob::PFParticle>> neutrinoOriginPFP, lar_pandora::PFParticleVector pfp_v);

  /**
   *  @brief Returns the origin of the TPC object (neutrino/cosmic)
   *
   *  @param e the ART event
   *  @param _particleLabel the PFP procuder module
   *  @param track_v the TPC object (vector of tracks)
   *  @param nhits_u number of hits in the U plane
   *  @param nhits_v number of hits in the V plane
   *  @param nhits_w number of hits in the W plane  */
  static void GetNumberOfHitsPerPlane(art::Event const & e, std::string _particleLabel, lar_pandora::TrackVector track_v, int & nhits_u, int & nhits_v, int & nhits_w );

  /**
   *  @brief Returns true if the track is crossing the FV border
   *
   *  @param track the track
   *  @param vtx_ok is 0 if the vtx is in the FV, 1 otherwise  */
  static bool IsCrossingBoundary(recob::Track track, int & vtx_ok);

  /**
   *  @brief Returns true if the track is crossing the FV border from the TOP
   *
   *  @param track the track
   *  @param vtx_ok is 0 if the vtx is in the FV, 1 otherwise  */
  static bool IsCrossingTopBoundary(recob::Track track, int & vtx_ok);

  /**
   *  @brief Gets the longest track from a TPC object, returns false if there are no tracks in the TPC object
   *
   *  @param track_v the TPC object (vector of tracks)
   *  @param out_track the longest track  */
  static bool GetLongestTrackFromTPCObj(lar_pandora::TrackVector track_v, recob::Track & out_track);

  /**
   *  @brief Returns true if the point passed is close to a dead region
   *
   *  @param reco_nu_vtx the point to check
   *  @param plane_no the plane we want to check with  */
  static bool PointIsCloseToDeadRegion(double *reco_nu_vtx, int plane_no);

};

#endif //  MY_PANDORA_HELPER_H
