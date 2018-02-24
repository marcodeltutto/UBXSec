/**
 * \class UBXSecHelper
 *
 * \ingroup UBXSec
 *
 * \brief An helper class for UBXSec
 * 
 *
 * \author $Author: Marco Del Tutto<marco.deltutto@physics.ox.ac.uk> $
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2017/03/02 $
 *
 * Contact: marco.deltutto@physics.ox.ac.uk
 *
 * Created on: Friday, February 03, 2017 at 17:24:25
 *
 */


#ifndef UBXSECHELPER_H
#define UBXSECHELPER_H

#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

typedef std::map< art::Ptr<recob::PFParticle>, unsigned int > RecoParticleToNMatchedHits;
typedef std::map< art::Ptr<simb::MCParticle>,  RecoParticleToNMatchedHits > ParticleMatchingMap;
typedef std::set< art::Ptr<recob::PFParticle> > PFParticleSet;
typedef std::set< art::Ptr<simb::MCParticle> >  MCParticleSet;


class UBXSecHelper {

  public:

   /**
   *  @brief Given a vector of hits associated to the PFP (or a track, indeed you are just passing the hits), returns the tracking eff and purity
   *
   *  @param recoHits the input vector of reconstructed hits
   *  @param trackPurity the output track (PFP, whatever) purity 
   *  @param trackEfficiency the output track (PFP, whatever) efficiency
   */
  static void GetTrackPurityAndEfficiency( lar_pandora::HitVector recoHits, double & trackPurity, double & trackEfficiency );

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

  static void GetRecoToTrueMatches(art::Event const & e, std::string _pfp_producer, std::string _spacepointLabel, std::string _hitfinderLabel, std::string _geantModuleLabel, lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits);

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
   *  @param e the ART event
   *  @param _particleLabel the PFP producer module
   *  @param _trackLabel the track producer module
   *  @param pfp_v_v output, a vector of vector of PFP (a vector of TPC objects)
   *  @param track_v_v output, a vector of vector of tracks (a vector of TPC objects)
   */
  static void GetTPCObjects(art::Event const & e, std::string _particleLabel, std::string _trackLabel, std::vector<lar_pandora::PFParticleVector> & pfp_v_v, std::vector<lar_pandora::TrackVector> & track_v_v);

  /**
   *  @brief Constructs TPC objects using Pandora PFP slices
   *
   *  @param e the ART event
   *  @param _particleLabel the PFP producer module
   *  @param _trackLabel the track producer module
   *  @param pfp_v_v output, a vector of vector of PFP (a vector of TPC objects)
   *  @param track_v_v output, a vector of vector of tracks (a vector of TPC objects)
   *  @param pfp_to_spacept ouput, a map from pfp to space points
   *  @param spacept_to_hits ouput, a map from space points to recon hits
   */
  static void GetTPCObjects(art::Event const & e, std::string _particleLabel, std::string _trackLabel, std::vector<lar_pandora::PFParticleVector> & pfp_v_v, std::vector<lar_pandora::TrackVector> & track_v_v, lar_pandora::PFParticlesToSpacePoints & pfp_to_spacept, lar_pandora::SpacePointsToHits & spacept_to_hits);

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
   *  @param reco_nu_vtx output, the nu vertex (three dimensional array: x, y, z) */
  static void GetNuVertexFromTPCObject(art::Event const & e, std::string _particleLabel, lar_pandora::PFParticleVector pfp_v, double *reco_nu_vtx);

  /**
   *  @brief Returns the nu recon vertex from a TPC object
   *
   *  @param e the ART event
   *  @param _particleLabel the PFP procuder module
   *  @param pfp_v the TPC object (vector of PFP)
   *  @param reco_nu_vtx output, the nu vertex (recob::Vertex)  */
  static void GetNuVertexFromTPCObject(art::Event const & e, std::string _particleLabel, lar_pandora::PFParticleVector pfp_v, recob::Vertex & reco_nu_vtx);

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
   *  @brief Returns the origin of the TPC object (neutrino/cosmic/mixed)
   *
   *  @param neutrinoOriginPFP list of PFP with neutrino origin
   *  @param cosmicOriginPFP list of PFP with cosmic origin
   *  @param pfp_v the TPC object (vector of PFP)  */
  static ubana::TPCObjectOrigin GetSliceOrigin(std::vector<art::Ptr<recob::PFParticle>> neutrinoOriginPFP, std::vector<art::Ptr<recob::PFParticle>> cosmicOriginPFP, lar_pandora::PFParticleVector pfp_v);

  /**
   *  @brief Returns extra origin information for a TPC object (for now if it is a stopping muon in the TPC)
   *
   *  @param protonNCOriginPFP list of PFP with nc origin with proton
   *  @param pionNCOriginPFP list of PFP with nc origin with pion
   *  @param pfp_v the TPC object (vector of PFP)  */
  static ubana::TPCObjectOriginExtra GetSliceOriginExtra_NC(lar_pandora::PFParticleVector protonNCOriginPFP, lar_pandora::PFParticleVector pionNCOriginPFP, lar_pandora::PFParticleVector pfp_v);

  /**
   *  @brief Returns extra origin information for a TPC object (for now only if it is a stopping muon in the TPC)
   *
   *  @param cosmicStoppingOriginPFP list of PFP with cosmic origin that stop in the TPC
   *  @param pfp_v the TPC object (vector of PFP)  */
  static ubana::TPCObjectOriginExtra GetSliceOriginExtra_Stopping(lar_pandora::PFParticleVector cosmicStoppingOriginPFP, lar_pandora::PFParticleVector pfp_v);

  /**
   *  @brief Returns number of hits on each plane for a TPC obj
   *
   *  @param e the ART event
   *  @param _particleLabel the PFP procuder module
   *  @param track_v the TPC object (vector of tracks)
   *  @param nhits_u number of hits in the U plane
   *  @param nhits_v number of hits in the V plane
   *  @param nhits_w number of hits in the W plane  */
  static void GetNumberOfHitsPerPlane(art::Event const & e, std::string _particleLabel, lar_pandora::TrackVector track_v, int & nhits_u, int & nhits_v, int & nhits_w );

  /**
   *  @brief Returns true if the track passes a minimum hit requirment in at least one plane
   *
   *  @param e the ART event
   *  @param _particleLabel the PFP procuder module
   *  @param trk the recob::Track
   *  @param nHitsReq the minimum number of hits */
  static bool TrackPassesHitRequirment(art::Event const & e, std::string _particleLabel, art::Ptr<recob::Track> trk, int nHitsReq);

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
   *  @brief Gets the longest shower from a TPC object, returns false if there are no showers in the TPC object
   *
   *  @param shower_v the TPC object (vector of showers)
   *  @param out_shower the longest shower  */
  static bool GetLongestShowerFromTPCObj(lar_pandora::ShowerVector shower_v, recob::Shower & out_shower);

  /**
   *  @brief Returns true if all the tracks passed are fully contained 
   *
   *  @param tracks a vector of recob::Track  */
  static bool TracksAreContained(std::vector<recob::Track> tracks);

  /**
   *  @brief Returns true if the track passed is fully contained 
   *
   *  @param track the recob::Track */
  static bool TrackIsContained(recob::Track track);

  /**
   *  @brief Returns the value of the phi angle for the track considering the nu vertex
   *
   *  @param px  
   *  @param py 
   *  @param pz */
  static double GetPhi(double px, double py, double pz);

  /**
   *  @brief Returns the value of the phi angle for the track considering the nu vertex
   *
   *  @param t the recob::Track 
   *  @param tpcobj_nu_vtx the recob::Vertex neutrino vertex */
  static double GetCorrectedPhi(recob::Track t, recob::Vertex tpcobj_nu_vtx);

  /**
   *  @brief Returns the value of the cos(theta) angle for the track considering the nu vertex
   *
   *  @param track the recob::Track 
   *  @param tpcobj_nu_vtx the recob::Vertex neutrino vertex */
  static double GetCorrectedCosTheta(recob::Track t, recob::Vertex tpcobj_nu_vtx);

  /**
   *  @brief Returns the value of the phi angle
   *
   *  @param dir the direction */
  static double GetPhi(TVector3 dir);

  /**
   *  @brief Returns the value of the cos(theta) angle
   *
   *  @param dir the direction */
  static double GetCosTheta(TVector3 dir);

  /**
   *  @brief Returns true if the point passed is close to a dead region
   *
   *  @param reco_nu_vtx the point to check
   *  @param plane_no the plane we want to check with  */
  static bool PointIsCloseToDeadRegion(double *reco_nu_vtx, int plane_no);

  /**
   *  @brief Returns the clostes opchannel to the 3D specified point
   *
   *  @param charge_center the 3D point (c array) */
  static int GetClosestPMT(double *charge_center);


  /**
   *  @brief Returns Z position of a flash
   *
   *  @param pe a vector of PEs per OpDet */
  static double GetFlashZCenter(std::vector<double> pe); 

  /**
   *  @brief Takes a dimension 3 array, interaction time and drift velocity, and returns a corrected point (shifting in x)
   *
   *  @param point_raw raw 3d point 
   *  @param point_corrected corrected 3d point
   *  @param interaction_time The interaction time (usually from flash)
   *  @param drift_velocity The drift velocity */
  static void GetTimeCorrectedPoint(double * point_raw, double * point_corrected, double interaction_time, double drift_velocity);

  static std::vector<double> GetDqDxVector(std::vector<art::Ptr<anab::Calorimetry>> calos, int plane_no = 2);

  static double GetDqDxTruncatedMean(std::vector<art::Ptr<anab::Calorimetry>> calos, int plane_no = 2);

  static double GetDqDxTruncatedMean(std::vector<double> dqdx_v);

  static double GetMean(std::vector<double>);

  static double GetMedian(std::vector<double>);

  static double GetVariance(std::vector<double>);

  static double GetSTD(std::vector<double>);

  /**
   *  @brief Returns the MCTruth given a GEANT track ID
   *
   *  @param e the Event
   *  @param _geant_producer the gean producer label
   *  @param geant_track_id the geant track id */
  static art::Ptr<simb::MCTruth> TrackIDToMCTruth(art::Event const & e, std::string _geant_producer, int geant_track_id);
};


#endif //  UBXSECHELPER_H
