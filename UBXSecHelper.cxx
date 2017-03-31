#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"


#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Wire.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "lardataobj/RecoBase/PFParticle.h"
//#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "uboone/UBXSec/UBXSecHelper.h"



//___________________________________________________________________________________________________
void UBXSecHelper::GetRecoToTrueMatches(art::Event const & e,
                                           std::string _pfp_producer,
                                           std::string _spacepointLabel,
                                           std::string _geantModuleLabel,
                                           std::string _hitfinderLabel,
                                           lar_pandora::MCParticlesToPFParticles &matchedParticles,
                                           lar_pandora::MCParticlesToHits &matchedParticleHits)
{

   bool _debug = true; 

  // --- Collect hits
  lar_pandora::HitVector hitVector;
  lar_pandora::LArPandoraHelper::CollectHits(e, _hitfinderLabel, hitVector);
 
   // --- Collect PFParticles and match Reco Particles to Hits
  lar_pandora::PFParticleVector  recoParticleVector;
  lar_pandora::PFParticleVector  recoNeutrinoVector;
  lar_pandora::PFParticlesToHits recoParticlesToHits;
  lar_pandora::HitsToPFParticles recoHitsToParticles;

  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, recoParticleVector);
  lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(recoParticleVector, recoNeutrinoVector);
  lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e, _pfp_producer, _spacepointLabel, recoParticlesToHits, recoHitsToParticles, lar_pandora::LArPandoraHelper::kAddDaughters);

  if (_debug) {
    std::cout << "  RecoNeutrinos: " << recoNeutrinoVector.size() << std::endl;
    std::cout << "  RecoParticles: " << recoParticleVector.size() << std::endl;
  }

  // --- Collect MCParticles and match True Particles to Hits
  lar_pandora::MCParticleVector     trueParticleVector;
  lar_pandora::MCTruthToMCParticles truthToParticles;
  lar_pandora::MCParticlesToMCTruth particlesToTruth;
  lar_pandora::MCParticlesToHits    trueParticlesToHits;
  lar_pandora::HitsToMCParticles    trueHitsToParticles;

  if (!e.isRealData()) {
    lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geantModuleLabel, trueParticleVector);
    lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geantModuleLabel, truthToParticles, particlesToTruth);
    lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(e, _geantModuleLabel, hitVector, trueParticlesToHits, trueHitsToParticles, lar_pandora::LArPandoraHelper::kAddDaughters);
  }

  if (_debug) {
    std::cout << "  TrueParticles: " << particlesToTruth.size() << std::endl;
    std::cout << "  TrueEvents: " << truthToParticles.size() << std::endl;
  }


  GetRecoToTrueMatches(recoParticlesToHits,
                       trueHitsToParticles,
                       matchedParticles,
                       matchedParticleHits);
}


//___________________________________________________________________________________________________
void UBXSecHelper::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits,
                                           const lar_pandora::HitsToMCParticles &trueHitsToParticles,
                                           lar_pandora::MCParticlesToPFParticles &matchedParticles,
                                           lar_pandora::MCParticlesToHits &matchedHits) 
{
    PFParticleSet recoVeto; MCParticleSet trueVeto;
    bool _recursiveMatching = false;

    GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedHits, recoVeto, trueVeto, _recursiveMatching);
}

//___________________________________________________________________________________________________
void UBXSecHelper::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits,
                                           const lar_pandora::HitsToMCParticles &trueHitsToParticles,
                                           lar_pandora::MCParticlesToPFParticles &matchedParticles,
                                           lar_pandora::MCParticlesToHits &matchedHits,
                                           PFParticleSet &vetoReco,
                                           MCParticleSet &vetoTrue,
                                           bool _recursiveMatching) 
{
    bool foundMatches(false);

    for (lar_pandora::PFParticlesToHits::const_iterator iter1 = recoParticlesToHits.begin(), iterEnd1 = recoParticlesToHits.end();
        iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> recoParticle = iter1->first;
        if (vetoReco.count(recoParticle) > 0)
            continue;

        const lar_pandora::HitVector &hitVector = iter1->second;

        lar_pandora::MCParticlesToHits truthContributionMap;

        for (lar_pandora::HitVector::const_iterator iter2 = hitVector.begin(), iterEnd2 = hitVector.end(); iter2 != iterEnd2; ++iter2)
        {
            const art::Ptr<recob::Hit> hit = *iter2;

            lar_pandora::HitsToMCParticles::const_iterator iter3 = trueHitsToParticles.find(hit);
            if (trueHitsToParticles.end() == iter3)
                continue;

            const art::Ptr<simb::MCParticle> trueParticle = iter3->second;
            if (vetoTrue.count(trueParticle) > 0)
                continue;

            truthContributionMap[trueParticle].push_back(hit);
        }

        lar_pandora::MCParticlesToHits::const_iterator mIter = truthContributionMap.end();

        for (lar_pandora::MCParticlesToHits::const_iterator iter4 = truthContributionMap.begin(), iterEnd4 = truthContributionMap.end();
            iter4 != iterEnd4; ++iter4)
        {
            if ((truthContributionMap.end() == mIter) || (iter4->second.size() > mIter->second.size()))
            {
                mIter = iter4;
            }
        }

        if (truthContributionMap.end() != mIter)
        {
            const art::Ptr<simb::MCParticle> trueParticle = mIter->first;

            lar_pandora::MCParticlesToHits::const_iterator iter5 = matchedHits.find(trueParticle);

            if ((matchedHits.end() == iter5) || (mIter->second.size() > iter5->second.size()))
            {
                matchedParticles[trueParticle] = recoParticle;
                matchedHits[trueParticle] = mIter->second;
                foundMatches = true;
            }
        }
    } // recoParticlesToHits loop ends

    if (!foundMatches)
        return;

    for (lar_pandora::MCParticlesToPFParticles::const_iterator pIter = matchedParticles.begin(), pIterEnd = matchedParticles.end();
        pIter != pIterEnd; ++pIter)
    {
        vetoTrue.insert(pIter->first);
        vetoReco.insert(pIter->second);
    }

    if (_recursiveMatching)
        GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedHits, vetoReco, vetoTrue, _recursiveMatching);
}


//______________________________________________________________________________________
void UBXSecHelper::GetTPCObjects(art::Event const & e, 
                                    std::string _particleLabel, 
                                    std::vector<lar_pandora::PFParticleVector> & pfp_v_v, 
                                    std::vector<lar_pandora::TrackVector> & track_v_v){

  //Vectors and maps we will use to store Pandora information
  lar_pandora::PFParticleVector pfParticleList;              //vector of PFParticles
  lar_pandora::PFParticlesToClusters pfParticleToClusterMap; //PFParticle-to-cluster map

  //Use LArPandoraHelper functions to collect Pandora information
  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _particleLabel, pfParticleList, pfParticleToClusterMap); //collect PFParticles and build map PFParticles to Clusters

  // Collect vertices and tracks
  lar_pandora::VertexVector           allPfParticleVertices;
  lar_pandora::PFParticlesToVertices  pfParticleToVertexMap;
  lar_pandora::LArPandoraHelper::CollectVertices(e, _particleLabel, allPfParticleVertices, pfParticleToVertexMap);
  lar_pandora::TrackVector            allPfParticleTracks;
  lar_pandora::PFParticlesToTracks    pfParticleToTrackMap;
  lar_pandora::LArPandoraHelper::CollectTracks(e, _particleLabel, allPfParticleTracks, pfParticleToTrackMap);

  GetTPCObjects(pfParticleList, pfParticleToTrackMap, pfParticleToVertexMap, pfp_v_v, track_v_v);
}



//____________________________________________________________
void UBXSecHelper::GetNuVertexFromTPCObject(art::Event const & e, 
                                               std::string _particleLabel,
                                               lar_pandora::PFParticleVector pfp_v, 
                                               double *reco_nu_vtx){

  reco_nu_vtx[0] = -9999;
  reco_nu_vtx[1] = -9999;
  reco_nu_vtx[2] = -9999;

  lar_pandora::VertexVector          vertexVector;
  lar_pandora::PFParticlesToVertices particlesToVertices;
  lar_pandora::LArPandoraHelper::CollectVertices(e, _particleLabel, vertexVector, particlesToVertices);

  for(unsigned int pfp = 0; pfp < pfp_v.size(); pfp++){

    if(lar_pandora::LArPandoraHelper::IsNeutrino(pfp_v.at(pfp))) {

      lar_pandora::VertexVector vertex_v = particlesToVertices.find(pfp_v.at(pfp))->second;
      if (vertex_v.size() > 1)
        std::cout << "More than one vertex associated to neutrino PFP!" << std::endl;
      else if (vertex_v.size() == 0)
        std::cout << "Zero vertices associated to neutrino PFP!" << std::endl;
      else {
       vertex_v[0]->XYZ(reco_nu_vtx);
       break;
      }
    }
  }

}



//_________________________________________________________________
art::Ptr<recob::PFParticle> UBXSecHelper::GetNuPFP(lar_pandora::PFParticleVector pfp_v){

  for (unsigned int pfp = 0; pfp < pfp_v.size(); pfp++) {

    if(lar_pandora::LArPandoraHelper::IsNeutrino(pfp_v.at(pfp))) {
      return pfp_v.at(pfp);
    }
  }
  std::cout << "No neutrino PFP found." << std::endl;

  art::Ptr<recob::PFParticle> temp;
  return temp;

}




//______________________________________________________________________________________
void UBXSecHelper::GetTPCObjects(lar_pandora::PFParticleVector pfParticleList, 
                                    lar_pandora::PFParticlesToTracks pfParticleToTrackMap, 
                                    lar_pandora::PFParticlesToVertices  pfParticleToVertexMap, 
                                    std::vector<lar_pandora::PFParticleVector> & pfp_v_v, 
                                    std::vector<lar_pandora::TrackVector> & track_v_v) {

  track_v_v.clear();
  pfp_v_v.clear();

  for (unsigned int n = 0; n < pfParticleList.size(); ++n) {
    const art::Ptr<recob::PFParticle> particle = pfParticleList.at(n);

    if(lar_pandora::LArPandoraHelper::IsNeutrino(particle)) {
      std::cout << "IS NEUTRINO, pfp id " << particle->Self() << std::endl;
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

      std::cout << "Start track_v.size() " << track_v.size() << std::endl;
      CollectTracksAndPFP(pfParticleToTrackMap, pfParticleList, particle, pfp_v, track_v);
      std::cout << "End   track_v.size() " << track_v.size() << std::endl;

      pfp_v_v.emplace_back(pfp_v);
      track_v_v.emplace_back(track_v);

      std::cout << "Number of pfp for this slice: "    << pfp_v.size()   << std::endl;
      std::cout << "Number of tracks for this slice: " << track_v.size() << std::endl;

      for (unsigned int i = 0; i < pfp_v.size(); i++) std::cout << "   pfp with ID " << pfp_v[i]->Self() << std::endl;
    } // end if neutrino
  } // end pfp loop

}



//______________________________________________________________________________________________________________________________________
void UBXSecHelper::CollectTracksAndPFP(lar_pandora::PFParticlesToTracks pfParticleToTrackMap,
                                          lar_pandora::PFParticleVector pfParticleList,
                                          art::Ptr<recob::PFParticle> particle,
                                          lar_pandora::PFParticleVector &pfp_v,
                                          lar_pandora::TrackVector &track_v) {

  pfp_v.emplace_back(particle);

  lar_pandora::PFParticlesToTracks::const_iterator trackMapIter = pfParticleToTrackMap.find(particle);
  if (trackMapIter != pfParticleToTrackMap.end()) {
    lar_pandora::TrackVector tracks = trackMapIter->second;
    for (unsigned int trk = 0; trk < tracks.size(); trk++) {
      track_v.emplace_back(tracks[trk]);
    }
  }
  std::cout << "Inter track_v.size() " << track_v.size() << std::endl;
  const std::vector<size_t> &daughterIDs = particle->Daughters();
  for (unsigned int d = 0; d < daughterIDs.size(); d++) std::cout << "daughter has id " << daughterIDs[d] << std::endl;
  if(daughterIDs.size() == 0) return;
  else {
    for (unsigned int m = 0; m < daughterIDs.size(); ++m) {
      const art::Ptr<recob::PFParticle> daughter = pfParticleList.at(daughterIDs.at(m));
      //pfp_v.emplace_back(daughter);
      CollectTracksAndPFP(pfParticleToTrackMap, pfParticleList, daughter, pfp_v, track_v);
    }
  }

}



//______________________________________________________________________________________________________________________________________
bool UBXSecHelper::InFV(double * nu_vertex_xyz){

  double x = nu_vertex_xyz[0];
  double y = nu_vertex_xyz[1];
  double z = nu_vertex_xyz[2];

  //This defines our current settings for the fiducial volume
  double FVx = 256.35;
  double FVy = 233;
  double FVz = 1036.8;
  double borderx = 10.;
  double bordery = 20.;
  double borderz = 10.;
  //double cryoradius = 191.61;
  //double cryoz = 1086.49 + 2*67.63;

  if(x < (FVx - borderx) && (x > borderx) && (y < (FVy/2. - bordery)) && (y > (-FVy/2. + bordery)) && (z < (FVz - borderz)) && (z > borderz)) return true;
  return false;

}



//__________________________________________________________________________
int UBXSecHelper::GetSliceOrigin(std::vector<art::Ptr<recob::PFParticle>> neutrinoOriginPFP, lar_pandora::PFParticleVector pfp_v) {

  bool isFromNu = false;
  int nuOrigin = 0;
  int cosmicOrigin = 0;

  for ( unsigned int i = 0; i < pfp_v.size(); i++) {

    for ( unsigned int j = 0; j < neutrinoOriginPFP.size(); j++) {

      if (neutrinoOriginPFP[j] == pfp_v[i]) {
        isFromNu = true;
        nuOrigin ++;
        break;
      }
    }

    if (!isFromNu) cosmicOrigin++;
    isFromNu = false;

  }

  if (nuOrigin > 0) return 0;
  else return 1;

}



//______________________________________________________________________________

void UBXSecHelper::GetNumberOfHitsPerPlane(art::Event const & e,
                                              std::string _particleLabel,
                                              lar_pandora::TrackVector track_v,
                                              int & nhits_u, 
                                              int & nhits_v, 
                                              int & nhits_w ) {
 
  nhits_u = 0;
  nhits_v = 0;
  nhits_w = 0;

  lar_pandora::TrackVector trackVector;
  lar_pandora::TracksToHits tracksToHits;
  lar_pandora::LArPandoraHelper::CollectTracks( e, _particleLabel, trackVector, tracksToHits );

  // Loop over the tracks in this TPC Object
  for (unsigned int t = 0; t < track_v.size(); t++) {

    // Get the hits associated with the track
    lar_pandora::HitVector hit_v = tracksToHits.at(track_v[t]);

    // Check where the hit is coming from
    for (unsigned int h = 0; h < hit_v.size(); h++){

      if (hit_v[h]->View() == 0) nhits_u++;
      if (hit_v[h]->View() == 1) nhits_v++;
      if (hit_v[h]->View() == 2) nhits_w++;

    }
  }

}

//_________________________________________________________________________________
bool UBXSecHelper::IsCrossingTopBoundary(recob::Track track, int & vtx_ok){

  double vtx[3];
  vtx[0] = track.Vertex().X();
  vtx[1] = track.Vertex().Y();
  vtx[2] = track.Vertex().Z();

  double end[3];
  end[0] = track.End().X();
  end[1] = track.End().Y();
  end[2] = track.End().Z();

  if ( InFV(vtx) && !InFV(end) && (end[1] > 233-20)) {
    vtx_ok = 0;
    return true;
  } else if (!InFV(vtx) && InFV(end) && (vtx[1] > 233-20))  {
    vtx_ok = 1;
    return true;
  }
  else {
    vtx_ok = -1;
    return false;
  }

}

//_________________________________________________________________________________
bool UBXSecHelper::IsCrossingBoundary(recob::Track track, int & vtx_ok){

  double vtx[3];
  vtx[0] = track.Vertex().X();
  vtx[1] = track.Vertex().Y();
  vtx[2] = track.Vertex().Z();

  double end[3];
  end[0] = track.End().X();
  end[1] = track.End().Y();
  end[2] = track.End().Z();

  if ( InFV(vtx) && !InFV(end) ) {
    vtx_ok = 0;
    return true;
  } else if (!InFV(vtx) && InFV(end))  {
    vtx_ok = 1;
    return true;
  }
  else {
    vtx_ok = -1;
    return false;
  }

}



//_________________________________________________________________________________
bool UBXSecHelper::GetLongestTrackFromTPCObj(lar_pandora::TrackVector track_v, recob::Track & out_track) {

  if (track_v.size() == 0) {
    return false;
  }

  int length = -1;
  int longest_track = -1;
  for (unsigned int t = 0; t < track_v.size(); t++){

    if (track_v[t]->Length() > length){

      length = track_v[t]->Length();
      longest_track = t;

    }

  }

  if (longest_track > -1){
    out_track = (*track_v[longest_track]);
    return true;
  } else {
    return false;
  }


}




//_________________________________________________________________________________
bool UBXSecHelper::PointIsCloseToDeadRegion(double *reco_nu_vtx, int plane_no){

  ::art::ServiceHandle<geo::Geometry> geo;

  // Get nearest channel
  raw::ChannelID_t ch = geo->NearestChannel(reco_nu_vtx, plane_no);
  std::cout << "nearest channel is " << ch << std::endl;

  // Get channel status
  const lariov::ChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
  std::cout << "ch status " << chanFilt.Status(ch) << std::endl;
  if( chanFilt.Status(ch) < 4) return true;

  // Now check close wires
  for(unsigned int new_ch = ch - 5; new_ch < ch - 5 + 11; new_ch++){
    std::cout << "now trying with channel " << new_ch << std::endl;
    if( chanFilt.Status(new_ch) < 4 ) return true;
  }

/*
  double new_point_in_tpc[3];
  for (int z_off = -5; z_off < 5; z_off += 2.5) {
    new_point_in_tpc[0] = reco_nu_vtx[0];
    new_point_in_tpc[1] = reco_nu_vtx[1];
    new_point_in_tpc[2] = reco_nu_vtx[2] + z_off;
    //std::cout << "trying with point " << new_point_in_tpc[0] << " " << new_point_in_tpc[1] << " " << new_point_in_tpc[2] << std::endl;

    raw::ChannelID_t ch = geo->NearestChannel(new_point_in_tpc, 2);
    if( chanFilt.Status(ch) < 4) return true;
  }
*/

  return false;

}









