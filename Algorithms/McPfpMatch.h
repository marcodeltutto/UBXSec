/**
 * \class McPfpMatch
 *
 * \ingroup UBXSec
 *
 * \brief Performs PFP to MC particles match
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


#ifndef MCPFPMATCH_H
#define MCPFPMATCH_H

#include <iostream>
#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "uboone/UBXSec/DataTypes/UBXSecFMWKInterface.h"

typedef std::map< art::Ptr<recob::PFParticle>, unsigned int > RecoParticleToNMatchedHits;
typedef std::map< art::Ptr<simb::MCParticle>,  RecoParticleToNMatchedHits > ParticleMatchingMap;
typedef std::set< art::Ptr<recob::PFParticle> > PFParticleSet;
typedef std::set< art::Ptr<simb::MCParticle> >  MCParticleSet;

namespace ubxsec {

  class McPfpMatch {
  
  public:

    /// Default constructor
    McPfpMatch();
  
    /// Default destructor
    ~McPfpMatch(){}

    /// Configure
    void _Configure_(const Config_t &pset);
     
    /// Configure function parameters
    void Configure(art::Event const & e, std::string _pfp_producer, std::string _spacepointLabel, std::string _hitfinderLabel, std::string _geantModuleLabel);
  
     /**
     *  @brief Returns matching between true and reconstructed particles
     *
     *  @param matchedParticles the output matches between reconstructed and true particles
     *  @param matchedHits the output matches between reconstructed particles and hits
     */  
    void GetRecoToTrueMatches(lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits);
  

  protected:

    /**
     *  @brief Perform matching between true and reconstructed particles
     *
     *  @param recoParticlesToHits the mapping from reconstructed particles to hits
     *  @param trueHitsToParticles the mapping from hits to true particles
     *  @param matchedParticles the output matches between reconstructed and true particles
     *  @param matchedHits the output matches between reconstructed particles and hits
     */
    void GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, const lar_pandora::HitsToMCParticles &trueHitsToParticles, lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits);

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
    void GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, const lar_pandora::HitsToMCParticles &trueHitsToParticles, lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits, PFParticleSet &recoVeto, MCParticleSet &trueVeto, bool _recursiveMatching);

    lar_pandora::HitsToMCParticles _trueHitsToParticles; ///< A map from recon hits to MCParticles
    lar_pandora::PFParticlesToHits _recoParticlesToHits; ///< A map from PFParticles to recon hits

  private:

    bool _debug             = false;
    bool _verbose           = false;
    bool _recursiveMatching = false;

  };
}

#endif //  MCPFPMATCH_H
