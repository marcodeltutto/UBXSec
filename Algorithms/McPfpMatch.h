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
 * Heavily updated on: Friday 14, 2017 at 16:52:41
 *
 */


#ifndef MCPFPMATCH_H
#define MCPFPMATCH_H

#include <iostream>
#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "uboone/UBXSec/DataTypes/UBXSecFMWKInterface.h"

namespace lar_pandora { 
  typedef std::map<art::Ptr<recob::PFParticle>, art::Ptr<simb::MCParticle>> PFParticlesToMCParticles;
}

namespace ubana {

  class McPfpMatch {
  
  public:

    /// Default constructor
    McPfpMatch();
  
    /// Default destructor
    ~McPfpMatch(){}

    /// Configure
    void _Configure_(const Config_t &pset);
     
    /// Configure function parameters
     /**
     *  @brief Configure function parameters (call this function first)
     *
     *  @param e the art::Event
     *  @param _pfp_producer the PFParticle producer label
     *  @param _spacepoint_producer the SpacePoint producer label
     *  @param _hitfinder_producer the Hit producer label
     *  @param _geant_producer The Geant4 producer label
     */
    void Configure(art::Event const & e, std::string _pfp_producer, std::string _spacepoint_producer, std::string _hitfinder_producer, std::string _geant_producer);

    /// Configure function parameters
     /**
     *  @brief Configure function parameters (call this function first)
     *
     *  @param e the art::Event
     *  @param _pfp_producer the PFParticle producer label
     *  @param _spacepoint_producer the SpacePoint producer label
     *  @param _hitfinder_producer the Hit producer label
     *  @param _geant_producer The Geant4 producer label
     *  @param _hit_mcp_producer The producer that created hit <-> MCP associations (the final news of MCC8)
     *  @param daughterMode
     */
    void Configure(art::Event const & e, std::string _pfp_producer, std::string _spacepoint_producer, std::string _hitfinder_producer, std::string _geant_producer, std::string _hit_mcp_producer, lar_pandora::LArPandoraHelper::DaughterMode daughterMode );
  
     /**
     *  @brief Returns matching between true and reconstructed particles
     *
     *  @param matchedParticles the output matches between reconstructed and true particles
     */  
    void GetRecoToTrueMatches(lar_pandora::PFParticlesToMCParticles & matchedParticles);

    /**
     *  @brief If calles
     *
     *  @param option It true, considers the event as MC always
     */
    void OverrideRealData(bool option) {_override_real_data = option;};
  

  protected:

    lar_pandora::HitsToMCParticles _hit_to_mcps_map; ///< A map from recon hits to MCParticles
    lar_pandora::PFParticlesToHits _pfp_to_hits_map; ///< A map from PFParticles to recon hits

  private:

    bool _configured = false;

    bool _debug      = false;
    bool _verbose    = false;

    bool _is_data    = false; ///< If true, we are running over a real data file.
    bool _override_real_data = false;

  };
}

#endif //  MCPFPMATCH_H
