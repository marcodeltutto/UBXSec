#ifndef NUMUCCEVENTSELECTION_CXX
#define NUMUCCEVENTSELECTION_CXX

#include "NuMuCCEventSelection.h"
#include <iostream>

namespace ubana {

  NuMuCCEventSelection::NuMuCCEventSelection() {
    _configured = false;
  }

  void NuMuCCEventSelection::Configure(fhicl::ParameterSet const& pset)
  {
    _deltax_cut_down          = pset.get< double > ( "DeltaXCutDown" );
    _deltax_cut_up            = pset.get< double > ( "DeltaXCutUp" );
    _deltaz_cut_down          = pset.get< double > ( "DeltaZCutDown" );
    _deltaz_cut_up            = pset.get< double > ( "DeltaZCutUp" );
    _flsmatch_score_cut       = pset.get< double > ( "FlsMatchScoreCut" );
    _vtxcheck_angle_cut_down  = pset.get< double > ( "VtxCheckAngleCutDown" );
    _vtxcheck_angle_cut_up    = pset.get< double > ( "VtxCheckAngleCutUp" );
    _ntrack_cut               = pset.get< double > ( "NTrackCut" );
    _pe_cut                   = pset.get< double > ( "PECut" );
    _beamSpillStarts          = pset.get< double > ( "BeamSpillStarts" );
    _beamSpillEnds            = pset.get< double > ( "BeamSpillEnds" );

    _verbose                  = pset.get< bool > ( "Verbose" );

    _configured = true;
  }

  void NuMuCCEventSelection::PrintConfig() {

    std::cout << "--- NuMuCCEventSelection configuration:" << std::endl;
    std::cout << "---   _deltax_cut_down          = " << _deltax_cut_down << std::endl;
    std::cout << "---   _deltax_cut_up            = " << _deltax_cut_up << std::endl;
    std::cout << "---   _deltaz_cut_down          = " << _deltaz_cut_down << std::endl;
    std::cout << "---   _deltaz_cut_up            = " << _deltaz_cut_up << std::endl;
    std::cout << "---   _flsmatch_score_cut       = " << _flsmatch_score_cut << std::endl;
    std::cout << "---   _vtxcheck_angle_cut_down  = " << _vtxcheck_angle_cut_down << std::endl;
    std::cout << "---   _vtxcheck_angle_cut_up    = " << _vtxcheck_angle_cut_up << std::endl;
    std::cout << "---   _ntrack_cut               = " << _ntrack_cut << std::endl;
    std::cout << "---   _pe_cut                   = " << _pe_cut << std::endl;
    std::cout << "---   _beamSpillStarts          = " << _beamSpillStarts << std::endl;
    std::cout << "---   _beamSpillEnds            = " << _beamSpillEnds << std::endl;
    std::cout << "---   _verbose                  = " << _verbose << std::endl;
  }

  void NuMuCCEventSelection::SetEvent(UBXSecEvent* ubxsec_event) {

    _ubxsec_event = ubxsec_event;
    _event_is_set = true;

  }

  bool NuMuCCEventSelection::IsSelected(size_t & slice_index, std::string & reason) {
    
    if (!_configured || !_event_is_set) {
      std::cout << "Call to NuMuCCEventSelection::IsSelected() without having done configuration or having set the UBXSecEvent" << std::endl;
      throw std::exception();
    }

    slice_index = -1;
    reason = "no_failure";

    // ************ 
    // Flashes
    // ************

    if (_ubxsec_event->nbeamfls == 0) {
      reason = "no_beam_disc_flashes";
      return false;
    }

    int flashInBeamSpill = -1;    
    double old_pe = -1;

    for (int fls = 0; fls < _ubxsec_event->nbeamfls; fls ++){

      if (_ubxsec_event->beamfls_time.at(fls) > _beamSpillStarts 
          && _ubxsec_event->beamfls_time.at(fls) < _beamSpillEnds
          && _ubxsec_event->beamfls_pe.at(fls) >= _pe_cut) {
        
        if (_ubxsec_event->beamfls_pe.at(fls) > old_pe) {
          flashInBeamSpill = fls;
          old_pe = _ubxsec_event->beamfls_pe.at(fls);
        }    
      }
    }
    
    if (flashInBeamSpill == -1) {
      reason = "no_beam_spill_flash";
      return false;   
    }
    
    if (_verbose) std::cout << "[NuMuCCEventSelection] Pass Flash" << std::endl;

    // *******************************
    // Find slice with maximum score
    // *******************************
    double score_max = -1;
    int scl_ll_max = -1;
    for (int slc = 0; slc < _ubxsec_event->nslices; slc ++){
      
      if(_ubxsec_event->slc_flsmatch_score.at(slc) > score_max){
        scl_ll_max = slc;
        score_max = _ubxsec_event->slc_flsmatch_score.at(slc);
      }

    }
    
    if (scl_ll_max == -1) {
      reason = "fail_flash_match";
      return false;
    }
    
    if (score_max <= _flsmatch_score_cut) {
      reason = "fail_flash_match_score";
      return false;
    }

    if (_verbose) std::cout << "[NuMuCCEventSelection] Pass FlashMatching Score" << std::endl;
 
    // Delta X
    if(_ubxsec_event->slc_flsmatch_qllx.at(scl_ll_max) - _ubxsec_event->slc_flsmatch_tpcx.at(scl_ll_max) > _deltax_cut_up) {
      reason = "fail_flash_match_deltax_up";
      return false;
    }
    if(_ubxsec_event->slc_flsmatch_qllx.at(scl_ll_max) - _ubxsec_event->slc_flsmatch_tpcx.at(scl_ll_max) < _deltax_cut_down) {
      reason = "fail_flash_match_deltax_down";
      return false;
    }
    if (_verbose) std::cout << "[NuMuCCEventSelection] Pass FlashMatching DeltaX" << std::endl;

    // DeltaZ
    if(_ubxsec_event->slc_flsmatch_hypoz.at(scl_ll_max) - _ubxsec_event->beamfls_z.at(flashInBeamSpill) > _deltaz_cut_up) {
      reason = "fail_flash_match_deltaz_up";
      return false;
    }
    if(_ubxsec_event->slc_flsmatch_hypoz.at(scl_ll_max) - _ubxsec_event->beamfls_z.at(flashInBeamSpill) < _deltaz_cut_down) {
      reason = "fail_flash_match_deltaz_down";
      return false;
    }
    if (_verbose) std::cout << "[NuMuCCEventSelection] Pass FlashMatching DeltaZ" << std::endl;

    // FV
    if(_ubxsec_event->slc_nuvtx_fv.at(scl_ll_max) == 0) {
      reason = "fail_fiducial_volume";
      return false;
    }
    if (_verbose) std::cout << "[NuMuCCEventSelection] Pass Fiducial Volume" << std::endl;

    // Vertex Check   
    if(_ubxsec_event->slc_vtxcheck_angle.at(scl_ll_max) > _vtxcheck_angle_cut_up) {
      reason = "fail_vertex_check_up";
      return false;
    }
    if(_ubxsec_event->slc_vtxcheck_angle.at(scl_ll_max) < _vtxcheck_angle_cut_down 
       && _ubxsec_event->slc_vtxcheck_angle.at(scl_ll_max) !=-9999 ) {
      reason = "fail_vertex_check_down";
      return false;
    } 
    if (_verbose) std::cout << "[NuMuCCEventSelection] Pass VertexCheck" << std::endl;

    // N Track
    if(_ubxsec_event->slc_ntrack.at(scl_ll_max) < _ntrack_cut) {
      reason = "fail_ntrack";
      return false;
    }
    if (_verbose) std::cout << "[NuMuCCEventSelection] Pass NTrack Requirement" << std::endl; 

    // Track Quality
    if(!_ubxsec_event->slc_passed_min_track_quality.at(scl_ll_max)) {
      reason = "fail_track quality";
      return false;
    }
    if (_verbose) std::cout << "[NuMuCCEventSelection] Pass Track Quality" << std::endl;
    
    // Vertex Quality
    if(!_ubxsec_event->slc_passed_min_vertex_quality.at(scl_ll_max)) {
      reason = "fail_vertex_quality";
      return false;
    }   
    if (_verbose) std::cout << "[NuMuCCEventSelection] Pass Vertex Quality" << std::endl;

    slice_index = scl_ll_max;
    return true;

  }
}


#endif
