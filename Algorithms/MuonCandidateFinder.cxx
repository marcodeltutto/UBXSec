#ifndef MUONCANDIDATEFINDER_CXX
#define MUONCANDIDATEFINDER_CXX

#include "MuonCandidateFinder.h"
#include <iostream>

namespace ubana {

  MuonCandidateFinder::MuonCandidateFinder()
  {
    _use_pida_cut   = false;
    _tracks_are_set = false;
    _tracktopidmap_is_set = false;
  }

  void MuonCandidateFinder::Configure(fhicl::ParameterSet const& pset)
  {
    _use_pida_cut      = pset.get<bool> ( "UsePIDACut", false );

    _gain_calib        = pset.get<double> ("GainCalib", 1);

    _svm_x             = pset.get<std::vector<double>> ("SVM_X");
    _alpha             = pset.get<std::vector<double>> ("SVM_alpha");
    _support_vectors_x = pset.get<std::vector<double>> ("SVM_SupportVectorX");
    _support_vectors_y = pset.get<std::vector<double>> ("SVM_SupportVectorY");
    _rho               = pset.get<double> ("SVM_Rho");
    _loop              = pset.get<double> ("SVM_Loop");
    _gamma             = pset.get<double> ("SVM_Gamma");
    _degree            = pset.get<double> ("SVM_Degree");
    _r                 = pset.get<double> ("SVM_R");
    _x1_mean           = pset.get<double> ("SVM_X1Mean");
    _x2_mean           = pset.get<double> ("SVM_X2Mean");
    _x1_range          = pset.get<double> ("SVM_X1Range");
    _x2_range          = pset.get<double> ("SVM_X2Range");

    if (_svm_x.size() != 1000) {
      std::cout << "[MuonCandidateFinder] _svm_x size error, is " << _svm_x.size() << ", should be 1000" << std::endl;
      throw std::exception();
    }

    if (_alpha.size() != _support_vectors_x.size() 
      && _support_vectors_x.size() != _support_vectors_y.size() 
      && _support_vectors_y.size() != _loop) {
      std::cout << "[MuonCandidateFinder] Vector size mismatch: " 
                << "alpha: " << _alpha.size()
                << "sv x: " << _support_vectors_x.size()
                << "sv y: " << _support_vectors_y.size() 
                << "loop: " << _loop << std::endl;
      throw std::exception();
    }
  }

  void MuonCandidateFinder::PrintConfig() {

    std::cout << "--- MuonCandidateFinder configuration:" << std::endl;
    std::cout << "---   _use_pida_cut  = " << _use_pida_cut << std::endl;
    std::cout << "---   _gain_calib    = " << _gain_calib << std::endl;
    std::cout << "---   _rho           = " << _rho << std::endl;
    std::cout << "---   _loop          = " << _loop << std::endl;
    std::cout << "---   _gamma         = " << _gamma << std::endl;
    std::cout << "---   _degree        = " << _degree << std::endl;
    std::cout << "---   _r             = " << _r << std::endl;
    std::cout << "---   _x1_mean       = " << _x1_mean << std::endl;
    std::cout << "---   _x2_mean       = " << _x2_mean << std::endl;
    std::cout << "---   _x1_range      = " << _x1_range << std::endl;
    std::cout << "---   _x2_range      = " << _x2_range << std::endl;

  }

  void MuonCandidateFinder::Reset() 
  {
    _tracks_are_set = false;
    _tracktopidmap_is_set = false;    
  }

  bool MuonCandidateFinder::GetCandidateTrack(art::Ptr<recob::Track> & out_track) {
   
    if (!_tracktopidmap_is_set || !_tracks_are_set) {
      std::cerr << "[MuonCandidateFinder] Call to GetCandidateTrack without having set input Tracks and map to PIDs. Exiting now." << std::endl;
      throw std::exception();
    }

    if (_tracks.size() == 0)
      return false;

    std::vector<art::Ptr<recob::Track>> possible_tracks;

    if (_use_pida_cut) {
      for (auto track : _tracks) {

        double length, pida;

        length = track->Length();

        auto iter = _track_to_pid.find(track);
        if(iter != _track_to_pid.end()) {
          pida = iter->second->PIDA();
        } else {
          possible_tracks.emplace_back(track);
          continue;
        }

        if (length < 90. && pida > 13) {
          continue;
        } else {
          possible_tracks.emplace_back(track);
        }
      }

      if (possible_tracks.size() == 0) {
        return false;
      }else if (possible_tracks.size() == 1) {
        out_track = possible_tracks.at(0);
        return true;
      }
    }
    else {
      possible_tracks = _tracks;
    }

    // Sort remaining tracks by length
    std::sort(possible_tracks.begin(), possible_tracks.end(),
              [](art::Ptr<recob::Track> a, art::Ptr<recob::Track> b) -> bool
              {
                return a->Length() > b->Length();
              });

    out_track = possible_tracks.at(0);

    return true;
  }



  bool MuonCandidateFinder::MIPConsistency(double dqds, double length) {

    if (length > 1000)
      return true;

    if (length <= 0) {
      std::cout << "[MuonCandidateFinder] Track length <= 0?!" << std::endl;
      return false;
    }

    dqds *= _gain_calib;

    int l = std::round(length);
    double dqds_cut = _svm_x.at(l);

    std::cout << "[MuonCandidateFinder] Track length is " << length << " cm, dqds_cut is " << dqds_cut << " e-/cm" << std::endl;

    if (dqds <= dqds_cut)
      return true;
  

    return false;

  }

  
  bool MuonCandidateFinder::SVMPredict(double dqds, double length) {

    dqds *= _gain_calib;

    // Feature scaling first
    std::pair<double, double> user_vector = this->_SVM_feature_scaling(dqds, length);

    // Check decision from classifier function
    double decision_argument = 0;

    for (size_t i = 0; i < _loop; i++) {

      decision_argument += _alpha.at(i) * this->_SVM_kernel(std::make_pair(_support_vectors_x.at(i), _support_vectors_y.at(i)), user_vector);
    }

    decision_argument += _rho;

    return (decision_argument > 0);

  }

  double MuonCandidateFinder::_SVM_kernel(std::pair<double, double> support_vector, std::pair<double, double> user_vector) {

    double result = support_vector.first * user_vector.first + support_vector.second * user_vector.second;

    result *= _gamma;

    result += _r;

    result = pow(result, _degree);

    return result;
  }

  std::pair<double, double>  MuonCandidateFinder::_SVM_feature_scaling(double dqds, double length) {


    double new_dqds = ((dqds - _x1_mean) / _x1_range + 1) * 10;
    double new_length = ((length - _x2_mean) / _x2_range + 1) * 10;

    return std::make_pair(new_dqds, new_length);

  }
  
}


#endif
