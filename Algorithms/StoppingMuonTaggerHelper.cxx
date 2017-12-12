#ifndef STOPPINGMUONTAGGERHELPER_CXX
#define STOPPINGMUONTAGGERHELPER_CXX

#include "StoppingMuonTaggerHelper.h"
#include <iostream>
#include <algorithm>
#include <numeric> 
#include "TVector3.h"

namespace ubana {

  //_________________________________________________________________________________ 
  StoppingMuonTaggerHelper::StoppingMuonTaggerHelper()
  {
    if (_debug) std::cout << "StoppingMuonTaggerHelper Instantiated." << std::endl;
  }

  //_________________________________________________________________________________ 
  void StoppingMuonTaggerHelper::Configure(fhicl::ParameterSet const& pset)
  {

    _w2cm = pset.get< double > ( "WireToCmConstant", 0.3 );
    _t2cm = pset.get< double > ( "TimeToCmConstant", 0.0557 ); 
    _dqds_calib = pset.get< double > ( "GainCalib", 198 );
    _slider_window = pset.get< size_t > ( "SliderWindow", 10 );
    _max_allowed_hit_distance = pset.get< double > ( "MaxAllowedHitDistance", 15 ); 
    _hits_to_remove = pset.get< int > ( "HitsToRemove", 3 );
    _pre_post_window = pset.get< int > ( "PrePostWindow", 5 ); 
    _perc_diff_cut = pset.get< double > ( "PercDiffCut", 20 );
    _debug = pset.get< bool > ( "DebugMode", false );

  }

  //_________________________________________________________________________________ 
  void StoppingMuonTaggerHelper::PrintConfig() {

    std::cout << "--- StoppingMuonTaggerHelper configuration:" << std::endl;
    //std::cout << "---   _max_distance  = " << _max_distance << std::endl;

  }

  //_________________________________________________________________________________ 
  void StoppingMuonTaggerHelper::Clear(){
     _s_hit_v.clear();
     _dqds_v.clear();
     _ds_v.clear();
     _dqds_slider.clear();

     _start_index = -1;
     _hits_ordered = false;
  }


  //_________________________________________________________________________________ 
  size_t StoppingMuonTaggerHelper::FilterByPlane(int plane_no) {

    std::vector<ubana::SimpleHit> new_vector;
    new_vector.clear();

    for (size_t element = 0; element < _s_hit_v.size(); element++) {

      if (_s_hit_v.at(element).plane == plane_no)
        new_vector.push_back(_s_hit_v.at(element));

    }

    std::swap(new_vector, _s_hit_v);

    return _s_hit_v.size();

  }


  //_________________________________________________________________________________
  void StoppingMuonTaggerHelper::SetStartHit(double time, 
                                             double wire_no, 
                                             int plane_no) {


    if (plane_no != 2) {
      std::cout << "Plane not supported." << std::endl;
      return;
    }

    TVector3 pt1(time, wire_no, 0);
    
    double min_dist = 1e9;
    int best_hit_id = -1;

    if (_debug) std::cout << "Simple hit vector size " << _s_hit_v.size() << std::endl;

    for (size_t i = 0; i < _s_hit_v.size(); i++) {

      auto sh = _s_hit_v.at(i);
      TVector3 pt2(sh.time, sh.wire, 0);
      double dist = (pt1-pt2).Mag();

      if (dist< min_dist) {
        best_hit_id = i;
        min_dist = dist;
      }
    }

    if (best_hit_id == -1) {
      std::cout << "Could not find start hit." << std::endl;
      return;
    }
   
 
    // The best hit may not be the start one, but one very close,
    // So let's find a border hit

    auto almost_best_hit = _s_hit_v.at(best_hit_id);

    if (_debug) std::cout << "Almost best hit has wire " << almost_best_hit.wire
                          << ", and time " << almost_best_hit.time*4 << std::endl;

    TVector3 pt0(almost_best_hit.time, almost_best_hit.wire, 0);

    // First create a map from wire to hit
    std::map<int,ubana::SimpleHit> wire_to_hit;

    for (size_t i = 0; i < _s_hit_v.size(); i++) {

      auto sh = _s_hit_v.at(i); 

      // If we never encountered this wire, just save it
      auto iter = wire_to_hit.find(sh.wire);
      if (iter == wire_to_hit.end()) {
        wire_to_hit[sh.wire] = sh;
        continue;
      }

      // Otherwise pick the one that is closer to the almost_best_hit
      auto previous_hit = iter->second;
      TVector3 pt1(previous_hit.time, previous_hit.wire, 0);
      double dist_previous = (pt0-pt1).Mag();
      TVector3 pt2(sh.time, sh.wire, 0);
      double dist_current = (pt0-pt2).Mag();

      if (dist_current < dist_previous)
        wire_to_hit[sh.wire] = sh;
      else 
        continue;
    }

    // Then try going to the left first
    int n_step_left = 0;
    int best_index_left = -1;
    for (int w = almost_best_hit.wire - 1; w > 0; w--) {

      auto iter = wire_to_hit.find(w);
      if (iter == wire_to_hit.end()) {
        break;
      }

      //ubana::SimpleHit temp = iter->second;
      auto it = std::find(_s_hit_v.begin(), _s_hit_v.end(), iter->second);
      best_index_left = it - _s_hit_v.begin();//w;//std::distance(_s_hit_v.begin(), it);
      TVector3 pt1(iter->second.time, iter->second.wire, 0);
      double dist = (pt0-pt1).Mag();
  
      if (dist > _max_allowed_hit_distance) // Check this number! 
        break;

      n_step_left ++;

      if (_debug) std::cout << "Found hit on the left, wire " << iter->second.wire
                            << ", time " << iter->second.time << std::endl;
    }

    // Then go rigth
    int n_step_right = 0;
    int best_index_right = -1;
    for (int w = almost_best_hit.wire + 1; w < 3456; w++) {

      auto iter = wire_to_hit.find(w);
      if (iter == wire_to_hit.end()) {
        break;
      }

      auto it = std::find(_s_hit_v.begin(), _s_hit_v.end(), iter->second);
      best_index_right = it - _s_hit_v.begin();//w;//std::distance(_s_hit_v.begin(), it);
      TVector3 pt1(iter->second.time, iter->second.wire, 0);
      double dist = (pt0-pt1).Mag(); 

      if (dist > _max_allowed_hit_distance) // Check this number!  
        break; 

      n_step_right ++;

      if (_debug) std::cout << "Found hit on the right, wire " << iter->second.wire 
                            << ", time " << iter->second.time*4 << std::endl; 
    }

    _start_index = 0;

    if (n_step_left == 0 || n_step_right == 0) {
      _start_index = best_hit_id;
    } else if (n_step_left > n_step_right) {
      _start_index = best_index_right;
    } else if (n_step_right >= n_step_left) {
      _start_index = best_index_left;
    }


    if (_debug) std::cout << "[StoppingMuonTaggerHelper] Start hit set to w: " 
                          << _s_hit_v.at(_start_index).wire << ", and t: " 
                          << _s_hit_v.at(_start_index).time*4 << std::endl;


    return;
  }

  

  //_________________________________________________________________________________ 
  void StoppingMuonTaggerHelper::OrderHits() {

    if (_start_index < 0) {
      std::cout << "Start hit not set." << std::endl;
      return;
    }

    if ((size_t)_start_index >= _s_hit_v.size()) {
      std::cout << "Start hit not compatible with current hit vector." 
       << "Start hit index is " << _start_index 
       << ", while hit vector size is " << _s_hit_v.size() << std::endl;
      return;
    }


    ubana::SimpleHitVector new_vector;
    new_vector.clear();
    new_vector.reserve(_s_hit_v.size());

    _ds_v.clear();
    _ds_v.reserve(_s_hit_v.size());

    new_vector.push_back(_s_hit_v.at(_start_index));

    if (_debug) {
      for (auto h : _s_hit_v) {
        std::cout << "BEFORE: " << h.wire << ", " << h.time*4 << std::endl;
      }
    }

    _s_hit_v.erase(_s_hit_v.begin() + _start_index);

    double min_dist = 1e9; 
    int min_index = -1;

    while (_s_hit_v.size() != 0) {

      min_dist = 1e9;
      min_index = -1; 

      for (size_t i = 0; i < _s_hit_v.size(); i++){

        TVector3 pt1(new_vector.back().time, new_vector.back().wire, 0);
        TVector3 pt2(_s_hit_v.at(i).time,    _s_hit_v.at(i).wire,    0);
        double dist = (pt1 - pt2).Mag();

        if (dist < min_dist) {
          min_index = i;
          min_dist = dist;
        }
      }

      if (min_index < 0) {
         std::cout << "Logic fail." << std::endl;
         throw std::exception();
      }

      // Emplace the next hit in the new vector...
      if (min_dist < _max_allowed_hit_distance) {
        if (_debug) std::cout << "min_dist: " << min_dist <<std::endl;
        new_vector.push_back(_s_hit_v.at(min_index));
        _ds_v.push_back(min_dist);
      }

      // ...and delete it from the old vector
      _s_hit_v.erase(_s_hit_v.begin() + min_index);

    }

    // For the last hit, use the last values of min_dist
    _ds_v.push_back(min_dist);

    // Now that the vector is ordered, reassing to the original one
    _s_hit_v = new_vector;

    if (_debug) {
      for (auto h : _s_hit_v) {
        std::cout << "AFTER: " << h.wire << ", " << h.time*4 << std::endl;
      }
    }

    _hits_ordered = true;

    return;
  }

  void StoppingMuonTaggerHelper::CalculatedQds() {

    if (!_hits_ordered) {
      std::cout << "Call to " << __FUNCTION__ << " without having ordered hits." << std::endl;
      throw std::exception();
    }
 
    if (_ds_v.size() != _s_hit_v.size()) {
      std::cout << __FUNCTION__ << ": ds size is different than hit vector size" << std::endl;
      throw std::exception();
    }
 
    _dqds_v.clear();
    _dqds_v.reserve(_s_hit_v.size());

    double ds = 1.;

    for (size_t i = 0; i < _s_hit_v.size()-1; i++) {

      TVector3 this_point(_s_hit_v.at(i).wire * _w2cm, _s_hit_v.at(i).time * _t2cm, 0);
      TVector3 next_point(_s_hit_v.at(i+1).wire * _w2cm, _s_hit_v.at(i+1).time * _t2cm, 0);
      ds = (this_point - next_point).Mag();

      _dqds_v.emplace_back(_s_hit_v.at(i).integral/ds * _dqds_calib);
      if (_debug) std::cout << "_dqds_v.back() " << _dqds_v.back() << std::endl;

    }

    // Finish with the last point
    _dqds_v.emplace_back(_s_hit_v.at(_s_hit_v.size()-1).integral/ds * _dqds_calib);

    return;
  }



  void StoppingMuonTaggerHelper::PerformdQdsSlider() {

    if (_dqds_v.size() < _slider_window * 2)
      return;

    _dqds_slider.clear();

    size_t window = _slider_window;

    for (size_t i = 0; i < _dqds_v.size() - window; i++) {

      std::vector<double> temp_v(_dqds_v.begin() + i, _dqds_v.begin() + i + window);
      double median_dqds = this->GetTruncMedian(temp_v);
      _dqds_slider.emplace_back(median_dqds);
      if (_debug) std::cout << "dqds_slider value " << median_dqds << std::endl;

    }

    return;

  }

  double StoppingMuonTaggerHelper::GetTruncMedian(std::vector<double> v) {

    // Find and erase max element
    auto it_max = std::max_element(v.begin(), v.end());
    v.erase(it_max);

    // Find and erase min element
    auto it_min = std::min_element(v.begin(), v.end());
    v.erase(it_min);

    double median = -1;

    size_t size = v.size();
    std::sort(v.begin(), v.end());
    if (size % 2 == 0){
      median = (v[size/2 - 1] + v[size/2]) / 2;
    }
    else{
      median = v[size/2];
    }

    return median;
  }


  bool StoppingMuonTaggerHelper::MakeDecision(::ubana::StopMuAlgoType algo) {

   switch (algo) {

     case kAlgoMichel:
       return this->IsStopMuMichel();
       break;

     case kAlgoBragg:
       return this->IsStopMuBragg();
       break;

     case kAlgoSimpleMIP:
       return false;
       break;

     case kAlgoCurvature:
       return false;
       break;

     default:
       std::cout << "Algo not found." << std::endl;
       throw std::exception();
    }

    return false;

  }


  bool StoppingMuonTaggerHelper::IsStopMuMichel() {

    if (_dqds_slider.size() < (unsigned int) (_hits_to_remove * 2 + _pre_post_window * 2)) {
      if (_debug) std::cout << "Can't make decision, number of simple hits is " << _dqds_slider.size() << ", which is less then " << _hits_to_remove * 2 + _pre_post_window * 2 << std::endl;
      return false;
    }

    std::vector<double> temp = _dqds_slider;

    // Remove first three and last three hits
    temp.erase(temp.begin(), temp.begin() + _hits_to_remove);
    temp.erase(temp.end() - _hits_to_remove, temp.end());

    // Get mean of first and last 5
    double start_mean = std::accumulate(temp.begin(), temp.begin() + _pre_post_window, 0);
    start_mean /= _pre_post_window;
    double end_mean = std::accumulate(temp.end() - _pre_post_window, temp.end(), 0);
    end_mean /= _pre_post_window;

    double perc_diff = (start_mean - end_mean) / start_mean * 100.;

    if (_debug) std::cout << "[IsStopMuBragg] Start mean: " << start_mean << ", end mean " << end_mean << ", Perc diff is " << perc_diff << std::endl;

    if (perc_diff > _perc_diff_cut) {
      return true;
    } 

    return false;   
  }




  bool StoppingMuonTaggerHelper::IsStopMuBragg() {

    if (_dqds_slider.size() < (unsigned int) (_pre_post_window * 2)) {
      if (_debug) std::cout << "Can't make decision, number of simple hits is " << _dqds_slider.size() << ", which is less then " << _pre_post_window * 2 << std::endl;
      return false;
    }

    // Vertex must not be in the FV
    if (_fv.InFV(_vertex))
      return false;

    // Take firsts hits, then lasts hits 
    double start_mean = std::accumulate(_dqds_slider.begin(), _dqds_slider.begin() + _pre_post_window, 0);
    start_mean /= _pre_post_window;
    double end_mean = std::accumulate(_dqds_slider.end() - _pre_post_window, _dqds_slider.end(), 0);
    end_mean /= _pre_post_window;

    double perc_diff = (start_mean - end_mean) / start_mean * 100.;

    if (_debug) std::cout << "[IsStopMuBragg] Start mean: " << start_mean << ", end mean " << end_mean << ", Perc diff is " << perc_diff << std::endl;

    if (perc_diff < -_perc_diff_cut) {
      return true;
    }

    return false;

  }


}


#endif
