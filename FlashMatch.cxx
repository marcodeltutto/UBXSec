#include "FlashMatch.h"
#include <vector>

namespace ubana {

  FlashMatch::FlashMatch() {
    fScore      = -9999;
    fTPCX       = -9999;
    fEstimatedX = -9999;

    fXFixedChi2 = -9999;
    fXFixedLl   = -9999;
  }

  FlashMatch::FlashMatch(double score = -9999) {
    fScore      = score;
    fTPCX       = -9999;
    fEstimatedX = -9999;

    fXFixedChi2 = -9999;
    fXFixedLl   = -9999;
  }
 
  FlashMatch::~FlashMatch(){
  }

  // Setter methoths
  void FlashMatch::SetScore               (double score)             { this->fScore = score;              }
  void FlashMatch::SetTPCX                (double tpcx)              { this->fTPCX = tpcx;                }
  void FlashMatch::SetEstimatedX          (double estx)              { this->fEstimatedX = estx;          }
  void FlashMatch::SetHypoFlashSpec       (std::vector<double> spec) { this->fHypoFlashSpec = spec;       } 
  void FlashMatch::SetRecoFlashSpec       (std::vector<double> spec) { this->fRecoFlashSpec = spec;       }
  void FlashMatch::SetMCFlashSpec         (std::vector<double> spec) { this->fMCFlashSpec = spec;         }
  void FlashMatch::SetXFixedHypoFlashSpec (std::vector<double> spec) { this->fXFixedHypoFlashSpec = spec; }
  void FlashMatch::SetXFixedChi2          (double chi2)              { this->fXFixedChi2 = chi2;          }
  void FlashMatch::SetXFixedLl            (double ll)                { this->fXFixedLl = ll;              }


  // Getter methods
  const double              & FlashMatch::GetScore()               const { return this->fScore;               }
  const double              & FlashMatch::GetTPCX()                const { return this->fTPCX;                }
  const double              & FlashMatch::GetEstimatedX()          const { return this->fEstimatedX;          }
  const std::vector<double> & FlashMatch::GetHypoFlashSpec()       const { return this->fHypoFlashSpec;       }
  const std::vector<double> & FlashMatch::GetRecoFlashSpec()       const { return this->fRecoFlashSpec;       }
  const std::vector<double> & FlashMatch::GetMCFlashSpec()         const { return this->fMCFlashSpec;         }
  const std::vector<double> & FlashMatch::GetXFixedHypoFlashSpec() const { return this->fXFixedHypoFlashSpec; }
  const double              & FlashMatch::GetXFixedChi2()          const { return this->fXFixedChi2;          }
  const double              & FlashMatch::GetXFixedLl()            const { return this->fXFixedLl;            }


}


