/**
 * \class ubana::FlashMatch
 *
 * \ingroup UBXSec
 *
 * \brief Data product to store a flash matching results
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
 * Created on: Friday, February 03, 2017 at 16:34:34
 *
 */

#ifndef FlashMatch_h
#define FlashMatch_h

#include <vector>

namespace ubana {

  class FlashMatch {

  public:

    FlashMatch();
    FlashMatch(double score);
    virtual ~FlashMatch();

    // Setter methods
    void SetScore(double);
    void SetEstimatedX(double);
    void SetTPCX(double);
    void SetT0(double);
    void SetHypoFlashSpec(std::vector<double>);
    void SetRecoFlashSpec(std::vector<double>);
    void SetMCFlashSpec(std::vector<double>);
    void SetXFixedHypoFlashSpec(std::vector<double>);
    void SetXFixedChi2(double);
    void SetXFixedLl(double);

    // Getter methods
    const double &              GetScore()                const;
    const double &              GetEstimatedX()           const;
    const double &              GetTPCX()                 const;
    const double &              GetT0()                   const;
    const std::vector<double> & GetHypoFlashSpec()        const;
    const std::vector<double> & GetRecoFlashSpec()        const;
    const std::vector<double> & GetMCFlashSpec()          const;
    const std::vector<double> & GetXFixedHypoFlashSpec()  const;
    const double &              GetXFixedChi2()           const;
    const double &              GetXFixedLl()             const;

  private:

    double fScore;
    double fTPCX;
    double fEstimatedX;
    double fT0;
    std::vector<double> fHypoFlashSpec;
    std::vector<double> fRecoFlashSpec;
    std::vector<double> fMCFlashSpec;
    std::vector<double> fXFixedHypoFlashSpec;
    double fXFixedChi2;
    double fXFixedLl;

 };
}

#endif /* FlashMatch_h */
