#ifndef ZeroDark_h
#define ZeroDark_h

/** This is free and unencumbered software released into the public domain.

The authors of ISIS do not claim copyright on the contents of this file.
For more details about the LICENSE terms and the AUTHORS, you will
find files of those names at the top level of this repository. **/

/* SPDX-License-Identifier: CC0-1.0 */

#include <cmath>
#include <string>
#include <vector>
#include <tuple>

#include "IString.h"
#include "HiCalTypes.h"
#include "HiCalUtil.h"
#include "HiCalConf.h"
#include "Module.h"
#include "FileName.h"
#include "LoadCSV.h"
#include "LowPassFilter.h"
#include "Statistics.h"
#include "IException.h"

namespace Isis {

  /**
   * @brief Computes a complex dark subtraction component (ZeroDark module)
   *
   * This class computes the HiRISE dark correction component using a
   * combination of the B matrix, slope/intercept components and temperature
   * profiles.
   *
   * @ingroup Utility
   *
   * @author 2008-01-10 Kris Becker
   * @internal
   * @history 2008-06-13 Kris Becker - Added PrintOn method to produce more
   *          detailed data dump;  Added computation of statistics
   * @history 2010-04-16 Kris Becker - Implemented standardized access to CSV
   *          files.
   * @history 2010-10-28 Kris Becker Renamed parameters replacing "Zb" with
   *            "ZeroDark"
   * @history 2026-07-31 Kris J. Becker Implemented ADC "hot" algorithm for 
   *            new B matrices and temperature conditions
   *
   */
  class ZeroDark : public Module {

    public:
      static inline double ADCFpaTemperatureThresholdDefault = 31.0;
      static inline QString ADCSettingsDefault = "54";

      //  Constructors and Destructor
      ZeroDark() : Module("ZeroDark") { }
      ZeroDark(const HiCalConf &conf) : Module("ZeroDark") {
        init(conf);
      }

      /** Destructor */
      virtual ~ZeroDark() { }

      /**
       * @brief Return statistics for filtered - raw Buffer
       *
       * @return const Statistics&  Statistics class with all stats
       */
      const Statistics &Stats() const { return (_stats); }

    private:
      int _tdi;
      int _bin;
      double _fpa_average_temp;


      HiVector _BM;
      HiVector _slope;
      HiVector _intercept;
      HiVector _tempProf;
      HiVector _dc;

      double _refTemp;

      Statistics _stats;

      /**
       * @brief Determine ADC FPA temperature threshold to invoke hot calibration correction
       * 
       * This method will fetch the effective FPA temperature minimum value to
       * determine which dark current correction algorithm to apply for this
       * image. 
       * 
       * The label keyword ADCTimingSetting contains a two element array that
       * specifies the ADC settings used for the image. These values are concatenated
       * and serves as an identifier for the setting that pertains to the hical
       * application.
       * 
       * This method first searches in the ZeroDark module contained in the
       * hical CONFIG file for a specific keyword of the form "ADC{ADC}_FpaTemperatureThreshold".
       * If it exists, it is extracted and return as well as the keyword name.
       * If it does not exist, then it will retrieve the default ADC temperature
       * from the ZeroDark module configuration called "ADC_FpaTemperatureThreshold".
       * If this keyword does not exist, an internal default specifed by this
       * class default defined in ADCFpaTemperatureThresholdDefault = 31 will be
       * used.
       * 
       * Both the keyword and the value is returned to the caller.
       * 
       * @param profile DbProfile that contains the combined contents of the
       *                  of the ISIS cube label, the ZeroDark module config profile
       *                  and some internal keys initialized at this object startup
       * @return std::tuple<QString, double> Name of the CONFIG keyword found in
       *                  the CONFIG file and the FPA temperature value (or internal
       *                  default)
       */
      inline std::tuple<QString, double> getADCFpaTemperatureThreshold( const HiCalConf &conf, 
                                                                        const DbProfile &profile ) 
                                                                        const {
        QString adc_temp_key = conf.parser( "ADC{ADC}_FpaTemperatureThreshold", 
                                            HiCalConf::ValueList( { "ADC" } ),
                                            profile );
        // Check for existance of specialized ADC FPA temperature threshold
        if ( !profile.exists( adc_temp_key ) ) {
          adc_temp_key = "ADC_FpaTemperatureThreshold";
        }

        double adc_fpa_thresh = toDouble(ConfKey(profile, adc_temp_key, 
                                         toString(ADCFpaTemperatureThresholdDefault)));

        return ( std::make_tuple( adc_temp_key, adc_fpa_thresh ) );
      }

      /**
       * @brief Determine which ADC setting to use for the matrix file
       * 
       * This method determines which ADC setting to use to get a valid ADC B
       * matrix for the current settings. 
       * 
       * Using the ADCTimingSetting keyword in the labels, it first constructs
       * an ADC B matrix filename using those settings. If it finds the ADC B
       * matrix, this indicates it has a unique ADC B matrix for these settings.
       * 
       * If there is no ADC B matrix file for the current settings, then the
       * a default ADC B matrix is resolved as specified by the internal variable,
       * ADCSettingsDefault. If that should fail, return the initial DbProfile
       * which is likely to error out and report the missing configuration.
       * 
       * The act of finding the highest version of a file that contains "????"
       * will cause an exception. This method catches these and will return the
       * initial condition if nither is found. It will error with a better error,
       * hopefully.
       * 
       * @param conf HiCalConf label/model configuration data
       * @param profile Current hical Module configuration 
       * @return std::tuple<QString, DbProfile> Returns the namne of the resolved
       *           ADC B matrix if it exists and the DbProfile to use to formally
       *           load the ADC B matrix. A blank filename indicates an error.
       */
      inline std::tuple<QString, DbProfile> getADCMatrixFile( const HiCalConf &conf, 
                                                              const DbProfile &profile ) 
                                                              const {
        
        // Retrieve current ADC settings filename and check if the file exists.
        // If not fetch the current default name
        QString adc_matrix_file;
        try {
           adc_matrix_file = conf.getMatrixSource("ADC_DarkCurrent", profile );
          if ( FileName( adc_matrix_file ).fileExists() ) {
            return ( std::make_tuple( adc_matrix_file, profile ) );
          }
        }
        catch ( ... ) {
          // noop in this case
        }


        // There is no ADC B matrix found for this setting, so use the default
        try {
          DbProfile profile_adc = profile;
          profile_adc.replace( "ADC", ADCSettingsDefault );
          adc_matrix_file = conf.getMatrixSource("ADC_DarkCurrent", profile_adc );
          if ( FileName( adc_matrix_file ).fileExists() ) {
            return ( std::make_tuple( adc_matrix_file, profile_adc ) );
          }
        }
        catch ( ... ) {
          // noop here too
        }
        

        // Neither found, so return this condition
        return ( std::make_tuple( "ADC_DarkCurrent", profile ) );
      }


      /**
       * @brief Main init method that determines sets up algorithm processing
       * 
       * This init method computes common data and determines which algorithm to
       * process the ZeroDark calibration.
       * 
       * @param conf hical calibration data
       */
      void init(const HiCalConf &conf) {

        // Get/record common data
        _history.clear();
        DbProfile profile = conf.getMatrixProfile();
        _history.add("Profile["+ profile.Name()+"]");
        _tdi = ToInteger(profile("Tdi"));
        _bin = ToInteger(profile("Summing"));

      //  Set average tempuratures
        double fpa_py_temp = ToDouble(profile("FpaPositiveYTemperature"));
        double fpa_my_temp = ToDouble(profile("FpaNegativeYTemperature"));
        _fpa_average_temp = (fpa_py_temp+fpa_my_temp) / 2.0;
        _history.add("FpaAverageTemperature[" + ToString(_fpa_average_temp) + "]");      

        // Compute the ADC FPA temperature threshold 
        auto [ adc_fpa_temp_key, adc_fpa_temp_threshold ] = getADCFpaTemperatureThreshold(conf, profile );

        // Now check for the proper algorithm to apply
        if ( _fpa_average_temp > adc_fpa_temp_threshold ) {
          applyADCHotCalibration( conf, profile, adc_fpa_temp_key, adc_fpa_temp_threshold );
        }
        else {
          applyColdCalibration( conf, profile );
        }

      }

      /**
       * @brief Apply the original cold image calibration routine
       * 
       * This method applies the original cold ZeroDark algorithm whem the average
       * FPA temperature is less than the ADC_FpaTemperatureThreshold. 
       * 
       * @param conf    hical system calibration data
       * @param profile ZeroDark calibration data
       */
      inline void applyColdCalibration ( const HiCalConf &conf,
                                         const DbProfile &profile ) {
 
        int samples = ToInteger(profile("Samples"));        

        //  Get dark current (B) matrix, slope and intercept CSV files
        _BM = loadCsv("DarkCurrent", conf, profile, samples);
        _slope = loadCsv("DarkSlope", conf, profile, 256);
        _intercept = loadCsv("DarkIntercept", conf, profile, 256);

        // Get temperation normalization factor
        _refTemp = toDouble(ConfKey(profile, "FpaReferenceTemperature", toString(21.0)));
        _history.add("FpaReferenceTemperature["+toString(_refTemp)+"])");

        //  Smooth/filter if requested
        int width =  toInt(ConfKey(profile,"ZeroDarkFilterWidth",toString(3)));
        int iters =  toInt(ConfKey(profile,"ZeroDarkFilterIterations",toString(0)));
        LowPassFilter smooth(width, iters);
        _history.add("Smooth(Width["+ToString(width)+"],Iters["+ToString(iters)+"])");

        //  Filter the slope/intercept
        smooth.Process(_slope);
        _slope = smooth.ref();

        smooth.Process(_intercept);
        _intercept = smooth.ref();

        HiVector t_prof(_slope.dim());
        for (int i = 0 ; i < _slope.dim() ; i++) {
          t_prof[i] = _intercept[i] + _slope[i] * _fpa_average_temp;
        }

        _tempProf = rebin(t_prof, samples);
        _history.add("Rebin(T_Profile," + ToString(t_prof.dim()) + "," +
                     ToString(samples) +")");

        _dc = HiVector(samples);
        double linetime = ToDouble(profile("ScanExposureDuration"));
        double scale = linetime * 1.0E-6 * (_bin*_bin) *
                       (20.0*103.0/89.0 + _tdi);
        double baseT = HiTempEqn(_refTemp);
        for (int j = 0 ; j < samples ; j++) {
          _dc[j] = _BM[j] * scale * HiTempEqn(_tempProf[j]) / baseT;
        }
        _history.add("HiTempEqn(FpaReferenceTemperature)[" + ToString(baseT) + "]");
        _history.add("scale[" + ToString(scale) + "]");

        //  Filter it yet again
        smooth.Process(_dc);
        _data = smooth.ref();

        //  Compute statistics and record to history
        _stats.Reset();
        for ( int i = 0 ; i < _data.dim() ; i++ ) {
          _stats.AddData(_data[i]);
        }
        _history.add("Statistics(Average["+ToString(_stats.Average())+
                     "],StdDev["+ToString(_stats.StandardDeviation())+"])");
        return;
      }

      /**
       * @brief Apply hot ZeroDark calibration algorithm
       * 
       * This method applies a special modified dark calibration algorithm designed
       * to process images taken during warmer camera operations. These data
       * require differnt B matrix that are specially derived for these conditions.
       * The matrices contain different data that are not consistent with the
       * original B matrix so requires a different implementation.
       * 
       * @param conf     hical calibration data
       * @param profile  ZeroDark configuration data
       * @param adc_reftemp_key ADC reference temperature key
       * @param adc_reftemp     ADC reference temperature that may be specially
       *                          set for the ADC settings
       */
      inline void applyADCHotCalibration( const HiCalConf &conf,
                                          const DbProfile &profile,
                                          const QString &adc_reftemp_key,
                                          const double &adc_reftemp ) {

        _history.add("ADC_FpaTemperatureThresholdKey["+adc_reftemp_key+"])");
        _history.add("ADC_FpaTemperatureThreshold["+ToString(adc_reftemp)+"])");

        int samples = ToInteger(profile("Samples"));        

        //  Get dark current (B) matrix
        auto [ adc_matrix_file, adc_profile ] = getADCMatrixFile( conf, profile );
        _BM = loadCsv("ADC_DarkCurrent", conf, adc_profile, samples);

        // These two data elements are unused for this algorithm 
        _slope = HiVector(samples, 0.0 );
        _intercept = HiVector(samples, 0.0 );
        
        // GetADX temperature normalization factor
        _refTemp = toDouble(ConfKey(profile, "ADC_FpaReferenceTemperature", toString(37.0)));
        _history.add("ADC_FpaReferenceTemperature["+toString(_refTemp)+"])");

        //  Smooth/filter if requested
        int width =  toInt(ConfKey(profile,"ZeroDarkFilterWidth",toString(3)));
        int iters =  toInt(ConfKey(profile,"ZeroDarkFilterIterations",toString(0)));
        LowPassFilter smooth(width, iters);
        _history.add("Smooth(Width["+ToString(width)+"],Iters["+ToString(iters)+"])");

        _dc = HiVector(samples);
        double linetime = ToDouble(profile("ScanExposureDuration"));
        double scale = linetime * 1.0E-6 * _tdi;
        double baseT = HiTempEqn(_refTemp);
        double hiTemp = HiTempEqn(_fpa_average_temp);
        _tempProf = HiVector(samples, hiTemp);

        _history.add("HiTempEqn(ADC_FpaReferenceTemperature)[" + ToString(baseT) + "]");
        _history.add("HiTempEqn(ADC_FpaAverageTemperature)[" + ToString(hiTemp) + "]");
        _history.add("scale[" + ToString(scale) + "]");
        
        for (int j = 0 ; j < samples ; j++) {
          _dc[j] = _BM[j] * scale * hiTemp / baseT;
        }

        //  Filter it yet again
        smooth.Process(_dc);
        _data = smooth.ref();

        //  Compute statistics and record to history
        _stats.Reset();
        for ( int i = 0 ; i < _data.dim() ; i++ ) {
          _stats.AddData(_data[i]);
        }
        _history.add("Statistics(Average["+ToString(_stats.Average())+
                     "],StdDev["+ToString(_stats.StandardDeviation())+"])");
        return;
      }

      /** Virtualized data dump method */
      virtual void printOn(std::ostream &o) const {
        o << "#  History = " << _history << std::endl;
        //  Write out the header
        o << std::setw(_fmtWidth)   << "DarkMatrix"
          << std::setw(_fmtWidth+1) << "TempNorm"
          << std::setw(_fmtWidth+1) << "ZeroDark\n";

        for (int i = 0 ; i < _data.dim() ; i++) {
          o << formatDbl(_BM[i]) << " "
            << formatDbl(_tempProf[i]) << " "
            << formatDbl(_data[i]) << std::endl;
        }
        return;
      }

  };

}     // namespace Isis
#endif
