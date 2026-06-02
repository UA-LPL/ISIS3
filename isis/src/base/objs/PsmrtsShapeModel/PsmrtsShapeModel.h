#ifndef PsmrtsShapeModel_h
#define PsmrtsShapeModel_h
/** This is free and unencumbered software released into the public domain.
The authors of ISIS do not claim copyright on the contents of this file.
For more details about the LICENSE terms and the AUTHORS, you will
find files of those names at the top level of this repository. **/

/* SPDX-License-Identifier: CC0-1.0 */
#include "ShapeModel.h"

#include <vector>

#include <QString>

#include "Latitude.h"
#include "Longitude.h"
#include "Pvl.h"
#include "SurfacePoint.h"

#include <psmrts/tracers/PsmrtsTracerSystem.hpp>
#include "IsisPsmrtsUtilities.hpp"


namespace Isis {
  class Target;

  /**
   * Implementation of the PSMRTS shape model class
   *
   * @author 2026-03-04 Kris J. Becker
   *
   * @internal
   *   @history 2026-03-04 - Kris J. Becker - Original Version
 
   */
  class PsmrtsShapeModel : public ShapeModel {
    public:
    /** Set default surface point tolerance to millimeter precision */
      static inline const double DefaultDistanceTolerance = 1.0e-6;

      // Constructors
      PsmrtsShapeModel();
      PsmrtsShapeModel(Target *target, Pvl &pvl);
      PsmrtsShapeModel(const psmrts::PsmrtsTracerSystem &tracer_s );

      // Destructor
      ~PsmrtsShapeModel();

      // Psmrts shape model creator methods */
      static PsmrtsShapeModel *create(Target *target, Pvl &pvl, 
                                      const bool throw_errors = true);

      // Intersect the shape model
      bool intersectSurface(std::vector<double> observerPos,
                            std::vector<double> lookDirection);
      virtual bool intersectSurface(const Latitude &lat, const Longitude &lon,
                                    const std::vector<double> &observerPos,
                                    const bool &checkOcclusion = true);
      virtual bool intersectSurface(const SurfacePoint &surfpt, 
                                    const std::vector<double> &observerPos,
                                    const bool &checkOcclusion = true);


      // Calculate the surface normal of the current intersection point
      void calculateDefaultlNormal();
      void calculateLocalNormal(QVector<double *> cornerNeighborPoints); // use default normal
      // void calculateSurfaceNormal();
      void setLocalNormalFromIntercept();


      // Determine if the internal intercept is occluded from the observer/lookdir
      virtual bool isVisibleFrom(const std::vector<double> observerPos,
                                 const std::vector<double> lookDirection);

      virtual void clearSurfacePoint();

      // CSMCamera uses all of these so be sure they are implemented!!
      // Calculate the emission angle of the current intersection point
      virtual double emissionAngle(const std::vector<double> & sB);

      // Calculate the incidence angle of the current intersection point
      virtual double incidenceAngle(const std::vector<double> &uB);

      Distance localRadius(const Latitude &lat, const Longitude &lon);
    
      bool isDEM() const;

      virtual void setSurfacePoint(const SurfacePoint &surfacePoint);

      inline const psmrts::PsmrtsTracerSystem &tracer_system() const {
        return ( m_tracer_s );
      }

      inline const psmrts::PRQRayTrace &get_shape_trace() const {
        return ( m_shape_ray_t );
      }

      inline const psmrts::PRQRayTrace &get_ellipsoid_trace() const {
        return ( m_ellipsoid_ray_t );
      }

      inline double get_tolerance() const {
        return ( m_tolerance );
      }

      inline void set_tolerance( const double tolerance = DefaultDistanceTolerance ) {
        m_tolerance = tolerance;
        return;
      }

      /**
       * @brief Compute the surface intersept point given a lat/lon coordinate
       * 
       * This method computes the surface point on the shape DEM given a
       * latitude/longitude mapping coordinate. This is useful for
       * orthorectified mapping operations to determined precision surface
       * intercept points and determining that radius at that coordinate. It
       * does not require an observer position as it is computed as the vector
       * from the body origin through the latitude/longitude coordinate extended
       * beyond the surface 1.5 times the maximum radius of the shape model. The
       * look direction is the negative of the this observer position. The
       * resulting trace of this configuration is returned as the surface
       * intercept point of this coordinate.
       * 
       * @param lat   Latitude of the mapping coordinate
       * @param lon   Longitude of the mapping coordinate
       * @return psmrts::PRQRayTrace Result of the shape model trace 
       */
      template <typename T>
        psmrts::PRQRayTrace getSurfaceLatLonRadius( const Latitude &lat, 
                                                    const Longitude &lon,
                                                    const T &tracer ) 
                                                    const {
          // Compute the observer position ensuring it is above the shape                                          
          double max_r = tracer.maximum_radius();
          Eigen::Vector3d llr( lon.degrees(), lat.degrees(), max_r * 1.5 );
          Eigen::Vector3d observer = psmrts::lonlatrad_to_xyz_d( llr );

          Eigen::Vector3d lookdir = -observer;  // Just negate the observer
          psmrts::PRQRayTrace ray_t( observer, lookdir );

          // Trace it to the surface for intersect potential
          (void) tracer.process( ray_t );
          return ( ray_t );
        }

      bool load_pvl_config( const QString &pvlconf_f, PvlFlatMap &flat_p ) const;
      bool load_shape_list( const QString &shapelist_f, PvlFlatMap &flat_p ) const;

    private:
      // Disallow copying because ShapeModel is not copyable
      Q_DISABLE_COPY(PsmrtsShapeModel)

      PvlFlatMap                 m_parameters;
      psmrts::PRQRayTrace        m_shape_ray_t;
      psmrts::PRQRayTrace        m_ellipsoid_ray_t;
      psmrts::PsmrtsTracerSystem m_tracer_s;
      double                     m_tolerance;

      inline bool psmrtsUpdateIsisLabel( Pvl &pvl, const PvlFlatMap &psmrts_data ) const; 
 
      inline void reset_all_rays( ) {
        m_shape_ray_t.reset();
        m_ellipsoid_ray_t.reset();
      }

      /**
       * @brief Update the internal state with the result of a ray trace 
       * 
       * This method updates the internal state of this shape model object with
       * the results of the ray trace object. If the trace has a hit the surface
       * point is updated with the surface point and normal data.
       * 
       * If the is no hit or an error occured, then the surface point is cleared.
       * 
       * The state of the result is returned to the caller.
       * 
       * @param ray Ray trace object to set 
       */
      inline bool updateTraceState( const psmrts::PRQRayTrace &ray,
                                    const psmrts::PRQRayTrace &ray_e = psmrts::PRQRayTrace( ) ) {
        if ( ray.hasHit() ) {

          setHasIntersection( ray.hasHit() ); 
          SurfacePoint point;
          point.FromNaifArray( ray.trace().xyz().data() );
          ShapeModel::setSurfacePoint( point );

          // Got the local normal so set it here
          Eigen::Vector3d normal_l = ray.trace().normal();
          setLocalNormal( normal_l[0], normal_l[1], normal_l[2] );

          // Set the ellipsoid normal as well
          if ( ray_e.hasHit() ) {
            Eigen::Vector3d normal_e = ray_e.trace().normal();
            setNormal( normal_e[0], normal_e[1], normal_e[2] );
          }
          else {
            setNormal( normal_l[0], normal_l[1], normal_l[2] );
          }
        }
        else {
          clearSurfacePoint();
        }

        return ( ray.hasHit() );
      }

  };
}

#endif
