#ifndef COORDINATE_TRANSFORM_H
#define COORDINATE_TRANSFORM_H
#include "tle.h"
#pragma once
#include "earth_model.h"
#include <cmath>

class CoordinateTransform {
public:
    struct GeodeticCoord {
        double latitude;   // degrees
        double longitude;  // degrees
        double altitude;   // kilometers
    };

    static GeodeticCoord ecefToGeodetic(const Vector3& ecef) {
        const double& a = EarthModel::SEMI_MAJOR_AXIS;
        const double& b = EarthModel::SEMI_MINOR_AXIS;
        const double e2 = EarthModel::ECCENTRICITY_SQ;

        const double p = std::sqrt(ecef.x*ecef.x + ecef.y*ecef.y);
        const double theta = std::atan2(ecef.z*a, p*b);

        const double sin_theta = std::sin(theta);
        const double cos_theta = std::cos(theta);

        const double lat = std::atan2(
            ecef.z + e2*(a/b)*a*std::pow(sin_theta, 3),
            p - e2*a*std::pow(cos_theta, 3)
            );

        const double lon = std::atan2(ecef.y, ecef.x);

        const double sin_lat = std::sin(lat);
        const double N = a / std::sqrt(1.0 - e2*sin_lat*sin_lat);
        const double alt = p/std::cos(lat) - N;

        constexpr double RAD_TO_DEG = 180.0/M_PI;
        return {
            lat * RAD_TO_DEG,
            lon * RAD_TO_DEG,
            alt
        };
    }

    static Vector3 applyPolarMotion(const Vector3& ecef, const EarthModel::PolarMotion& pm) {
        const double xp_rad = pm.xp * (M_PI/(180.0*3600.0));
        const double yp_rad = pm.yp * (M_PI/(180.0*3600.0));

        Vector3 result;
        result.x = ecef.x - yp_rad*ecef.z;
        result.y = ecef.y + xp_rad*ecef.z;
        result.z = ecef.z + yp_rad*ecef.x - xp_rad*ecef.y;

        return result;
    }
};
#endif // COORDINATE_TRANSFORM_H
