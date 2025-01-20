#ifndef EARTH_MODEL_H
#define EARTH_MODEL_H
#pragma once

namespace EarthModel {
// WGS84 Constants
constexpr double SEMI_MAJOR_AXIS = 6378.137;          // Earth's semi-major axis (km)
constexpr double FLATTENING = 1.0/298.257223563;      // WGS84 flattening
constexpr double SEMI_MINOR_AXIS = SEMI_MAJOR_AXIS * (1.0 - FLATTENING);
constexpr double ECCENTRICITY_SQ = FLATTENING * (2.0 - FLATTENING);

// Earth's Gravitational Model
constexpr double GM = 398600.4418;                    // Earth's gravitational constant (km³/s²)
constexpr double J2 = 1.082629E-3;                    // J2 perturbation coefficient
constexpr double J3 = -2.53881E-6;                    // J3 perturbation coefficient
constexpr double J4 = -1.65597E-6;                    // J4 perturbation coefficient

// Earth Rotation
constexpr double OMEGA_E = 7.292115E-5;               // Earth's rotation rate (rad/sec)

struct PolarMotion {
    double xp;    // X pole position (arcseconds)
    double yp;    // Y pole position (arcseconds)
    double dut1;  // UT1-UTC difference (seconds)

    // Default modern values (should be updated with IERS data)
    PolarMotion() : xp(0.162), yp(0.358), dut1(0.1789) {}
};
}
#endif // EARTH_MODEL_H
