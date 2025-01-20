// sgp4.h
#pragma once

#include "tle.h"
#include <cmath>

class SGP4 {
public:
    explicit SGP4(const TLE& tle);

    OrbitalState getPosition(const QDateTime& time) const;

private:
    struct PolarMotion {
        double xp; // в угловых секундах
        double yp; // в угловых секундах
    };
    static constexpr double XKE = 0.743669161E-1;
    static constexpr double AE = 1.0;
    static constexpr double XKMPER = 6378.137;    // Earth's radius in km
    static constexpr double CK2 = 5.413080E-4;
    static constexpr double CK4 = 0.62098875E-6;

    void initializeParameters();
    void calculateConstants();

    TLE tle_;

    // Orbital parameters
    double aodp_;
    double cosio_;
    double sinio_;
    double omgdot_;
    double xmdot_;
    double xlldot_;
    double xnodot_;
    double xnodp_;
    double calculateGMST(const QDateTime &time) const;
    Vector3 applyPolarMotion(const Vector3 &ecef, const PolarMotion &pm) const;
    void calculatePerturbations(double tsince, Vector3 &pos, Vector3 &vel) const;
    double solveKepler(double M, double e, double tolerance = 1e-12) const;
};
