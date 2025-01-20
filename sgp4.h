// sgp4.h
#pragma once

#include "tle.h"
#include <cmath>

class SGP4 {
public:
    explicit SGP4(const TLE& tle);

    OrbitalState getPosition(const QDateTime& time) const;

private:
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
};
