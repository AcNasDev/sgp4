// sgp4.cpp
#include "sgp4.h"

SGP4::SGP4(const TLE& tle) : tle_(tle) {
    initializeParameters();
    calculateConstants();
}

void SGP4::initializeParameters() {
    // Convert TLE mean motion to internal units (rad/min)
    tle_.no = tle_.meanMotion * 2.0 * M_PI / 1440.0;

    // Semi-major axis (earth radii)
    tle_.a = pow(XKE / tle_.no, 2.0/3.0);

    // Compute perigee and apogee
    double e = tle_.eccentricity;
    tle_.alta = tle_.a * (1.0 + e) - 1.0;
    tle_.altp = tle_.a * (1.0 - e) - 1.0;
}

void SGP4::calculateConstants() {
    // Convert degrees to radians
    double rad_per_deg = M_PI / 180.0;
    double temp = tle_.inclination * rad_per_deg;
    cosio_ = cos(temp);
    sinio_ = sin(temp);

    // Calculate Earth's gravity field
    double theta2 = cosio_ * cosio_;
    double x3thm1 = 3.0 * theta2 - 1.0;
    double eosq = tle_.eccentricity * tle_.eccentricity;
    double betao2 = 1.0 - eosq;
    double betao = sqrt(betao2);

    // Calculate secular rates
    double temp1 = 1.5 * CK2 * x3thm1 / (betao * betao2);
    double temp2 = -0.5 * temp1;
    double temp3 = -0.5 * CK2 * x3thm1;

    xmdot_ = tle_.no + temp1 * tle_.no;
    omgdot_ = temp2 * tle_.no;
    xnodot_ = temp3 * tle_.no;
}

OrbitalState SGP4::getPosition(const QDateTime& time) const {
    // Calculate time since epoch in minutes
    QDateTime epoch = QDateTime::fromString(
        QString("%1%2").arg(tle_.epochYear).arg(tle_.epochDay, 0, 'f', 8),
        "yy.dddddddd");

    double tsince = epoch.msecsTo(time) / (1000.0 * 60.0);  // Convert to minutes

    // Implementation of the SGP4 propagator
    // This is a simplified version - full implementation would be much longer

    double xmo = tle_.meanAnomaly * M_PI / 180.0 + xmdot_ * tsince;
    double omegao = tle_.argumentPerigee * M_PI / 180.0 + omgdot_ * tsince;
    double xno = tle_.rightAscension * M_PI / 180.0 + xnodot_ * tsince;

    // Kepler's equation (simplified)
    double e = tle_.eccentricity;
    double a = tle_.a * XKMPER;  // Convert to km

    double M = xmo;
    double E = M;

    // Iterate to solve Kepler's equation
    for(int i = 0; i < 10; i++) {
        double E_new = M + e * sin(E);
        if(std::abs(E_new - E) < 1e-12) break;
        E = E_new;
    }

    // Calculate true anomaly
    double nu = 2.0 * atan(sqrt((1.0 + e)/(1.0 - e)) * tan(E/2.0));

    // Calculate position in orbital plane
    double r = a * (1.0 - e * cos(E));
    double x = r * cos(nu);
    double y = r * sin(nu);

    // Rotate to geocentric coordinates
    double xh = x * cos(omegao) - y * sin(omegao);
    double yh = x * sin(omegao) + y * cos(omegao);

    Vector3 pos;
    pos.x = xh * cos(xno) - yh * cosio_ * sin(xno);
    pos.y = xh * sin(xno) + yh * cosio_ * cos(xno);
    pos.z = yh * sinio_;

    // Velocity calculations would go here
    Vector3 vel;  // Simplified - actual implementation needed

    return OrbitalState{pos, vel};
}
