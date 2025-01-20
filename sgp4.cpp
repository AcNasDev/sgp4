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

    // Initialize aodp_ and xnodp_
    aodp_ = tle_.a;
    xnodp_ = tle_.no;
}

void SGP4::calculateConstants() {
    // Convert degrees to radians
    double rad_per_deg = M_PI / 180.0;

    // Calculate orientation angles in radians
    double xincl = tle_.inclination * rad_per_deg;
    cosio_ = cos(xincl);
    sinio_ = sin(xincl);

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

    xmdot_ = tle_.no * (1.0 + temp1);
    omgdot_ = temp2 * tle_.no;
    xnodot_ = temp3 * tle_.no;
}

OrbitalState SGP4::getPosition(const QDateTime& time) const {
    // Calculate time since epoch in minutes
    int year = tle_.epochYear;
    if (year < 57) year += 2000;
    else year += 1900;

    QDateTime epoch = QDateTime(QDate(year, 1, 1), QTime(0, 0), Qt::UTC);
    epoch = epoch.addSecs(static_cast<qint64>((tle_.epochDay - 1.0) * 86400.0));

    double tsince = epoch.secsTo(time) / 60.0;  // Convert to minutes

    // Implementation of the SGP4 propagator
    double xmo = tle_.meanAnomaly * M_PI / 180.0 + xmdot_ * tsince;
    double omegao = tle_.argumentPerigee * M_PI / 180.0 + omgdot_ * tsince;
    double xno = tle_.rightAscension * M_PI / 180.0 + xnodot_ * tsince;

    // Kepler's equation
    double e = tle_.eccentricity;
    double a = tle_.a * XKMPER;  // Convert to km

    double M = fmod(xmo, 2.0 * M_PI);
    double E = M;

    // Solve Kepler's equation using Newton-Raphson method
    for(int i = 0; i < 10; i++) {
        double E_new = E - (E - e * sin(E) - M) / (1.0 - e * cos(E));
        if(std::abs(E_new - E) < 1e-12) {
            E = E_new;
            break;
        }
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

    // Calculate velocity (simplified)
    double xn = xnodp_;
    double rdot = -xn * a * e * sin(E) / sqrt(1.0 - e * e);
    double rfdot = xn * a * sqrt(1.0 - e * e) * cos(E) / (1.0 - e * cos(E));

    Vector3 vel;
    vel.x = rdot * cos(nu) - r * rfdot * sin(nu);
    vel.y = rdot * sin(nu) + r * rfdot * cos(nu);
    vel.z = 0.0;  // Simplified - needs proper implementation

    return OrbitalState{pos, vel};
}
