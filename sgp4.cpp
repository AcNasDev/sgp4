// sgp4.cpp
#include "sgp4.h"
#include <QDebug>

SGP4::SGP4(const TLE& tle) : tle_(tle) {
    initializeParameters();
    calculateConstants();
}

void SGP4::initializeParameters() {
    // Constants
    const double MINUTES_PER_DAY = 1440.0;
    const double TWO_PI = 2.0 * M_PI;
    const double XKE = 0.743669161E-1;
    const double GM = 398600.4418; // Earth's gravitational constant (km³/s²)

    // Convert mean motion from revs/day to rad/min
    tle_.no = tle_.meanMotion * TWO_PI / MINUTES_PER_DAY;

    // Calculate semi-major axis
    double n = tle_.meanMotion * TWO_PI / 86400.0; // Convert to rad/sec
    double a3 = GM / (n * n);
    double a = std::pow(a3, 1.0/3.0); // Semi-major axis in km
    tle_.a = a / XKMPER; // Convert to Earth radii

    // Compute perigee and apogee (in Earth radii)
    tle_.alta = tle_.a * (1.0 + tle_.eccentricity) - 1.0;
    tle_.altp = tle_.a * (1.0 - tle_.eccentricity) - 1.0;

    // Initialize additional parameters
    aodp_ = tle_.a;
    xnodp_ = tle_.no;
}

void SGP4::calculateConstants() {
    // Convert degrees to radians
    double rad_per_deg = M_PI / 180.0;
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
    if (year < 57) {
        year += 2000;
    } else {
        year += 1900;
    }

    QDateTime epoch = QDateTime(QDate(year, 1, 1), QTime(0, 0), Qt::UTC);
    epoch = epoch.addSecs(static_cast<qint64>((tle_.epochDay - 1.0) * 86400.0));

    double tsince = epoch.secsTo(time) / 60.0;  // Convert to minutes

    // Compute angles in radians
    double xmo = fmod(tle_.meanAnomaly * M_PI / 180.0 + xmdot_ * tsince, 2.0 * M_PI);
    double omegao = tle_.argumentPerigee * M_PI / 180.0 + omgdot_ * tsince;
    double xno = tle_.rightAscension * M_PI / 180.0 + xnodot_ * tsince;

    // Solve Kepler's equation
    double e = tle_.eccentricity;
    double E = xmo;

    for (int i = 0; i < 10; i++) {
        double E_new = E - (E - e * sin(E) - xmo) / (1.0 - e * cos(E));
        if (std::abs(E_new - E) < 1e-12) {
            E = E_new;
            break;
        }
        E = E_new;
    }

    // Calculate true anomaly
    double sinE = sin(E);
    double cosE = cos(E);
    double sinv = sqrt(1.0 - e * e) * sinE / (1.0 - e * cosE);
    double cosv = (cosE - e) / (1.0 - e * cosE);
    double nu = atan2(sinv, cosv);

    // Calculate position in orbital plane
    double r = tle_.a * XKMPER * (1.0 - e * cosE);  // Convert to km
    double u = omegao + nu;

    double sin_u = sin(u);
    double cos_u = cos(u);
    double sin_i = sinio_;
    double cos_i = cosio_;
    double sin_node = sin(xno);
    double cos_node = cos(xno);

    // Convert position to geocentric coordinates
    Vector3 pos;
    pos.x = r * (cos_node * cos_u - sin_node * sin_u * cos_i);
    pos.y = r * (sin_node * cos_u + cos_node * sin_u * cos_i);
    pos.z = r * sin_u * sin_i;

    // Calculate velocity (simplified)
    double xn = xnodp_;
    double rdot = xn * tle_.a * XKMPER * e * sinE / sqrt(1.0 - e * e);
    double rfdot = xn * tle_.a * XKMPER * sqrt(1.0 - e * e) / (1.0 - e * cosE);

    Vector3 vel;
    vel.x = rdot * (cos_node * cos_u - sin_node * sin_u * cos_i) -
            rfdot * (cos_node * sin_u + sin_node * cos_u * cos_i);
    vel.y = rdot * (sin_node * cos_u + cos_node * sin_u * cos_i) -
            rfdot * (sin_node * sin_u - cos_node * cos_u * cos_i);
    vel.z = rdot * sin_u * sin_i + rfdot * cos_u * sin_i;

    return OrbitalState{pos, vel};
}
