// sgp4.cpp
#include "sgp4.h"
#include <QDebug>

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

    // Calculate secular rates with perturbations
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

    // Compute mean motion
    double xn = xnodp_;

    // Update for secular gravity and atmospheric drag
    double xmdf = tle_.meanAnomaly * M_PI / 180.0 + xmdot_ * tsince;
    double omgadf = tle_.argumentPerigee * M_PI / 180.0 + omgdot_ * tsince;
    double xnode = tle_.rightAscension * M_PI / 180.0 + xnodot_ * tsince;

    // Solve Kepler's equation
    double e = tle_.eccentricity;
    double E = xmdf;
    double M = xmdf;

    for (int i = 0; i < 10; i++) {
        double E_new = M + e * sin(E);
        if (std::abs(E_new - E) < 1e-12) {
            E = E_new;
            break;
        }
        E = E_new;
    }

    // Calculate true anomaly
    double nu = 2.0 * atan2(sqrt(1.0 + e) * sin(E/2.0), sqrt(1.0 - e) * cos(E/2.0));

    // Calculate position in orbital plane
    double r = aodp_ * XKMPER * (1.0 - e * cos(E));
    double u = omgadf + nu;

    // Calculate position in orbital plane
    double x = r * cos(u);
    double y = r * sin(u);
    double z = 0.0;

    // Rotate to geocentric coordinates
    double sinNode = sin(xnode);
    double cosNode = cos(xnode);
    double sinI = sinio_;
    double cosI = cosio_;
    double sinU = sin(u);
    double cosU = cos(u);

    Vector3 pos;
    pos.x = x * (cosNode * cosU - sinNode * cosI * sinU) -
            y * (cosNode * sinU + sinNode * cosI * cosU);
    pos.y = x * (sinNode * cosU + cosNode * cosI * sinU) +
            y * (cosNode * cosI * cosU - sinNode * sinU);
    pos.z = x * sinI * sinU + y * sinI * cosU;

    // Calculate velocity (simplified)
    double xdot = -xn * r * sin(nu);
    double ydot = xn * r * (e + cos(nu));

    Vector3 vel;
    vel.x = xdot * (cosNode * cosU - sinNode * cosI * sinU) -
            ydot * (cosNode * sinU + sinNode * cosI * cosU);
    vel.y = xdot * (sinNode * cosU + cosNode * cosI * sinU) +
            ydot * (cosNode * cosI * cosU - sinNode * sinU);
    vel.z = xdot * sinI * sinU + ydot * sinI * cosU;

    return OrbitalState{pos, vel};
}
