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
    qDebug() << "\n=== SGP4 Position Calculation Debug ===";
    qDebug() << "Input UTC time:" << time.toString("yyyy-MM-dd HH:mm:ss");

    // Calculate epoch time
    int year = tle_.epochYear;
    if (year < 57) year += 2000;
    else year += 1900;

    QDateTime epoch = QDateTime(QDate(year, 1, 1), QTime(0, 0), Qt::UTC);
    epoch = epoch.addSecs(static_cast<qint64>((tle_.epochDay - 1.0) * 86400.0));

    qDebug() << "Epoch time:" << epoch.toString("yyyy-MM-dd HH:mm:ss");

    // Calculate time since epoch
    double tsince = epoch.secsTo(time) / 60.0;  // Convert to minutes
    qDebug() << "Time since epoch (minutes):" << tsince;

    // Compute angles
    double xmo = fmod(tle_.meanAnomaly * M_PI / 180.0 + xmdot_ * tsince, 2.0 * M_PI);
    double omegao = tle_.argumentPerigee * M_PI / 180.0 + omgdot_ * tsince;
    double xno = tle_.rightAscension * M_PI / 180.0 + xnodot_ * tsince;

    qDebug() << "\n--- Orbital Elements at Time t ---";
    qDebug() << "Mean Anomaly (rad):" << xmo;
    qDebug() << "Argument of Perigee (rad):" << omegao;
    qDebug() << "Right Ascension (rad):" << xno;

    // Solve Kepler's equation
    double e = tle_.eccentricity;
    double E = xmo;
    int iter = 0;

    qDebug() << "\n--- Kepler's Equation Solution ---";
    for (int i = 0; i < 10; i++) {
        double E_new = E - (E - e * sin(E) - xmo) / (1.0 - e * cos(E));
        double delta = std::abs(E_new - E);
        E = E_new;
        iter = i + 1;
        qDebug() << "Iteration" << iter << "- E:" << E << "Delta:" << delta;
        if (delta < 1e-12) break;
    }

    // Calculate true anomaly
    double sinE = sin(E);
    double cosE = cos(E);
    double sinv = sqrt(1.0 - e * e) * sinE / (1.0 - e * cosE);
    double cosv = (cosE - e) / (1.0 - e * cosE);
    double nu = atan2(sinv, cosv);

    qDebug() << "\n--- Position Calculation ---";
    qDebug() << "True Anomaly (rad):" << nu;

    // Calculate position in orbital plane
    double r = tle_.a * XKMPER * (1.0 - e * cosE);
    double u = omegao + nu;

    qDebug() << "Radius vector (km):" << r;
    qDebug() << "Argument of latitude (rad):" << u;

    // Convert to geocentric coordinates
    double sin_u = sin(u);
    double cos_u = cos(u);
    double sin_i = sinio_;
    double cos_i = cosio_;
    double sin_node = sin(xno);
    double cos_node = cos(xno);

    Vector3 pos;
    pos.x = r * (cos_node * cos_u - sin_node * sin_u * cos_i);
    pos.y = r * (sin_node * cos_u + cos_node * sin_u * cos_i);
    pos.z = r * sin_u * sin_i;

    qDebug() << "\n--- Final Position (km) ---";
    qDebug() << "X:" << pos.x;
    qDebug() << "Y:" << pos.y;
    qDebug() << "Z:" << pos.z;

    // Calculate velocity
    double xn = xnodp_;
    xn = xn / 60.0;  // Добавляем эту строку!
    double rdot = xn * tle_.a * XKMPER * e * sinE / sqrt(1.0 - e * e);
    double rfdot = xn * tle_.a * XKMPER * sqrt(1.0 - e * e) / (1.0 - e * cosE);

    qDebug() << "\n--- Velocity Components ---";
    qDebug() << "Radial velocity (km/s):" << rdot;
    qDebug() << "Cross-track velocity (km/s):" << rfdot;

    Vector3 vel;
    vel.x = rdot * (cos_node * cos_u - sin_node * sin_u * cos_i) -
            rfdot * (cos_node * sin_u + sin_node * cos_u * cos_i);
    vel.y = rdot * (sin_node * cos_u + cos_node * sin_u * cos_i) -
            rfdot * (sin_node * sin_u - cos_node * cos_u * cos_i);
    vel.z = rdot * sin_u * sin_i + rfdot * cos_u * sin_i;

    qDebug() << "\n--- Final Velocity (km/s) ---";
    qDebug() << "VX:" << vel.x;
    qDebug() << "VY:" << vel.y;
    qDebug() << "VZ:" << vel.z;

    qDebug() << "=== End SGP4 Position Calculation ===\n";

    QDateTime j2000(QDate(2000, 1, 1), QTime(12, 0), Qt::UTC);
    double days_since_j2000 = j2000.daysTo(time) + time.time().msecsSinceStartOfDay() / (1000.0 * 86400.0);

    // GMST в радианах
    double gmst = fmod(280.4606 + 360.9856473 * days_since_j2000, 360.0) * M_PI / 180.0;

    // На этом этапе pos и vel находятся в системе ECI
    // Преобразуем их в ECEF, учитывая вращение Земли

    // Преобразование позиции из ECI в ECEF
    Vector3 pos_ecef;
    pos_ecef.x = pos.x * cos(gmst) + pos.y * sin(gmst);
    pos_ecef.y = -pos.x * sin(gmst) + pos.y * cos(gmst);
    pos_ecef.z = pos.z;

    // Угловая скорость вращения Земли (рад/с)
    const double earth_rotation_rate = 7.2921150e-5;

    // Преобразование скорости из ECI в ECEF
    Vector3 vel_ecef;
    vel_ecef.x = vel.x * cos(gmst) + vel.y * sin(gmst) - earth_rotation_rate * pos_ecef.y;
    vel_ecef.y = -vel.x * sin(gmst) + vel.y * cos(gmst) + earth_rotation_rate * pos_ecef.x;
    vel_ecef.z = vel.z;

    qDebug() << "ECI to ECEF conversion:";
    qDebug() << "GMST (degrees):" << gmst * 180.0 / M_PI;
    qDebug() << "ECI position (km):" << pos.x << pos.y << pos.z;
    qDebug() << "ECEF position (km):" << pos_ecef.x << pos_ecef.y << pos_ecef.z;
    qDebug() << "ECI velocity (km/s):" << vel.x << vel.y << vel.z;
    qDebug() << "ECEF velocity (km/s):" << vel_ecef.x << vel_ecef.y << vel_ecef.z;

    return OrbitalState{pos_ecef, vel_ecef};

    return OrbitalState{pos, vel};
}
