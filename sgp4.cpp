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

    // Расчет времени с эпохи
    int year = tle_.epochYear;
    if (year < 57) {
        year += 2000;
    } else {
        year += 1900;
    }

    QDateTime epoch = QDateTime(QDate(year, 1, 1), QTime(0, 0), Qt::UTC);
    epoch = epoch.addSecs(static_cast<qint64>((tle_.epochDay - 1.0) * 86400.0));
    double tsince = epoch.secsTo(time) / 60.0;  // время в минутах

    qDebug() << "Epoch time:" << epoch.toString("yyyy-MM-dd HH:mm:ss");
    qDebug() << "Time since epoch (minutes):" << tsince;

    // Вычисление средней аномалии, аргумента перигея и прямого восхождения
    double xmo = fmod(tle_.meanAnomaly * M_PI / 180.0 + xmdot_ * tsince, 2.0 * M_PI);
    double omegao = tle_.argumentPerigee * M_PI / 180.0 + omgdot_ * tsince;
    double xno = tle_.rightAscension * M_PI / 180.0 + xnodot_ * tsince;

    qDebug() << "\n--- Orbital Elements at Time t ---";
    qDebug() << "Mean Anomaly (rad):" << xmo;
    qDebug() << "Argument of Perigee (rad):" << omegao;
    qDebug() << "Right Ascension (rad):" << xno;

    // Улучшенное решение уравнения Кеплера
    double e = tle_.eccentricity;
    double E = xmo;  // начальное приближение
    double delta;
    int iter = 0;
    const int maxIter = 30;
    const double tolerance = 1e-12;

    qDebug() << "\n--- Kepler's Equation Solution ---";
    do {
        double E_old = E;
        double sinE = sin(E);
        double cosE = cos(E);

        // Метод Ньютона-Рафсона с регуляризацией
        delta = (E - e * sinE - xmo) / (1.0 - e * cosE);
        if (fabs(delta) > 0.5) {
            delta = copysign(0.5, delta);
        }
        E = E_old - delta;

        iter++;
        qDebug() << "Iteration" << iter << "- E:" << E << "Delta:" << fabs(delta);
    } while (fabs(delta) > tolerance && iter < maxIter);

    // Вычисление истинной аномалии
    double sinE = sin(E);
    double cosE = cos(E);
    double sinv = sqrt(1.0 - e * e) * sinE / (1.0 - e * cosE);
    double cosv = (cosE - e) / (1.0 - e * cosE);
    double nu = atan2(sinv, cosv);

    // Вычисление позиции в орбитальной плоскости
    double r = tle_.a * XKMPER * (1.0 - e * cosE);  // радиус в км
    double u = omegao + nu;  // аргумент широты

    qDebug() << "\n--- Position Calculation ---";
    qDebug() << "True Anomaly (rad):" << nu;
    qDebug() << "Radius vector (km):" << r;
    qDebug() << "Argument of latitude (rad):" << u;

    // Тригонометрические функции для преобразования координат
    double sin_u = sin(u);
    double cos_u = cos(u);
    double sin_i = sinio_;
    double cos_i = cosio_;
    double sin_node = sin(xno);
    double cos_node = cos(xno);

    // Позиция в ECI координатах
    Vector3 pos_eci;
    pos_eci.x = r * (cos_node * cos_u - sin_node * sin_u * cos_i);
    pos_eci.y = r * (sin_node * cos_u + cos_node * sin_u * cos_i);
    pos_eci.z = r * sin_u * sin_i;

    // Вычисление скорости в ECI
    double xn = xnodp_;
    double rdot = xn * tle_.a * XKMPER * e * sinE / sqrt(1.0 - e * e);
    double rfdot = xn * tle_.a * XKMPER * sqrt(1.0 - e * e) / (1.0 - e * cosE);

    Vector3 vel_eci;
    vel_eci.x = rdot * (cos_node * cos_u - sin_node * sin_u * cos_i) -
                rfdot * (cos_node * sin_u + sin_node * cos_u * cos_i);
    vel_eci.y = rdot * (sin_node * cos_u + cos_node * sin_u * cos_i) -
                rfdot * (sin_node * sin_u - cos_node * cos_u * cos_i);
    vel_eci.z = rdot * sin_u * sin_i + rfdot * cos_u * sin_i;

    // Учет возмущений
    calculatePerturbations(tsince, pos_eci, vel_eci);

    // Вычисление GMST
    double gmst = calculateGMST(time);

    // Преобразование в ECEF координаты
    Vector3 pos_ecef;
    pos_ecef.x = pos_eci.x * cos(gmst) + pos_eci.y * sin(gmst);
    pos_ecef.y = -pos_eci.x * sin(gmst) + pos_eci.y * cos(gmst);
    pos_ecef.z = pos_eci.z;

    // Преобразование скорости в ECEF
    const double earth_rotation_rate = 7.2921150e-5;  // рад/с
    Vector3 vel_ecef;
    vel_ecef.x = vel_eci.x * cos(gmst) + vel_eci.y * sin(gmst) -
                 earth_rotation_rate * pos_ecef.y;
    vel_ecef.y = -vel_eci.x * sin(gmst) + vel_eci.y * cos(gmst) +
                 earth_rotation_rate * pos_ecef.x;
    vel_ecef.z = vel_eci.z;

    // Учет полярного движения
    PolarMotion pm = {0.1, 0.2};  // Примерные значения, нужно получать актуальные
    pos_ecef = applyPolarMotion(pos_ecef, pm);

    qDebug() << "\n--- Final Position and Velocity ---";
    qDebug() << "ECEF Position (km):" << pos_ecef.x << pos_ecef.y << pos_ecef.z;
    qDebug() << "ECEF Velocity (km/s):" << vel_ecef.x << vel_ecef.y << vel_ecef.z;

    return OrbitalState{pos_ecef, vel_ecef};
}

// Добавим новую функцию в класс SGP4
double SGP4::calculateGMST(const QDateTime& time) const {
    QDateTime j2000(QDate(2000, 1, 1), QTime(12, 0), Qt::UTC);
    double days_since_j2000 = j2000.secsTo(time) / 86400.0;
    double T = days_since_j2000 / 36525.0;  // юлианские столетия

    // GMST в градусах, учитывая высокоточные члены
    double gmst = 280.46061837 +
                  360.98564736629 * days_since_j2000 +
                  0.000387933 * T * T -
                  T * T * T / 38710000.0;

    // Учитываем прецессию
    double zeta = (2306.2181 * T + 0.30188 * T * T + 0.017998 * T * T * T) / 3600.0;
    double z = (2306.2181 * T + 1.09468 * T * T + 0.018203 * T * T * T) / 3600.0;
    double theta = (2004.3109 * T - 0.42665 * T * T - 0.041833 * T * T * T) / 3600.0;

    gmst += zeta + z;

    // Добавляем вклад от нутации
    double omega = 125.04452222 - 1934.136261111 * T;  // долгота восходящего узла Луны
    double L = 280.4665 + 36000.7698 * T;             // средняя долгота Солнца
    double Lp = 218.3165 + 481267.8813 * T;           // средняя долгота Луны

    double dpsi = -17.20 * sin(omega * M_PI/180.0) - 1.32 * sin(2*L * M_PI/180.0)
                  - 0.23 * sin(2*Lp * M_PI/180.0) + 0.21 * sin(2*omega * M_PI/180.0);
    dpsi /= 3600.0;  // переводим из угловых секунд в градусы

    double eps = 23.43929111 - 0.013004167 * T;  // наклон эклиптики
    gmst += dpsi * cos(eps * M_PI/180.0);

    // Нормализуем в диапазон [0, 360]
    gmst = fmod(gmst, 360.0);
    if (gmst < 0) gmst += 360.0;

    qDebug() << "\n=== GMST Calculation Debug ===";
    qDebug() << "UTC time:" << time.toString("yyyy-MM-dd HH:mm:ss");
    qDebug() << "Julian centuries T:" << T;
    qDebug() << "Base GMST (degrees):" << gmst - zeta - z - dpsi * cos(eps * M_PI/180.0);
    qDebug() << "Precession (degrees):" << zeta + z;
    qDebug() << "Nutation (arcsec):" << dpsi;
    qDebug() << "Final GMST (degrees):" << gmst;

    return gmst * M_PI / 180.0;  // возвращаем в радианах
}

Vector3 SGP4::applyPolarMotion(const Vector3& ecef, const PolarMotion& pm) const {
    // Преобразование угловых секунд в радианы
    double xp = pm.xp * M_PI / (180.0 * 3600.0);
    double yp = pm.yp * M_PI / (180.0 * 3600.0);

    // Матрица поворота для учета полярного движения
    Vector3 result;
    result.x = ecef.x - ecef.y*xp + ecef.z*yp;
    result.y = ecef.x*xp + ecef.y - ecef.z*xp*yp;
    result.z = -ecef.x*yp + ecef.y*xp*yp + ecef.z;

    return result;
}
void SGP4::calculatePerturbations(double tsince, Vector3& pos, Vector3& vel) const {
    // Возмущения от J2
    double r = sqrt(pos.x*pos.x + pos.y*pos.y + pos.z*pos.z);
    double r2 = r * r;
    double r3 = r2 * r;
    double r4 = r3 * r;

    double z2 = pos.z * pos.z / r2;
    double j2Effect = 1.5 * CK2 / r2;

    // Возмущения положения
    pos.x *= (1.0 - j2Effect * (1.0 - 5.0 * z2));
    pos.y *= (1.0 - j2Effect * (1.0 - 5.0 * z2));
    pos.z *= (1.0 - j2Effect * (3.0 - 5.0 * z2));

    // Возмущения скорости
    double vDotR = (vel.x*pos.x + vel.y*pos.y + vel.z*pos.z) / r;
    double j2Vel = 1.5 * CK2 / r3;

    vel.x -= j2Vel * (pos.x/r * (1.0 - 5.0*z2) + 2.0*pos.z*pos.z/r2);
    vel.y -= j2Vel * (pos.y/r * (1.0 - 5.0*z2) + 2.0*pos.z*pos.z/r2);
    vel.z -= j2Vel * (pos.z/r * (3.0 - 5.0*z2) + 2.0*pos.z);
}

double SGP4::solveKepler(double M, double e, double tolerance) const {
    double E = M;
    double delta;
    int iter = 0;
    const int maxIter = 30;

    do {
        double E_old = E;
        double sinE = sin(E);
        double cosE = cos(E);

        // Метод Ньютона-Рафсона с регуляризацией
        delta = (E - e * sinE - M) / (1.0 - e * cosE);
        E = E_old - delta;

        // Предотвращение расходимости
        if (fabs(delta) > 0.5) {
            delta = copysign(0.5, delta);
            E = E_old - delta;
        }

        iter++;
    } while (fabs(delta) > tolerance && iter < maxIter);

    return E;
}
