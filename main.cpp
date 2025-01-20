// main.cpp
#include <QApplication>
#include <QMainWindow>
#include <QTimer>
#include <QVBoxLayout>
#include <QLabel>
#include <QDateTime>
#include <QGroupBox>
#include <cmath>
#include "sgp4.h"
#include "tle_parser.h"

class SatelliteTracker : public QMainWindow {
    Q_OBJECT

public:
    SatelliteTracker(const TLE& tle, QWidget* parent = nullptr)
        : QMainWindow(parent), propagator_(tle) {
        setupUi();

        // Создаем таймер для обновления каждую секунду
        timer_ = new QTimer(this);
        connect(timer_, &QTimer::timeout, this, &SatelliteTracker::updatePosition);
        timer_->start(1000); // Обновление каждую секунду

        updatePosition(); // Первоначальное обновление
    }

private:
    void setupUi() {
        QWidget* centralWidget = new QWidget(this);
        setCentralWidget(centralWidget);

        QVBoxLayout* mainLayout = new QVBoxLayout(centralWidget);

        // Группа для времени
        QGroupBox* timeGroup = new QGroupBox("Time", centralWidget);
        QVBoxLayout* timeLayout = new QVBoxLayout(timeGroup);
        utcTimeLabel_ = new QLabel(timeGroup);
        timeLayout->addWidget(utcTimeLabel_);
        mainLayout->addWidget(timeGroup);

        // Группа для координат в ECEF
        QGroupBox* ecefGroup = new QGroupBox("ECEF Coordinates", centralWidget);
        QVBoxLayout* ecefLayout = new QVBoxLayout(ecefGroup);
        posXLabel_ = new QLabel(ecefGroup);
        posYLabel_ = new QLabel(ecefGroup);
        posZLabel_ = new QLabel(ecefGroup);
        ecefLayout->addWidget(posXLabel_);
        ecefLayout->addWidget(posYLabel_);
        ecefLayout->addWidget(posZLabel_);
        mainLayout->addWidget(ecefGroup);

        // Группа для географических координат
        QGroupBox* geoGroup = new QGroupBox("Geographic Coordinates", centralWidget);
        QVBoxLayout* geoLayout = new QVBoxLayout(geoGroup);
        latitudeLabel_ = new QLabel(geoGroup);
        longitudeLabel_ = new QLabel(geoGroup);
        altitudeLabel_ = new QLabel(geoGroup);
        geoLayout->addWidget(latitudeLabel_);
        geoLayout->addWidget(longitudeLabel_);
        geoLayout->addWidget(altitudeLabel_);
        mainLayout->addWidget(geoGroup);

        // Группа для скорости
        QGroupBox* velocityGroup = new QGroupBox("Velocity", centralWidget);
        QVBoxLayout* velocityLayout = new QVBoxLayout(velocityGroup);
        velocityLabel_ = new QLabel(velocityGroup);
        velocityLayout->addWidget(velocityLabel_);
        mainLayout->addWidget(velocityGroup);

        setWindowTitle("Satellite Tracker");
        resize(400, 600);
    }

    // Преобразование ECEF координат в географические (широта, долгота, высота)
    struct GeographicCoords {
        double latitude;  // градусы
        double longitude; // градусы
        double altitude;  // км
    };

    GeographicCoords ecefToGeographic(const Vector3& ecef) {
        qDebug() << "\n=== Geographic Coordinates Calculation Debug ===";
        qDebug() << "Input ECEF coordinates (km):";
        qDebug() << "X:" << ecef.x;
        qDebug() << "Y:" << ecef.y;
        qDebug() << "Z:" << ecef.z;

        // WGS84 ellipsoid parameters
        const double a = 6378.137; // semi-major axis (km)
        const double f = 1.0/298.257223563; // flattening
        const double b = a * (1.0 - f); // semi-minor axis (km)
        const double e = sqrt(1.0 - (b*b)/(a*a)); // eccentricity

        qDebug() << "\n--- WGS84 Parameters ---";
        qDebug() << "Semi-major axis (km):" << a;
        qDebug() << "Flattening:" << f;
        qDebug() << "Semi-minor axis (km):" << b;
        qDebug() << "Eccentricity:" << e;

        // Calculate intermediate values
        const double p = sqrt(ecef.x*ecef.x + ecef.y*ecef.y);
        const double theta = atan2(a*ecef.z, b*p);

        qDebug() << "\n--- Intermediate Values ---";
        qDebug() << "p (distance from Z axis) (km):" << p;
        qDebug() << "theta (parametric latitude) (rad):" << theta;

        // Calculate longitude
        const double longitude = atan2(ecef.y, ecef.x);

        // Calculate latitude using improved algorithm
        const double latitude = atan2(
            ecef.z + (e*e*(a/b))*a*pow(sin(theta), 3),
            p - (e*e)*a*pow(cos(theta), 3)
            );

        // Calculate altitude
        const double N = a / sqrt(1.0 - e*e*sin(latitude)*sin(latitude));
        const double altitude = p/cos(latitude) - N;

        qDebug() << "\n--- Final Geographic Coordinates ---";
        qDebug() << "Latitude (rad):" << latitude;
        qDebug() << "Longitude (rad):" << longitude;
        qDebug() << "Altitude (km):" << altitude;

        // Convert to degrees
        GeographicCoords geo;
        geo.latitude = latitude * 180.0 / M_PI;
        geo.longitude = longitude * 180.0 / M_PI;
        geo.altitude = altitude;

        qDebug() << "\n--- Final Coordinates (Degrees) ---";
        qDebug() << "Latitude (deg):" << geo.latitude;
        qDebug() << "Longitude (deg):" << geo.longitude;
        qDebug() << "Altitude (km):" << geo.altitude;

        // Validation checks
        qDebug() << "\n--- Validation ---";
        qDebug() << "Latitude in range [-90,90]:" << (geo.latitude >= -90 && geo.latitude <= 90);
        qDebug() << "Longitude in range [-180,180]:" << (geo.longitude >= -180 && geo.longitude <= 180);
        qDebug() << "Reasonable altitude (100-2000 km):" << (geo.altitude >= 100 && geo.altitude <= 2000);

        // Calculate ground track speed
        double r = sqrt(ecef.x*ecef.x + ecef.y*ecef.y + ecef.z*ecef.z);
        qDebug() << "Orbital radius (km):" << r;
        qDebug() << "Height above Earth's surface (km):" << (r - a);

        qDebug() << "=== End Geographic Coordinates Calculation ===\n";

        return geo;
    }

private slots:
    void updatePosition() {
        QDateTime currentTime = QDateTime::currentDateTimeUtc();
        utcTimeLabel_->setText("UTC: " + currentTime.toString("yyyy-MM-dd HH:mm:ss"));

        // Получаем текущее положение спутника
        OrbitalState state = propagator_.getPosition(currentTime);

        // Обновляем ECEF координаты
        posXLabel_->setText(QString("X: %1 km").arg(state.position.x, 0, 'f', 3));
        posYLabel_->setText(QString("Y: %1 km").arg(state.position.y, 0, 'f', 3));
        posZLabel_->setText(QString("Z: %1 km").arg(state.position.z, 0, 'f', 3));

        // Преобразуем в географические координаты
        GeographicCoords geo = ecefToGeographic(state.position);

        // Обновляем метки с географическими координатами
        latitudeLabel_->setText(QString("Latitude: %1°").arg(geo.latitude, 0, 'f', 6));
        longitudeLabel_->setText(QString("Longitude: %1°").arg(geo.longitude, 0, 'f', 6));
        altitudeLabel_->setText(QString("Altitude: %1 km").arg(geo.altitude, 0, 'f', 3));

        // Вычисляем скорость
        double velocity = sqrt(
            state.velocity.x*state.velocity.x +
            state.velocity.y*state.velocity.y +
            state.velocity.z*state.velocity.z
            );
        velocityLabel_->setText(QString("Velocity: %1 km/s").arg(velocity, 0, 'f', 3));
    }

private:
    SGP4 propagator_;
    QTimer* timer_;
    QLabel* utcTimeLabel_;
    QLabel* posXLabel_;
    QLabel* posYLabel_;
    QLabel* posZLabel_;
    QLabel* latitudeLabel_;
    QLabel* longitudeLabel_;
    QLabel* altitudeLabel_;
    QLabel* velocityLabel_;
};

void debugPrintTLE(const TLE& tle) {
    qDebug() << "\n=== TLE Values After Parsing ===";

    qDebug() << "\n--- First Line Values ---";
    qDebug() << "Satellite Number:" << tle.satelliteNumber;
    qDebug() << "Classification:" << tle.classification;
    qDebug() << "International Designator:" << QString::fromStdString(tle.internationalDesignator);
    qDebug() << "Epoch Year:" << tle.epochYear;
    qDebug() << "Epoch Day:" << tle.epochDay;
    qDebug() << "First Derivative Mean Motion:" << tle.firstDerivativeMeanMotion;
    qDebug() << "Second Derivative Mean Motion:" << tle.secondDerivativeMeanMotion;
    qDebug() << "B* drag term:" << tle.bstar;
    qDebug() << "Ephemeris Type:" << tle.ephemerisType;
    qDebug() << "Element Number:" << tle.elementNumber;

    qDebug() << "\n--- Second Line Values ---";
    qDebug() << "Inclination (degrees):" << tle.inclination;
    qDebug() << "Right Ascension (degrees):" << tle.rightAscension;
    qDebug() << "Eccentricity:" << tle.eccentricity;
    qDebug() << "Argument of Perigee (degrees):" << tle.argumentPerigee;
    qDebug() << "Mean Anomaly (degrees):" << tle.meanAnomaly;
    qDebug() << "Mean Motion (revs per day):" << tle.meanMotion;
    qDebug() << "Revolution Number:" << tle.revolutionNumber;

    qDebug() << "\n--- Computed Values ---";
    qDebug() << "no (rad/min):" << tle.no;
    qDebug() << "a (earth radii):" << tle.a;
    qDebug() << "alta (earth radii):" << tle.alta;
    qDebug() << "altp (earth radii):" << tle.altp;
    qDebug() << "jo (rad/min):" << tle.jo;

    qDebug() << "\n--- Original TLE Lines ---";
    qDebug() << "Line 1:" << "1 57890U 23145E   25020.52408874  .00019534  00000-0  81581-3 0  9999";
    qDebug() << "Line 2:" << "2 57890  34.9931 184.3980 0003522 251.7861 108.2470 15.22764851 74908";

    qDebug() << "=== End TLE Values ===\n";
}

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    // Пример TLE данных (МКС)
    QString line1 = "1 57890U 23145E   25020.52408874  .00019534  00000-0  81581-3 0  9999";
    QString line2 = "2 57890  34.9931 184.3980 0003522 251.7861 108.2470 15.22764851 74908";

    auto tle_opt = TLEParser::parseFromTxt(line1, line2);
    if (!tle_opt) {
        qDebug() << "Failed to parse TLE";
        return 1;
    }

    // Добавьте эту строку
    debugPrintTLE(*tle_opt);

    // Создаем и показываем окно трекера
    SatelliteTracker tracker(*tle_opt);
    tracker.show();

    return app.exec();
}

// Добавьте в конец файла для поддержки Q_OBJECT
#include "main.moc"
