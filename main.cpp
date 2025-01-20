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
        const double a = 6378.137; // большая полуось эллипсоида WGS84 (км)
        const double f = 1.0/298.257223563; // сжатие WGS84
        const double b = a * (1.0 - f); // малая полуось
        const double e2 = 2*f - f*f; // квадрат эксцентриситета

        const double p = sqrt(ecef.x*ecef.x + ecef.y*ecef.y);
        const double theta = atan2(ecef.z*a, p*b);

        // Вычисление широты
        const double latitude = atan2(
            ecef.z + e2 * (a/b) * a * pow(sin(theta), 3),
            p - e2 * a * pow(cos(theta), 3)
            );

        // Вычисление долготы
        double longitude = atan2(ecef.y, ecef.x);

        // Преобразование долготы в диапазон [-180, 180]
        longitude = longitude * 180.0 / M_PI;
        while (longitude > 180.0) longitude -= 360.0;
        while (longitude < -180.0) longitude += 360.0;

        // Вычисление высоты
        const double sin_lat = sin(latitude);
        const double N = a / sqrt(1.0 - e2 * sin_lat * sin_lat);
        const double altitude = p/cos(latitude) - N;

        GeographicCoords geo;
        geo.latitude = latitude * 180.0 / M_PI;  // широта в градусах
        geo.longitude = longitude;                // долгота уже в градусах
        geo.altitude = altitude;

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
        QString latStr = QString("Latitude: %1° %2").arg(fabs(geo.latitude), 0, 'f', 6)
                             .arg(geo.latitude < 0 ? "S" : "N");
        QString lonStr = QString("Longitude: %1° %2").arg(fabs(geo.longitude), 0, 'f', 6)
                             .arg(geo.longitude < 0 ? "W" : "E");

        latitudeLabel_->setText(latStr);
        longitudeLabel_->setText(lonStr);
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
    qDebug() << "Line 1:" << "1 57890U 23145E   25020.58965078  .00018142  00000-0  75795-3 0  9993";
    qDebug() << "Line 2:" << "2 57890  34.9931 183.9858 0003532 252.3040 107.7289 15.22765489 74915";

    qDebug() << "=== End TLE Values ===\n";
}

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    // Пример TLE данных (МКС)
    QString line1 = "1 57890U 23145E   25020.58965078  .00018142  00000-0  75795-3 0  9993";
    QString line2 = "2 57890  34.9931 183.9858 0003532 252.3040 107.7289 15.22765489 74915";

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
