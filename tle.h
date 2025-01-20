// tle.h
#pragma once

#include <string>
#include <array>
#include <chrono>
#include <optional>
#include <QDateTime>

struct TLE {
    // First line
    int satelliteNumber;
    char classification;
    std::string internationalDesignator;
    int epochYear;
    double epochDay;
    double firstDerivativeMeanMotion;
    double secondDerivativeMeanMotion;
    double bstar;
    int ephemerisType;
    int elementNumber;

    // Second line
    double inclination;        // degrees
    double rightAscension;    // degrees
    double eccentricity;
    double argumentPerigee;   // degrees
    double meanAnomaly;       // degrees
    double meanMotion;        // revolutions per day
    int revolutionNumber;

    // Computed values
    double no;                // rad/min
    double a;                 // earth radii
    double alta;             // earth radii
    double altp;             // earth radii
    double jo;               // rad/min
};

struct Vector3 {
    double x;
    double y;
    double z;

    Vector3(double x_ = 0.0, double y_ = 0.0, double z_ = 0.0)
        : x(x_), y(y_), z(z_) {}
};

struct OrbitalState {
    Vector3 position;     // km
    Vector3 velocity;     // km/s
};
