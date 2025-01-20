// tle_parser.cpp
#include "tle_parser.h"
#include <QRegularExpression>

std::optional<TLE> TLEParser::parseFromTxt(const QString& line1, const QString& line2) {
    if (line1.length() != 69 || line2.length() != 69 || !checksum(line1) || !checksum(line2)) {
        return std::nullopt;
    }

    TLE tle;

    // First line parsing
    tle.satelliteNumber = line1.mid(2, 5).toInt();
    tle.classification = line1.at(7).toLatin1();
    tle.internationalDesignator = line1.mid(9, 8).trimmed().toStdString();
    tle.epochYear = line1.mid(18, 2).toInt();
    tle.epochDay = parseDecimal(line1.mid(20, 12));

    // Fix First Derivative parsing
    tle.firstDerivativeMeanMotion = parseDecimal(line1.mid(33, 10));

    // Fix Second Derivative parsing
    QString secondDerivStr = line1.mid(44, 6).trimmed();
    int secondDerivExp = line1.mid(50, 2).toInt();
    tle.secondDerivativeMeanMotion = parseDecimal(secondDerivStr) * std::pow(10.0, secondDerivExp);

    // Fix B* drag term parsing
    QString bstarStr = line1.mid(53, 6).trimmed();
    int bstarExp = line1.mid(59, 2).toInt();
    // Учитываем, что значение должно быть представлено как 0.81581 × 10^-3
    tle.bstar = (parseDecimal("0." + bstarStr) * std::pow(10.0, bstarExp));

    tle.ephemerisType = line1.at(62).digitValue();
    tle.elementNumber = line1.mid(64, 4).toInt();

    // Second line parsing
    tle.inclination = parseDecimal(line2.mid(8, 8).trimmed());
    tle.rightAscension = parseDecimal(line2.mid(17, 8).trimmed());
    tle.eccentricity = parseDecimal("0." + line2.mid(26, 7).trimmed());
    tle.argumentPerigee = parseDecimal(line2.mid(34, 8).trimmed());
    tle.meanAnomaly = parseDecimal(line2.mid(43, 8).trimmed());
    tle.meanMotion = parseDecimal(line2.mid(52, 11).trimmed());
    tle.revolutionNumber = line2.mid(63, 5).toInt();

    // Initialize computed values
    const double TWO_PI = 2.0 * M_PI;
    const double MINUTES_PER_DAY = 1440.0;

    // Convert mean motion to radians per minute
    tle.no = tle.meanMotion * TWO_PI / MINUTES_PER_DAY;

    // Calculate semi-major axis
    const double XKE = 0.743669161E-1;
    tle.a = std::pow(XKE / tle.no, 2.0/3.0);

    // Calculate apogee and perigee
    tle.alta = tle.a * (1.0 + tle.eccentricity) - 1.0;
    tle.altp = tle.a * (1.0 - tle.eccentricity) - 1.0;

    // Initialize jo (mean motion in radians per minute)
    tle.jo = tle.no;

    return tle;
}

std::optional<TLE> TLEParser::parseFromJson(const QJsonObject& json) {
    if (!json.contains("tle_line1") || !json.contains("tle_line2")) {
        return std::nullopt;
    }

    return parseFromTxt(json["tle_line1"].toString(), json["tle_line2"].toString());
}

std::optional<TLE> TLEParser::parseFromXml(const QDomElement& element) {
    QString line1 = element.firstChildElement("line1").text();
    QString line2 = element.firstChildElement("line2").text();

    if (line1.isEmpty() || line2.isEmpty()) {
        return std::nullopt;
    }

    return parseFromTxt(line1, line2);
}

bool TLEParser::checksum(const QString& line) {
    int sum = 0;
    for (int i = 0; i < 68; i++) {
        QChar c = line[i];
        if (c.isDigit()) {
            sum += c.digitValue();
        } else if (c == '-') {
            sum += 1;
        }
    }
    return (sum % 10) == (line[68].digitValue());
}

double TLEParser::parseDecimal(const QString& str) {
    QString trimmed = str.trimmed();
    if (trimmed.isEmpty()) {
        return 0.0;
    }

    // Заменяем несколько пробелов одним
    trimmed = trimmed.simplified();

    // Обработка знака, если он отделен пробелом
    if (trimmed.startsWith(" -") || trimmed.startsWith(" +")) {
        trimmed = trimmed.mid(1);
    }

    bool ok;
    double value = trimmed.toDouble(&ok);

    if (!ok) {
        return 0.0;
    }

    return value;
}
