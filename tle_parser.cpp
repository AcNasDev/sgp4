// tle_parser.cpp
#include "tle_parser.h"
#include <QRegularExpression>

std::optional<TLE> TLEParser::parseFromTxt(const QString& line1, const QString& line2) {
    if (line1.length() != 69 || line2.length() != 69 || !checksum(line1) || !checksum(line2)) {
        return std::nullopt;
    }

    TLE tle;

    // Parse first line
    tle.satelliteNumber = line1.mid(2, 5).toInt();
    tle.classification = line1.at(7).toLatin1();
    tle.internationalDesignator = line1.mid(9, 8).trimmed().toStdString();
    tle.epochYear = line1.mid(18, 2).toInt();
    tle.epochDay = parseDecimal(line1.mid(20, 12));
    tle.firstDerivativeMeanMotion = parseDecimal(line1.mid(33, 10));
    tle.secondDerivativeMeanMotion = parseDecimal("0." + line1.mid(44, 6));
    tle.bstar = parseDecimal("0." + line1.mid(53, 6)) * std::pow(10, line1.mid(59, 2).toInt());

    // Parse second line
    tle.inclination = parseDecimal(line2.mid(8, 8));
    tle.rightAscension = parseDecimal(line2.mid(17, 8));
    tle.eccentricity = parseDecimal("0." + line2.mid(26, 7));
    tle.argumentPerigee = parseDecimal(line2.mid(34, 8));
    tle.meanAnomaly = parseDecimal(line2.mid(43, 8));
    tle.meanMotion = parseDecimal(line2.mid(52, 11));
    tle.revolutionNumber = line2.mid(63, 5).toInt();

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
    static QRegularExpression re("^[\\-\\+]?\\d*\\.?\\d*$");
    if (!re.match(str).hasMatch()) {
        return 0.0;
    }
    return str.toDouble();
}
