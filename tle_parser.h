// tle_parser.h
#pragma once

#include "tle.h"
#include <QString>
#include <QJsonObject>
#include <QDomElement>
#include <optional>

class TLEParser {
public:
    static std::optional<TLE> parseFromTxt(const QString& line1, const QString& line2);
    static std::optional<TLE> parseFromJson(const QJsonObject& json);
    static std::optional<TLE> parseFromXml(const QDomElement& element);

private:
    static bool checksum(const QString& line);
    static double parseDecimal(const QString& str);
};
