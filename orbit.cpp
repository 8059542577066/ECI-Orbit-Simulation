#include "orbit.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>


Vector3D::Vector3D() {}

Vector3D::Vector3D(double x, double y, double z)
    : x(x), y(y), z(z) {}

bool Vector3D::OrderByX::operator()(
    std::vector<Vector3D>::const_iterator ptr1,
    std::vector<Vector3D>::const_iterator ptr2)
{
    return ptr1->x < ptr2->x;
}

bool Vector3D::OrderByY::operator()(
    std::vector<Vector3D>::const_iterator ptr1,
    std::vector<Vector3D>::const_iterator ptr2)
{
    return ptr1->y < ptr2->y;
}

bool Vector3D::OrderByZ::operator()(
    std::vector<Vector3D>::const_iterator ptr1,
    std::vector<Vector3D>::const_iterator ptr2)
{
    return ptr1->z < ptr2->z;
}

Vector3D Vector3D::operator-() const
{
    Vector3D result;
    result.x = -this->x;
    result.y = -this->y;
    result.z = -this->z;

    return result;
}

Vector3D Vector3D::operator+(const Vector3D &other) const
{
    Vector3D result;
    result.x = this->x + other.x;
    result.y = this->y + other.y;
    result.z = this->z + other.z;

    return result;
}

Vector3D Vector3D::operator-(const Vector3D &other) const
{
    Vector3D result;
    result.x = this->x - other.x;
    result.y = this->y - other.y;
    result.z = this->z - other.z;

    return result;
}

Vector3D Vector3D::operator*(double value) const
{
    Vector3D result;
    result.x = this->x * value;
    result.y = this->y * value;
    result.z = this->z * value;

    return result;
}

Vector3D Vector3D::operator/(double value) const
{
    Vector3D result;
    result.x = this->x / value;
    result.y = this->y / value;
    result.z = this->z / value;

    return result;
}

Vector3D &Vector3D::operator+=(const Vector3D &other)
{
    this->x += other.x;
    this->y += other.y;
    this->z += other.z;

    return *this;
}

Vector3D &Vector3D::operator-=(const Vector3D &other)
{
    this->x -= other.x;
    this->y -= other.y;
    this->z -= other.z;

    return *this;
}

Vector3D &Vector3D::operator*=(double value)
{
    this->x *= value;
    this->y *= value;
    this->z *= value;

    return *this;
}

Vector3D &Vector3D::operator/=(double value)
{
    this->x /= value;
    this->y /= value;
    this->z /= value;

    return *this;
}

double Vector3D::dot(const Vector3D &other) const
{
    return this->x * other.x + this->y * other.y + this->z * other.z;
}

double Vector3D::size() const
{
    return std::sqrt(this->dot(*this));
}

double Vector3D::distance(const Vector3D &point1,
                          const Vector3D &point2) const
{
    if ((point2 - point1).dot(*this - point1) < 0)
        return (*this - point1).size();
    else if ((point1 - point2).dot(*this - point2) < 0)
        return (*this - point2).size();
    else
        return (point1 + (*this - point1).project(point2 - point1) -
                *this)
            .size();
}

Vector3D Vector3D::cross(const Vector3D &other) const
{
    Vector3D result;
    result.x = this->y * other.z - this->z * other.y;
    result.y = this->z * other.x - this->x * other.z;
    result.z = this->x * other.y - this->y * other.x;

    return result;
}

Vector3D Vector3D::project(const Vector3D &other) const
{
    return other * (this->dot(other) / other.dot(other));
}

Vector3D Vector3D::rotateX(double angle) const
{
    Vector3D result;
    result.x = this->x;
    result.y = this->y * std::cos(angle) - this->z * std::sin(angle);
    result.z = this->z * std::cos(angle) + this->y * std::sin(angle);

    return result;
}

Vector3D Vector3D::rotateY(double angle) const
{
    Vector3D result;
    result.x = this->x * std::cos(angle) + this->z * std::sin(angle);
    result.y = this->y;
    result.z = this->z * std::cos(angle) - this->x * std::sin(angle);

    return result;
}

Vector3D Vector3D::rotateZ(double angle) const
{
    Vector3D result;
    result.x = this->x * std::cos(angle) - this->y * std::sin(angle);
    result.y = this->y * std::cos(angle) + this->x * std::sin(angle);
    result.z = this->z;

    return result;
}


std::vector<Vector3D>::const_iterator
Orbit::findClosestPoint(const Vector3D &point, double margin)
{
    std::vector<Vector3D> lower, upper;
    lower.emplace_back(point), upper.emplace_back(point);
    auto lowerZ = this->pointsByZ.cend(),
         upperZ = this->pointsByZ.cend();

    while (true)
    {
        this->pointsByY.clear(), this->pointsByZ.clear();

        lower[0].x -= margin, upper[0].x += margin;
        auto lowerX = this->pointsByX.lower_bound(lower.cbegin()),
             upperX = this->pointsByX.upper_bound(upper.cbegin());

        for (auto iterX = lowerX; iterX != upperX; ++iterX)
            this->pointsByY.emplace(*iterX);

        lower[0].y -= margin, upper[0].y += margin;
        auto lowerY = this->pointsByY.lower_bound(lower.cbegin()),
             upperY = this->pointsByY.upper_bound(upper.cbegin());

        for (auto iterY = lowerY; iterY != upperY; ++iterY)
            this->pointsByZ.emplace(*iterY);

        lower[0].z -= margin, upper[0].z += margin;
        lowerZ = this->pointsByZ.lower_bound(lower.cbegin());
        upperZ = this->pointsByZ.upper_bound(upper.cbegin());

        if (lowerZ != upperZ)
            break;
    }

    auto result = *lowerZ;
    double min = (point - *result).size();

    for (auto iterZ = lowerZ; iterZ != upperZ; ++iterZ)
    {
        double size = (point - **iterZ).size();

        if (size < min)
            result = *iterZ, min = size;
    }

    return result;
}

double Orbit::distance(const Vector3D &point, double margin)
{
    auto closest = this->findClosestPoint(point, margin);
    Vector3D prev, next;

    if (closest == this->points.cbegin())
        prev = *this->points.crbegin(), next = *(closest + 1);
    else if (closest + 1 == this->points.cend())
        prev = *(closest - 1), next = *this->points.cbegin();
    else
        prev = *(closest - 1), next = *(closest + 1);

    double distPrev = point.distance(*closest, prev),
           distNext = point.distance(*closest, next);

    return distPrev < distNext ? distPrev : distNext;
}

double Orbit::meanAnomRight(double eccAnom) const
{
    return eccAnom - this->eccen * std::sin(eccAnom);
}

double Orbit::eccAnom(double time) const
{
    double meanAnom = 2 * std::acos(-1) / this->period * time,
           x = 0, m = std::acos(-1);

    while (x + m != x)
    {
        if (this->meanAnomRight(x) < meanAnom)
            x += m;
        else
            x -= m;

        m /= 2;
    }

    return x;
}

Vector3D Orbit::accel(const Vector3D &pos) const
{
    return -pos * (this->stdGrvPrm / std::pow(pos.size(), 3));
}

Orbit::Orbit(double stdGrvPrm,
             Vector3D pos, Vector3D vel,
             double tps, bool lastOnly)
    : stdGrvPrm(stdGrvPrm)
{
    double posSize = pos.size(), velSize = vel.size();
    this->semMajAxis = stdGrvPrm * posSize;
    this->semMajAxis /= 2 * stdGrvPrm - posSize * velSize * velSize;

    this->period = 2 * std::acos(-1);
    this->period *= std::sqrt(std::pow(semMajAxis, 3) / stdGrvPrm);

    std::size_t ticks = this->period * 1.1 * tps;
    double deltaTime = 1.0 / tps;

    if (ticks == 0)
        throw NOT_ENOUGH_TICKS;

    for (std::size_t i = 0; i < ticks - 1; ++i)
    {
        if (i % 1000 == 0)
            std::cout << "Integrating Step " << i
                      << " / " << ticks << std::endl;

        if (!lastOnly)
            this->points.emplace_back(pos);

        Vector3D posRK4_1 = vel * deltaTime,
                 velRK4_1 = this->accel(pos) * deltaTime,
                 posRK4_2 = (vel + velRK4_1 / 2) * deltaTime,
                 velRK4_2 = this->accel(pos + posRK4_1 / 2) * deltaTime,
                 posRK4_3 = (vel + velRK4_2 / 2) * deltaTime,
                 velRK4_3 = this->accel(pos + posRK4_2 / 2) * deltaTime,
                 posRK4_4 = (vel + velRK4_3) * deltaTime,
                 velRK4_4 = this->accel(pos + posRK4_3) * deltaTime;
        posRK4_2 *= 2, velRK4_2 *= 2, posRK4_3 *= 2, velRK4_3 *= 2;

        Vector3D posRK4 = (posRK4_1 + posRK4_2 + posRK4_3 + posRK4_4),
                 velRK4 = (velRK4_1 + velRK4_2 + velRK4_3 + velRK4_4);

        pos += posRK4 / 6, vel += velRK4 / 6;
    }

    this->points.emplace_back(pos);
}

Orbit::Orbit(double stdGrvPrm,
             const Vector3D &pos, const Vector3D &vel,
             std::size_t fragments)
    : stdGrvPrm(stdGrvPrm)
{
    if (fragments < 16)
        throw FRAGS_TOO_LITTLE;

    double posSize = pos.size(), velSize = vel.size();
    this->semMajAxis = stdGrvPrm * posSize;
    this->semMajAxis /= 2 * stdGrvPrm - posSize * velSize * velSize;

    this->spcAngMnt = pos.cross(vel);
    double spcAngMntSize = this->spcAngMnt.size(),
           spcAngMntSizeSqr = spcAngMntSize * spcAngMntSize;
    this->eccen = std::sqrt(1 - spcAngMntSizeSqr /
                                    (stdGrvPrm * semMajAxis));
    this->period = 2 * std::acos(-1);
    this->period *= std::sqrt(std::pow(semMajAxis, 3) / stdGrvPrm);

    for (std::size_t i = 0; i < fragments; ++i)
    {
        double dist = spcAngMntSizeSqr,
               angle = 2 * std::acos(-1) * i / fragments;
        dist /= stdGrvPrm * (1 + this->eccen * std::cos(angle));

        this->points.emplace_back(dist * std::cos(angle),
                                  dist * std::sin(angle),
                                  0);
    }

    this->incl = std::acos(this->spcAngMnt.z / spcAngMntSize);
    this->longAscNode = std::atan2(this->spcAngMnt.y, this->spcAngMnt.x);
    Vector3D posPeriFoc = pos.rotateZ(-this->longAscNode)
                              .rotateY(-this->incl);
    double angle = std::acos(
        (spcAngMntSizeSqr / (stdGrvPrm * posSize) - 1) /
        this->eccen);
    this->argPeri = pos.dot(vel) >= 0
                        ? std::atan2(posPeriFoc.y, posPeriFoc.x) - angle
                        : std::atan2(posPeriFoc.y, posPeriFoc.x) + angle;

    for (std::size_t i = 0; i < this->points.size(); ++i)
    {
        this->points[i] = this->points[i]
                              .rotateZ(this->argPeri)
                              .rotateY(this->incl)
                              .rotateZ(this->longAscNode);

        this->pointsByX.emplace(&this->points[i]);
    }
}

double Orbit::getMeanDist(const Orbit &other)
{
    if (this->pointsByX.size() == 0)
        throw NOT_ANALYTICAL;

    double summed = 0, max = 0,
           margin = (this->points[0] - this->points[1]).size() * 4;

    std::ofstream fout;
    fout.open("errors.txt");
    fout << std::setprecision(14);
    std::cout << std::setprecision(14);

    for (std::size_t i = 0; i < other.points.size(); ++i)
    {
        double error = this->distance(other.points[i], margin);
        fout << "Step " << i << ": " << error << std::endl;
        std::cout << "Step " << i << ": " << error << std::endl;
        summed += error, max = error > max ? error : max;
    }

    fout << "\nAverage over " << other.points.size()
         << " steps: " << summed / other.points.size();
    fout << "\nMaximum Error: " << max;
    fout.close();

    return summed / other.points.size();
}

double Orbit::getFinalDist(const Orbit &other)
{
    if (this->pointsByX.size() == 0)
        throw NOT_ANALYTICAL;

    double margin = (this->points[0] - this->points[1]).size() * 4,
           error = this->distance(*other.points.crbegin(), margin);

    return error;
}

Vector3D Orbit::getPosVector(double time) const
{
    double eccAnom = this->eccAnom(std::fmod(time, this->period));

    Vector3D result;
    result.x = std::cos(eccAnom) - this->eccen;
    result.y = std::sqrt(1 - this->eccen * this->eccen) *
               std::sin(eccAnom);
    result.z = 0;
    result *= this->semMajAxis;

    return result.rotateZ(this->argPeri)
        .rotateY(this->incl)
        .rotateZ(this->longAscNode);
}

Vector3D Orbit::getVelVector(double time) const
{
    double eccAnom = this->eccAnom(std::fmod(time, this->period)),
           dist = this->semMajAxis *
                  (1 - this->eccen * std::cos(eccAnom));

    Vector3D result;
    result.x = -std::sin(eccAnom);
    result.y = std::sqrt(1 - this->eccen * this->eccen) *
               std::cos(eccAnom);
    result.z = 0;
    result *= std::sqrt(this->stdGrvPrm * this->semMajAxis) / dist;

    return result.rotateZ(this->argPeri)
        .rotateY(this->incl)
        .rotateZ(this->longAscNode);
}

void Orbit::saveElements(const char *name) const
{
    std::stringstream ss;
    ss << name << ".txt";
    std::ofstream fout;
    fout.open(ss.str());
    fout << std::setprecision(14);

    fout << "Semi-major Axis: " << this->semMajAxis;
    fout << "\nEccentricity: " << this->eccen;
    fout << "\nOrbital Period: " << this->period;
    fout << "\nInclination: " << this->incl;
    fout << "\nLongitude of Ascending Node: " << this->longAscNode;
    fout << "\nArgument of Periapsis: " << this->argPeri;

    fout.close();
}

void Orbit::savePoints(const char *name) const
{
    std::stringstream ss;
    ss << name << ".py";
    std::ofstream fout;
    fout.open(ss.str());
    fout << std::setprecision(14) << name << " = []\n";

    for (std::size_t i = 0; i < this->points.size(); ++i)
        fout << name << ".append(("
             << this->points[i].x << ", "
             << this->points[i].y << ", "
             << this->points[i].z << "))\n";

    fout.close();
}