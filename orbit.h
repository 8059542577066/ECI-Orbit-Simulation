#ifndef ORBIT_H
#define ORBIT_H

#ifdef __WIN32__
#ifdef BUILD_LIB
#define LIB_CLASS __declspec(dllexport)
#else
#define LIB_CLASS __declspec(dllimport)
#endif
#else
#define LIB_CLASS
#endif

#include <vector>
#include <set>


#ifndef NOT_ANALYTICAL
#define NOT_ANALYTICAL 1
#endif

#ifndef NOT_ENOUGH_TICKS
#define NOT_ENOUGH_TICKS 2
#endif

#ifndef FRAGS_TOO_LITTLE
#define FRAGS_TOO_LITTLE 3
#endif


struct LIB_CLASS Vector3D
{
    double x, y, z;

    Vector3D();
    Vector3D(double, double, double);

    struct OrderByX
    {
        bool operator()(
            std::vector<Vector3D>::const_iterator,
            std::vector<Vector3D>::const_iterator);
    };
    struct OrderByY
    {
        bool operator()(
            std::vector<Vector3D>::const_iterator,
            std::vector<Vector3D>::const_iterator);
    };
    struct OrderByZ
    {
        bool operator()(
            std::vector<Vector3D>::const_iterator,
            std::vector<Vector3D>::const_iterator);
    };

    Vector3D operator-() const;

    Vector3D operator+(const Vector3D &) const;
    Vector3D operator-(const Vector3D &) const;
    Vector3D operator*(double) const;
    Vector3D operator/(double) const;

    Vector3D &operator+=(const Vector3D &);
    Vector3D &operator-=(const Vector3D &);
    Vector3D &operator*=(double);
    Vector3D &operator/=(double);

    double dot(const Vector3D &) const;
    double size() const;
    double distance(const Vector3D &,
                    const Vector3D &) const;

    Vector3D cross(const Vector3D &) const;
    Vector3D project(const Vector3D &) const;

    Vector3D rotateX(double) const;
    Vector3D rotateY(double) const;
    Vector3D rotateZ(double) const;
};


class LIB_CLASS Orbit
{
    double stdGrvPrm,
        semMajAxis, eccen, period,
        incl, longAscNode, argPeri;
    Vector3D spcAngMnt;

    std::vector<Vector3D> points;
    std::set<std::vector<Vector3D>::const_iterator,
             Vector3D::OrderByX>
        pointsByX;
    std::set<std::vector<Vector3D>::const_iterator,
             Vector3D::OrderByY>
        pointsByY;
    std::set<std::vector<Vector3D>::const_iterator,
             Vector3D::OrderByZ>
        pointsByZ;

    std::vector<Vector3D>::const_iterator
    findClosestPoint(const Vector3D &, double);

    double distance(const Vector3D &, double);
    double meanAnomRight(double) const;
    double eccAnom(double) const;

    Vector3D accel(const Vector3D &) const;

public:
    Orbit(double,
          Vector3D, Vector3D,
          double, bool);
    Orbit(double,
          const Vector3D &, const Vector3D &,
          std::size_t);

    double getMeanDist(const Orbit &);
    double getFinalDist(const Orbit &);

    Vector3D getPosVector(double) const;
    Vector3D getVelVector(double) const;

    void saveElements(const char *) const;
    void savePoints(const char *) const;
};


#endif