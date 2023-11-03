#ifndef UTIL_MATH_H
#define UTIL_MATH_H

#define T_STEP 0.1
#define SIMILARITY_TOLERANCE 0.05;

#include <math.h>
#include <iostream>

struct point_t
{
public:
    point_t() { x = 0.0; y = 0.0; z = 0.0; }
    point_t(double i, double j) { x = i; y = j; z = 0; }
    point_t(double i, double j, double k) { x = i; y = j; z = k; }
    virtual ~point_t() {}

    bool is_zero() const { return ( x == 0 && y == 0 && z == 0); }

    double L2_norm3D() const { return sqrt(x * x + y * y + z * z);}
    double L2_norm2D() const { return sqrt(x * x + y * y); }
    void normalize() { double norm = sqrt(x * x + y * y + z * z); x /= norm; y /= norm; z /= norm; }
    bool similar_to(const point_t& P) const { return (*this - P).L2_norm3D() < SIMILARITY_TOLERANCE; }

    const point_t operator -(const point_t P) const { return point_t(x - P.x, y - P.y, z - P.z); }
    const point_t operator +(const point_t P) const { return point_t(x + P.x, y + P.y, z + P.z); }
    const point_t operator *(double a) const { return point_t(a * x, a * y, a * z); }
    const point_t operator /(double a) const { return a != 0 ? point_t(x/a, y/a, z/a) : point_t(x,y,z); }
    const point_t operator -() const { return point_t(-x, -y, -z); }
    void operator +=(const point_t P) { x += P.x; y += P.y; z += P.z; }
    void operator -=(const point_t P) { x -= P.x; y -= P.y; z -= P.z; }

    bool operator !=(point_t& P) { if (x != P.x || y != P.y || z != P.z) return true; return false; }
    bool operator ==(point_t& P) { if (x == P.x && y == P.y && z == P.z) return true; return false; }

    friend std::ostream& operator<<(std::ostream& out, const point_t& P) { out << "(" << P.x << ", " << P.y << ", " << P.z << ")"; return out; }

    double x, y, z;
};

struct vector_t : public point_t
{
public:
    vector_t() { x = 0.0; y = 0.0; z = 0.0; }
    vector_t(double i, double j, double k) { x = i; y = j; z = k; }
    vector_t(point_t p) { x = p.x; y = p.y; z = p.z; }

    bool similar_to(const vector_t& V) const { return (*this - V).L2_norm3D() < SIMILARITY_TOLERANCE; }

    const vector_t operator -(const vector_t V) const { return vector_t(x - V.x, y - V.y, z - V.z); }
    const vector_t operator +(const vector_t V) const { return vector_t(x + V.x, y + V.y, z + V.z); }
    const vector_t operator *(double a) const { return vector_t(a * x, a * y, a * z); }
    const vector_t operator /(double a) const {return a != 0 ? vector_t(x/a, y/a, z/a) : vector_t(x,y,z); }
    const vector_t operator -() const { return vector_t(-x, -y, -z); }
    void operator +=(const vector_t V) { x += V.x; y += V.y; z += V.z; }
    void operator -=(const vector_t V) { x -= V.x; y -= V.y; z -= V.z; }

    double operator *(vector_t u) const { return x * u.x + y * u.y + z * u.z; }
    const vector_t operator ^(vector_t V) const { vector_t U; U.x = y * V.z - z * V.y; U.y = z * V.x - x * V.z; U.z = x * V.y - y * V.x; return U; }

    bool operator !=(vector_t& V) { if (x != V.x || y != V.y || z != V.z) return true; return false; }
    bool operator ==(vector_t& V) { if (x == V.x && y == V.y && z == V.z) return true; return false; }

    bool is_zero() { if (x == 0.0 && y == 0.0 && z == 0.0) return true; return false; }
};

// Rotates v in the direction and angle of a -> b
static vector_t rotate(const vector_t& a, const vector_t& b, const vector_t& v)
{
    // Guard: If a and b are the same (0.0000001 is arbitrary theshold)
    // don't attempt because the cross product will be 0
    if ((a - b).L2_norm3D() < 0.0000001) return v;

    // Axis of rotation
    vector_t axis = a^b;
    axis.normalize();

    // Angle of rotation
    double angle = acos((a * b)/(a.L2_norm3D() * b.L2_norm3D()));

    // The rotated component
    vector_t f = axis^v;

    // cos(theta) of the original component
    // sin(theta) of the rotated component
    return v * cos(angle) + f * sin(angle);
}

static vector_t project_into_plane(vector_t vector, vector_t plane)
{
    vector_t projected_vector = vector;
    projected_vector.z = -(plane.x * projected_vector.x + plane.y * projected_vector.y) / plane.z;
    return projected_vector;
}

#endif // UTIL_MATH_H
