/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#pragma once


#include "pbrt.h"


class
__attribute__ ((aligned(16)))
Vector
{
public:
    union
    {
        struct { float x, y, z; };
        float4_t _vec;
    };


    Vector()
    {
        x = y = z = 0.f;
    }
    
    Vector(float xx, float yy, float zz):
        x(xx), y(yy), z(zz)
    {
        Assert(!HasNaNs());
    }
    
    Vector(const Vector &v)
    {
        Assert(!v.HasNaNs());
        x = v.x;  y = v.y;  z = v.z;
    }

    Vector(float4_t v):
        _vec(v)
    {
    }

    
    explicit Vector(const Point &p);
    explicit Vector(const Normal &n);

    
    Vector &operator=(const Vector &v)
    {
        Assert(!v.HasNaNs());
        x = v.x;  y = v.y;  z = v.z;
        return *this;
    }

    bool operator==(const Vector &v) const
    {
        return x == v.x && y == v.y && z == v.z;
    }
    
    bool operator!=(const Vector &v) const
    {
        return x != v.x || y != v.y || z != v.z;
    }
    
    Vector operator+(const Vector &v) const
    {
        Assert(!v.HasNaNs());
#ifdef USE_SIMD_AVX
        return Vector(_mm_add_ps(_vec, v._vec));
#else
        return Vector(x + v.x, y + v.y, z + v.z);
#endif
    }
    
    Vector& operator+=(const Vector &v)
    {
        Assert(!v.HasNaNs());
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    
    Vector operator-(const Vector &v) const
    {
        Assert(!v.HasNaNs());
        return Vector(x - v.x, y - v.y, z - v.z);
    }
    
    Vector& operator-=(const Vector &v)
    {
        Assert(!v.HasNaNs());
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }
    
    Vector operator*(float f) const
    {
        return Vector(f*x, f*y, f*z);
    }
    
    Vector &operator*=(float f)
    {
        Assert(!isnan(f));
        x *= f; y *= f; z *= f;
        return *this;
    }
    
    Vector operator/(float f) const
    {
        Assert(f != 0);
        float inv = 1.f / f;
        return Vector(x * inv, y * inv, z * inv);
    }
    
    Vector &operator/=(float f)
    {
        Assert(f != 0);
        float inv = 1.f / f;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }
    
    Vector operator-() const
    {
        return Vector(-x, -y, -z);
    }
    
    float operator[](int i) const
    {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }
    
    float &operator[](int i)
    {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }
    

    bool HasNaNs() const
    {
        return isnan(x) || isnan(y) || isnan(z);
    }
        
    float LengthSquared() const
    {
#ifdef USE_SIMD_AVX
        return _mm_cvtss_f32(_mm_dp_ps(_vec, _vec, 0x71));
#else
        return x*x + y*y + z*z;
#endif
    }
    
    float Length() const
    {
#ifdef USE_SIMD_AVX
        return _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(_vec, _vec, 0x71)));
#else
        return sqrtf(LengthSquared());
#endif
    }
};


class Point3
{
public:
    float x, y, z;
    
    Point3()
    {
        x = y = z = 0.f;
    }
    
    Point3(float xx, float yy, float zz):
        x(xx), y(yy), z(zz)
    {
    }
    
    Point3(const Point3 &p)
    {
        x = p.x;  y = p.y;  z = p.z;
    }
    
    explicit Point3(const Point &p);
};


class
__attribute__ ((aligned(16)))
Point
{
public:
    union
    {
        struct { float x, y, z; };
        float4_t _vec;
    };

    
    Point()
    {
        x = y = z = 0.f;
    }
    
    Point(float xx, float yy, float zz):
        x(xx), y(yy), z(zz)
    {
        Assert(!HasNaNs());
    }
    
    Point(const Point &p)
    {
        Assert(!p.HasNaNs());
        x = p.x;  y = p.y;  z = p.z;
    }
    
    explicit Point(const Point3 &p)
    {
        x = p.x;  y = p.y;  z = p.z;
        Assert(!HasNaNs());
    }

    
    Point &operator=(const Point &p)
    {
        Assert(!p.HasNaNs());
        x = p.x;  y = p.y;  z = p.z;
        return *this;
    }
    

    bool operator==(const Point &p) const
    {
        return x == p.x && y == p.y && z == p.z;
    }

    bool operator!=(const Point &p) const
    {
        return x != p.x || y != p.y || z != p.z;
    }


    Point operator+(const Vector &v) const
    {
        Assert(!v.HasNaNs());
        return Point(x + v.x, y + v.y, z + v.z);
    }
    
    Point &operator+=(const Vector &v)
    {
        Assert(!v.HasNaNs());
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    
    Vector operator-(const Point &p) const
    {
        Assert(!p.HasNaNs());
        return Vector(x - p.x, y - p.y, z - p.z);
    }
    
    Point operator-(const Vector &v) const
    {
        Assert(!v.HasNaNs());
        return Point(x - v.x, y - v.y, z - v.z);
    }
    
    Point &operator-=(const Vector &v)
    {
        Assert(!v.HasNaNs());
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }

    Point &operator+=(const Point &p)
    {
        Assert(!p.HasNaNs());
        x += p.x; y += p.y; z += p.z;
        return *this;
    }

    Point operator+(const Point &p) const
    {
        Assert(!p.HasNaNs());
        return Point(x + p.x, y + p.y, z + p.z);
    }

    Point operator* (float f) const
    {
        return Point(f*x, f*y, f*z);
    }

    Point &operator*=(float f)
    {
        x *= f; y *= f; z *= f;
        return *this;
    }

    Point operator/ (float f) const
    {
        float inv = 1.f/f;
        return Point(inv*x, inv*y, inv*z);
    }

    Point &operator/=(float f)
    {
        float inv = 1.f/f;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }

    float operator[](int i) const
    {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }
    
    float &operator[](int i)
    {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }

    
    bool HasNaNs() const
    {
        return isnan(x) || isnan(y) || isnan(z);
    }
};


class
__attribute__ ((aligned(16)))
Normal
{
public:
    union
    {
        struct { float x, y, z; };
        float4_t _vec;
    };

    
    Normal()
    {
        x = y = z = 0.f;
    }

    Normal(float xx, float yy, float zz):
        x(xx), y(yy), z(zz)
    {
        Assert(!HasNaNs());
    }

    Normal(const Normal &n)
    {
        Assert(!n.HasNaNs());
        x = n.x; y = n.y; z = n.z;
    }

    explicit Normal(const Vector &v):
        x(v.x), y(v.y), z(v.z)
    {
        Assert(!v.HasNaNs());
    }

    
    Normal &operator=(const Normal &n)
    {
        Assert(!n.HasNaNs());
        x = n.x; y = n.y; z = n.z;
        return *this;
    }
        
    bool operator==(const Normal &n) const
    {
        return x == n.x && y == n.y && z == n.z;
    }

    bool operator!=(const Normal &n) const
    {
        return x != n.x || y != n.y || z != n.z;
    }

    Normal operator-() const
    {
        return Normal(-x, -y, -z);
    }

    Normal operator+ (const Normal &n) const
    {
        Assert(!n.HasNaNs());
        return Normal(x + n.x, y + n.y, z + n.z);
    }
    
    Normal& operator+=(const Normal &n)
    {
        Assert(!n.HasNaNs());
        x += n.x; y += n.y; z += n.z;
        return *this;
    }

    Normal operator- (const Normal &n) const
    {
        Assert(!n.HasNaNs());
        return Normal(x - n.x, y - n.y, z - n.z);
    }
    
    Normal& operator-=(const Normal &n)
    {
        Assert(!n.HasNaNs());
        x -= n.x; y -= n.y; z -= n.z;
        return *this;
    }

    Normal operator*(float f) const
    {
        return Normal(f*x, f*y, f*z);
    }
    
    Normal &operator*=(float f)
    {
        x *= f; y *= f; z *= f;
        return *this;
    }

    Normal operator/(float f) const
    {
        Assert(f != 0);
        float inv = 1.f/f;
        return Normal(x * inv, y * inv, z * inv);
    }
    
    Normal &operator/=(float f)
    {
        Assert(f != 0);
        float inv = 1.f/f;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }

    float operator[](int i) const
    {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }
    
    float &operator[](int i)
    {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }

    
    bool HasNaNs() const
    {
        return isnan(x) || isnan(y) || isnan(z);
    }

    float LengthSquared() const
    {
        return x*x + y*y + z*z;
    }

    float Length() const
    {
        return sqrtf(LengthSquared());
    }
};


class Ray
{
public:
    Point o;
    Vector d;
    mutable float mint, maxt;
    float time;
    int depth;


    Ray():
        mint(0.f), maxt(INFINITY), time(0.f), depth(0)
    {
    }
    
    Ray(const Point &origin, const Vector &direction,
        float start, float end = INFINITY, float t = 0.f, int d = 0):
        o(origin), d(direction), mint(start), maxt(end), time(t), depth(d)
    {
    }

    Ray(const Point &origin, const Vector &direction, const Ray &parent,
        float start, float end = INFINITY):
        o(origin), d(direction), mint(start), maxt(end),
        time(parent.time), depth(parent.depth+1)
    {
    }

    
    Point operator()(float t) const
    {
        return o + d * t;
    }

    
    bool HasNaNs() const
    {
        return (o.HasNaNs() || d.HasNaNs() ||
                isnan(mint) || isnan(maxt));
    }
};


class RayDifferential :
    public Ray
{
public:
    bool hasDifferentials;
    Point rxOrigin, ryOrigin;
    Vector rxDirection, ryDirection;

    
    RayDifferential()
    {
        hasDifferentials = false;
    }
    
    RayDifferential(const Point &org, const Vector &dir, float start,
                    float end = INFINITY, float t = 0.f, int d = 0):
        Ray(org, dir, start, end, t, d)
    {
        hasDifferentials = false;
    }

    RayDifferential(const Point &org, const Vector &dir, const Ray &parent,
                    float start, float end = INFINITY):
        Ray(org, dir, start, end, parent.time, parent.depth+1)
    {
        hasDifferentials = false;
    }

    explicit RayDifferential(const Ray &ray):
        Ray(ray)
    {
        hasDifferentials = false;
    }

    
    bool HasNaNs() const
    {
        return Ray::HasNaNs() ||
            (hasDifferentials && (rxOrigin.HasNaNs() || ryOrigin.HasNaNs() ||
                                  rxDirection.HasNaNs() || ryDirection.HasNaNs()));
    }

    void ScaleDifferentials(float s)
    {
        rxOrigin = o + (rxOrigin - o) * s;
        ryOrigin = o + (ryOrigin - o) * s;
        rxDirection = d + (rxDirection - d) * s;
        ryDirection = d + (ryDirection - d) * s;
    }
};


class BBox
{
public:
    Point pMin, pMax;

    
    BBox()
    {
        pMin = Point( INFINITY,  INFINITY,  INFINITY);
        pMax = Point(-INFINITY, -INFINITY, -INFINITY);
    }

    BBox(const Point &p):
        pMin(p), pMax(p)
    {
    }

    BBox(const Point &p1, const Point &p2)
    {
        pMin = Point(min(p1.x, p2.x), min(p1.y, p2.y), min(p1.z, p2.z));
        pMax = Point(max(p1.x, p2.x), max(p1.y, p2.y), max(p1.z, p2.z));
    }

    
    bool operator==(const BBox &b) const
    {
        return b.pMin == pMin && b.pMax == pMax;
    }

    bool operator!=(const BBox &b) const
    {
        return b.pMin != pMin || b.pMax != pMax;
    }

    const Point &operator[](int i) const;
    
    Point &operator[](int i);

    
    friend BBox Union(const BBox &b, const Point &p);
    friend BBox Union(const BBox &b, const BBox &b2);
    
    bool Overlaps(const BBox &b) const
    {
        bool x = (pMax.x >= b.pMin.x) && (pMin.x <= b.pMax.x);
        bool y = (pMax.y >= b.pMin.y) && (pMin.y <= b.pMax.y);
        bool z = (pMax.z >= b.pMin.z) && (pMin.z <= b.pMax.z);
        return (x && y && z);
    }

    bool Inside(const Point &pt) const
    {
        return (pt.x >= pMin.x && pt.x <= pMax.x &&
                pt.y >= pMin.y && pt.y <= pMax.y &&
                pt.z >= pMin.z && pt.z <= pMax.z);
    }

    void Expand(float delta)
    {
        pMin -= Vector(delta, delta, delta);
        pMax += Vector(delta, delta, delta);
    }

    float SurfaceArea() const
    {
        Vector d = pMax - pMin;
        return 2.f * (d.x * d.y + d.x * d.z + d.y * d.z);
    }

    float Volume() const
    {
        Vector d = pMax - pMin;
        return d.x * d.y * d.z;
    }

    int MaximumExtent() const
    {
        Vector diag = pMax - pMin;
        if (diag.x > diag.y && diag.x > diag.z)
            return 0;
        else if (diag.y > diag.z)
            return 1;
        else
            return 2;
    }
        
    Point Lerp(float tx, float ty, float tz) const
    {
        return Point(::Lerp(tx, pMin.x, pMax.x), ::Lerp(ty, pMin.y, pMax.y),
                     ::Lerp(tz, pMin.z, pMax.z));
    }

    Vector Offset(const Point &p) const
    {
        return Vector((p.x - pMin.x) / (pMax.x - pMin.x),
                      (p.y - pMin.y) / (pMax.y - pMin.y),
                      (p.z - pMin.z) / (pMax.z - pMin.z));
    }

    void BoundingSphere(Point *c, float *rad) const;
    
    bool IntersectP(const Ray &ray, float *hitt0 = NULL, float *hitt1 = NULL) const;
};


inline Vector::Vector(const Point &p):
    x(p.x), y(p.y), z(p.z)
{
    Assert(!HasNaNs());
}


inline Vector::Vector(const Normal &n):
    x(n.x), y(n.y), z(n.z)
{
    Assert(!n.HasNaNs());
}


inline Point3::Point3(const Point &p)
{
    x = p.x;  y = p.y;  z = p.z;
}


inline const Point &BBox::operator[](int i) const
{
    Assert(i == 0 || i == 1);
    return (&pMin)[i];
}


inline Point &BBox::operator[](int i)
{
    Assert(i == 0 || i == 1);
    return (&pMin)[i];
}


inline Vector operator*(float f, const Vector &v)
{
    return v*f;
}

inline float Dot(const Vector &v1, const Vector &v2)
{
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
#ifdef USE_SIMD_AVX
    return _mm_cvtss_f32(_mm_dp_ps(v1._vec, v2._vec, 0x71));
#else
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
#endif
}


inline float AbsDot(const Vector &v1, const Vector &v2)
{
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
    return fabsf(Dot(v1, v2));
}

inline __m128 Cross(const __m128* v1, const __m128* v2)
{
    return _mm_sub_ps(
        _mm_mul_ps(_mm_shuffle_ps(*v1, *v1, _MM_SHUFFLE(3, 0, 2, 1)),
                   _mm_shuffle_ps(*v2, *v2, _MM_SHUFFLE(3, 1, 0, 2))),
        _mm_mul_ps(_mm_shuffle_ps(*v1, *v1, _MM_SHUFFLE(3, 1, 0, 2)),
                   _mm_shuffle_ps(*v2, *v2, _MM_SHUFFLE(3, 0, 2, 1))));
}

inline Vector Cross(const Vector &v1, const Vector &v2)
{
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
#ifdef USE_SIMD_AVX
    return Cross(&v1._vec, &v2._vec);
#else
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector((v1y * v2z) - (v1z * v2y),
                  (v1z * v2x) - (v1x * v2z),
                  (v1x * v2y) - (v1y * v2x));
#endif
}


inline Vector Cross(const Vector &v1, const Normal &v2)
{
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
#ifdef USE_SIMD_AVX
    return Cross(&v1._vec, &v2._vec);
#else
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector((v1y * v2z) - (v1z * v2y),
                  (v1z * v2x) - (v1x * v2z),
                  (v1x * v2y) - (v1y * v2x));
#endif
}


inline Vector Cross(const Normal &v1, const Vector &v2)
{
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
#ifdef USE_SIMD_AVX
    return Cross(&v1._vec, &v2._vec);
#else
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector((v1y * v2z) - (v1z * v2y),
                  (v1z * v2x) - (v1x * v2z),
                  (v1x * v2y) - (v1y * v2x));
#endif
}


inline Vector Normalize(const Vector &v)
{
    return v / v.Length();
}

inline void CoordinateSystem(const Vector &v1, Vector *v2, Vector *v3)
{
    if (fabsf(v1.x) > fabsf(v1.y)) {
        float invLen = 1.f / sqrtf(v1.x*v1.x + v1.z*v1.z);
        *v2 = Vector(-v1.z * invLen, 0.f, v1.x * invLen);
    }
    else {
        float invLen = 1.f / sqrtf(v1.y*v1.y + v1.z*v1.z);
        *v2 = Vector(0.f, v1.z * invLen, -v1.y * invLen);
    }
    *v3 = Cross(v1, *v2);
}


inline float Distance(const Point &p1, const Point &p2)
{
    return (p1 - p2).Length();
}


inline float DistanceSquared(const Point &p1, const Point &p2)
{
    return (p1 - p2).LengthSquared();
}


inline Point operator*(float f, const Point &p)
{
    Assert(!p.HasNaNs());
    return p*f;
}


inline Normal operator*(float f, const Normal &n)
{
    return Normal(f*n.x, f*n.y, f*n.z);
}


inline Normal Normalize(const Normal &n)
{
    return n / n.Length();
}


inline float Dot(const Normal &n1, const Vector &v2)
{
    Assert(!n1.HasNaNs() && !v2.HasNaNs());
#ifdef USE_SIMD_AVX
    return _mm_cvtss_f32(_mm_dp_ps(n1._vec, v2._vec, 0x71));
    //return n1.x * v2.x + n1.y * v2.y + n1.z * v2.z;
#else
    return n1.x * v2.x + n1.y * v2.y + n1.z * v2.z;
#endif
}


inline float Dot(const Vector &v1, const Normal &n2)
{
    Assert(!v1.HasNaNs() && !n2.HasNaNs());
#ifdef USE_SIMD_AVX
    return _mm_cvtss_f32(_mm_dp_ps(v1._vec, n2._vec, 0x71));
    //return v1.x * n2.x + v1.y * n2.y + v1.z * n2.z;
#else
    return v1.x * n2.x + v1.y * n2.y + v1.z * n2.z;
#endif
}


inline float Dot(const Normal &n1, const Normal &n2)
{
    Assert(!n1.HasNaNs() && !n2.HasNaNs());
#ifdef USE_SIMD_AVX
    return _mm_cvtss_f32(_mm_dp_ps(n1._vec, n2._vec, 0x71));
    //return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
#else
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
#endif
}


inline float AbsDot(const Normal &n1, const Vector &v2)
{
    Assert(!n1.HasNaNs() && !v2.HasNaNs());
    return fabsf(Dot(n1, v2));
}


inline float AbsDot(const Vector &v1, const Normal &n2)
{
    Assert(!v1.HasNaNs() && !n2.HasNaNs());
    return fabsf(Dot(v1, n2));
}


inline float AbsDot(const Normal &n1, const Normal &n2)
{
    Assert(!n1.HasNaNs() && !n2.HasNaNs());
    return fabsf(Dot(n1, n2));
}


inline Normal Faceforward(const Normal &n, const Vector &v)
{
    return (Dot(n, v) < 0.f) ? -n : n;
}


inline Normal Faceforward(const Normal &n, const Normal &n2)
{
    return (Dot(n, n2) < 0.f) ? -n : n;
}



inline Vector Faceforward(const Vector &v, const Vector &v2)
{
    return (Dot(v, v2) < 0.f) ? -v : v;
}


inline Vector Faceforward(const Vector &v, const Normal &n2)
{
    return (Dot(v, n2) < 0.f) ? -v : v;
}


inline Vector SphericalDirection(float sintheta,
                                 float costheta, float phi)
{
    return Vector(sintheta * cosf(phi),
                  sintheta * sinf(phi),
                  costheta);
}


inline Vector SphericalDirection(float sintheta, float costheta,
                                 float phi, const Vector &x,
                                 const Vector &y, const Vector &z)
{
    return sintheta * cosf(phi) * x +
           sintheta * sinf(phi) * y + costheta * z;
}


inline float SphericalTheta(const Vector &v)
{
    return acosf(Clamp(v.z, -1.f, 1.f));
}


inline float SphericalPhi(const Vector &v)
{
    float p = atan2f(v.y, v.x);
    return (p < 0.f) ? p + 2.f*M_PI : p;
}
