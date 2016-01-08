#include "catch.hpp"
#include "geometry.h"

TEST_CASE("Basic vector ops", "[vector]")
{
    Vector v1(2, 2, 2);
    Vector v2(3, 3, 3);

    auto vRes = v1+v2;
    CHECK(vRes == Vector(5, 5, 5));

    vRes = v1-v2;
    CHECK(vRes == Vector(-1, -1, -1));

    vRes = -v1;
    CHECK(vRes == Vector(-2, -2, -2));
}

TEST_CASE("Length", "[vector]")
{
    Vector v(0, 0, 0);
    
    auto len = v.Length();
    auto lenSq = v.LengthSquared();

    CHECK(len == 0);
    CHECK(lenSq == 0);

    v = Vector(1, 1, 1);
    len = v.Length();
    lenSq = v.LengthSquared();

    CHECK(len == Approx(1.73205));
    CHECK(lenSq == Approx(3));

    v = Vector(-1, -1, -1);
    len = v.Length();
    lenSq = v.LengthSquared();

    CHECK(len == Approx(1.73205));
    CHECK(lenSq == Approx(3));
}

TEST_CASE("Basic BBox ops", "[bbox]")
{
    BBox b(Point(1, 2, 3), Point(4, 5, 6));
    CHECK(b.pMin.x == 1);
    CHECK(b.pMin.y == 2);
    CHECK(b.pMin.z == 3);
    CHECK(b.pMax.x == 4);
    CHECK(b.pMax.y == 5);
    CHECK(b.pMax.z == 6);
}

TEST_CASE("DotProduct", "[geom]")
{
    auto dp = Dot(Vector(1, 2, 3), Vector(0, 0, 0));
    CHECK(dp == 0);

    dp = Dot(Vector(1, 2, 3), Vector(10, 20, 30));
    CHECK(dp == 140);

    dp = Dot(Vector(1, 2, 3), Vector(-10, -20, -30));
    CHECK(dp == -140);

    dp = AbsDot(Vector(1, 2, 3), Vector(-10, -20, -30));
    CHECK(dp == 140);
}
