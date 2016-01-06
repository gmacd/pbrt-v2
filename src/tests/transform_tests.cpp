#include "catch.hpp"

#include "transform.h"
#include "geometry.h"


TEST_CASE("Transform identity", "[transforms]")
{
    Transform transform;

    BBox boxA(Point(-3, -3, -3), Point(3, 3, 3));
    auto boxB = transform(boxA);
    
    REQUIRE(boxB.pMin == boxA.pMin);
    REQUIRE(boxB.pMax == boxA.pMax);
}

TEST_CASE("Transform bounding box", "[transforms]")
{
    // Example taken from killeroo scene
    auto transform = Translate(Vector(150, 120, 20));

    BBox boxA(Point(-3, -3, -3), Point(3, 3, 3));
    auto boxB = transform(boxA);
    
    REQUIRE(boxB.pMin.x == 147);
    REQUIRE(boxB.pMin.y == 117);
    REQUIRE(boxB.pMin.z == 17);
    REQUIRE(boxB.pMax.x == 153);
    REQUIRE(boxB.pMax.y == 123);
    REQUIRE(boxB.pMax.z == 23);
}
