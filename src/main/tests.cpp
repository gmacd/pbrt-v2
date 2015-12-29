#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
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
    
    REQUIRE(boxB.pMin == Point(147, 117, 17));
    REQUIRE(boxB.pMax == Point(153, 123, 23));
}
