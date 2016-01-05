#include "catch.hpp"
#include <iostream>
#include "geometry.h"

TEST_CASE("Scratch", "[scratch]")
{
    CHECK(__alignof__(Point) == 16);
    CHECK(sizeof(Point) == 16);
}
