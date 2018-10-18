#include "Common.h"

TEST_CASE("bernulli") {
    SECTION("odd") {
        CHECK(0 == bernulli(1));
        CHECK(0 == bernulli(3));
        CHECK(0 == bernulli(5));
        CHECK(0 == bernulli(7));
        CHECK(0 == bernulli(9));
        CHECK(0 == bernulli(11));
        CHECK(0 == bernulli(13));
        CHECK(0 == bernulli(15));
        CHECK(0 == bernulli(17));
        CHECK(0 == bernulli(19));
        CHECK(0 == bernulli(21));
        CHECK(0 == bernulli(23));
    }
    SECTION("even") {
        CHECK(_1 / 6 == bernulli(2));
        CHECK(-_1 / 30 == bernulli(4));
        CHECK(_1 / 42 == bernulli(6));
        CHECK(-_1 / 30 == bernulli(8));
        CHECK(_5 / 66 == bernulli(10));
        CHECK(-_691 / 2730 == bernulli(12));
        CHECK(_7 / 6 == bernulli(14));
        CHECK(-_3617 / 510 == bernulli(16));
        CHECK(_43867 / 798 == bernulli(18));
        CHECK(-_174611 / 330 == bernulli(20));
        CHECK(_854513 / 138 == bernulli(22));
        CHECK(-_236364091 / 2730 == bernulli(24));
    }
}

TEST_CASE("dzeta") {
    SECTION("odd") {
        CHECK(0 == dzetaFunction(1));
        CHECK(0 == dzetaFunction(3));
        CHECK(0 == dzetaFunction(5));
        CHECK(0 == dzetaFunction(7));
        CHECK(0 == dzetaFunction(9));
        CHECK(0 == dzetaFunction(11));
        CHECK(0 == dzetaFunction(13));
        CHECK(0 == dzetaFunction(15));
        CHECK(0 == dzetaFunction(17));
        CHECK(0 == dzetaFunction(19));
        CHECK(0 == dzetaFunction(21));
        CHECK(0 == dzetaFunction(23));
    }
    SECTION("even") {
        CHECK(_1 / 6 == dzetaFunction(2) / pow(pi(), 2));
        CHECK(_1 / 90 == dzetaFunction(4) / pow(pi(), 4));
        CHECK(_1 / 945 == dzetaFunction(6) / pow(pi(), 6));
        CHECK(_1 / 9450 == dzetaFunction(8) / pow(pi(), 8));
        CHECK(_1 / 93550 == dzetaFunction(10) / pow(pi(), 10));
        CHECK(_691 / 638512875 == dzetaFunction(12) / pow(pi(), 12));
        CHECK(_2 / 18243225 == dzetaFunction(14) / pow(pi(), 14));
        CHECK(_3617 / 325641566250 == dzetaFunction(16) / pow(pi(), 16));
        CHECK(_43867 / 38979295480125 == dzetaFunction(18) / pow(pi(), 18));
        CHECK(_174611 / 1531329465290625 == dzetaFunction(20) / pow(pi(), 20));
        CHECK(_77683 / 13447856940643124 == dzetaFunction(22) / pow(pi(), 22));
        CHECK(_236364091 / BmpReal(201919571963756511232.0) == dzetaFunction(24) / pow(pi(), 24));
    }
}
