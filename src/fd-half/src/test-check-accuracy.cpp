#include "Common.h"

namespace {

    void checkAccuracy(BmpReal k, BmpReal left, BmpReal right) {
        setPreciseOutput();
        std::cout << "k = " << k << ", d = " << abs(left/right - 1) << std::endl;
    }
}

TEST_CASE("check_integer_index") {
    SECTION("0") {
        // TODO: fix
        BmpReal I_fcs;
        BmpReal I_prec = log(2);
        //checkAccuracy(0, I_fcs, I_prec);
    }
    SECTION("1") {
        BmpReal I_fcs = 0.82246703342411309;
        BmpReal I_prec = pi()*pi() / 12;
        checkAccuracy(1, I_fcs, I_prec);
    }
    SECTION("3") {
        BmpReal I_fcs = 5.6821969769834748;
        BmpReal I_prec = 7 * pow(pi(), 4) / 120;
        checkAccuracy(3, I_fcs, I_prec);
    }
}

TEST_CASE("check_fcs_and_quad") {
    SECTION("m12") {
        BmpReal I_fcs = 1.0721549299401913;
        BmpReal I_quad = 1.0721549299401913;
        BmpReal k = -0.5;
        checkAccuracy(k, I_fcs, I_quad);
    }
    SECTION("12") {
        BmpReal I_fcs = 0.67809389515310081;
        BmpReal I_quad = 0.6780938951531011;
        BmpReal k = 0.5;
        checkAccuracy(k, I_fcs, I_quad);
    }
    SECTION("32") {
        BmpReal I_fcs =  1.1528038370883613;
        BmpReal I_quad = 1.1528038370883615;
        BmpReal k = 1.5;
        checkAccuracy(k, I_fcs, I_quad);
    }
    SECTION("52") {
        BmpReal I_fcs =  3.0825860828374179;
        BmpReal I_quad = 3.0825860828374188;
        BmpReal k = 2.5;
        checkAccuracy(k, I_fcs, I_quad);
    }
    SECTION("72") {
        BmpReal I_fcs =  11.183716751693318;
        BmpReal I_quad = 11.18371675169332;
        BmpReal k = 3.5;
        checkAccuracy(k, I_fcs, I_quad);
    }
}
