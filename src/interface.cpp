//: src:interface.cpp
// interface between cpp and R
// Copyright 2011 Gao Wang

#include <R.h>
#include "kbac_interface.hpp"

extern "C" {
    void R_set_up_kbac_test(int* nn, int* qq, double* aa, double* mafUpper, double* xdatIn, double* ydatIn, double* mafIn, int* xcol, int* ylen) {
        set_up_kbac_test(nn, qq, aa, mafUpper, xdatIn, ydatIn, mafIn, xcol, ylen);
    }

    void R_do_kbac_test(double* pvalue, int* sided) {
        do_kbac_test(pvalue, sided);
    }

    void R_clear_kbac_test() {
        clear_kbac_test();
    }
} // extern "C"
