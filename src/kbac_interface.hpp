//: src:kbac_interface.hpp
// kbac class interface
// Copyright 2011 Gao Wang
#ifndef _KBAC_INTERFACE_HPP 
#define _KBAC_INTERFACE_HPP
#include "kbac.hpp"
static KbacTest* Ktest = NULL;
void set_up_kbac_test(int* nn, int* qq, double* aa, double* mafUpper, double* xdatIn, double* ydatIn, double* mafIn, int* xcol, int* ylen) {
    Ktest = new KbacTest(nn, qq, aa, mafUpper, xdatIn, ydatIn, mafIn, xcol, ylen);
    return;
}
void do_kbac_test(double* pvalue, int* sided) {
    Ktest->calcKbacP(pvalue, sided);
    return;
}
void clear_kbac_test() {
    delete Ktest;
    return;
}
#endif
