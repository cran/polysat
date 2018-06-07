// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// G
int G(int q, int n);
RcppExport SEXP _polysat_G(SEXP qSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(G(q, n));
    return rcpp_result_gen;
END_RCPP
}
// INDEXG
int INDEXG(IntegerVector ag1, int na1, int m2);
RcppExport SEXP _polysat_INDEXG(SEXP ag1SEXP, SEXP na1SEXP, SEXP m2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type ag1(ag1SEXP);
    Rcpp::traits::input_parameter< int >::type na1(na1SEXP);
    Rcpp::traits::input_parameter< int >::type m2(m2SEXP);
    rcpp_result_gen = Rcpp::wrap(INDEXG(ag1, na1, m2));
    return rcpp_result_gen;
END_RCPP
}
// GENLIST
IntegerMatrix GENLIST(int ng, int na1, int m2);
RcppExport SEXP _polysat_GENLIST(SEXP ngSEXP, SEXP na1SEXP, SEXP m2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ng(ngSEXP);
    Rcpp::traits::input_parameter< int >::type na1(na1SEXP);
    Rcpp::traits::input_parameter< int >::type m2(m2SEXP);
    rcpp_result_gen = Rcpp::wrap(GENLIST(ng, na1, m2));
    return rcpp_result_gen;
END_RCPP
}
// RANMUL
List RANMUL(int ng, int na1, IntegerMatrix ag, int m2);
RcppExport SEXP _polysat_RANMUL(SEXP ngSEXP, SEXP na1SEXP, SEXP agSEXP, SEXP m2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ng(ngSEXP);
    Rcpp::traits::input_parameter< int >::type na1(na1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ag(agSEXP);
    Rcpp::traits::input_parameter< int >::type m2(m2SEXP);
    rcpp_result_gen = Rcpp::wrap(RANMUL(ng, na1, ag, m2));
    return rcpp_result_gen;
END_RCPP
}
// SELFMAT
IntegerMatrix SELFMAT(int ng, int na1, IntegerMatrix ag, int m2);
RcppExport SEXP _polysat_SELFMAT(SEXP ngSEXP, SEXP na1SEXP, SEXP agSEXP, SEXP m2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ng(ngSEXP);
    Rcpp::traits::input_parameter< int >::type na1(na1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ag(agSEXP);
    Rcpp::traits::input_parameter< int >::type m2(m2SEXP);
    rcpp_result_gen = Rcpp::wrap(SELFMAT(ng, na1, ag, m2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_polysat_G", (DL_FUNC) &_polysat_G, 2},
    {"_polysat_INDEXG", (DL_FUNC) &_polysat_INDEXG, 3},
    {"_polysat_GENLIST", (DL_FUNC) &_polysat_GENLIST, 3},
    {"_polysat_RANMUL", (DL_FUNC) &_polysat_RANMUL, 4},
    {"_polysat_SELFMAT", (DL_FUNC) &_polysat_SELFMAT, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_polysat(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
