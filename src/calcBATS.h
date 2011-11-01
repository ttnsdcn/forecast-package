#ifndef _forecast_CALCBATS
#define _forecast_CALCBATS

#include <RcppArmadillo.h>

RcppExport SEXP calcBATS(SEXP ys, SEXP yHats, SEXP wTransposes, SEXP Fs, SEXP xs, SEXP gs, SEXP es ) ;

RcppExport SEXP calcWTilda(SEXP wTildaTransposes, SEXP Ds) ;

RcppExport SEXP makeBATSWMatrix(SEXP smallPhi_s, SEXP sPeriods_s, SEXP arCoefs_s, SEXP maCoefs_s) ;

RcppExport SEXP makeBATSGMatrix(SEXP alpha_s, SEXP beta_s, SEXP gammaVector_s, SEXP seasonalPeriods_s, SEXP p_s, SEXP q_s) ;

#endif
