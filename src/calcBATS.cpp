#include "calcBATS.h"

using namespace Rcpp ;

SEXP calcBATS(SEXP ys, SEXP yHats, SEXP wTransposes, SEXP Fs, SEXP xs, SEXP gs, SEXP es ){
	BEGIN_RCPP
	

	NumericMatrix yr(ys);
	NumericMatrix yHatr(yHats);
	NumericMatrix wTransposer(wTransposes);
	NumericMatrix Fr(Fs);
	NumericMatrix xr(xs);
	NumericMatrix gr(gs);
	NumericMatrix er(es);

	int t;

	arma::mat y(yr.begin(), yr.nrow(), yr.ncol(), false);
	arma::mat yHat(yHatr.begin(), yHatr.nrow(), yHatr.ncol(), false);
	arma::mat wTranspose(wTransposer.begin(), wTransposer.nrow(), wTransposer.ncol(), false);
	arma::mat F(Fr.begin(), Fr.nrow(), Fr.ncol(), false);
	arma::mat x(xr.begin(), xr.nrow(), xr.ncol(), false);
	arma::mat g(gr.begin(), gr.nrow(), gr.ncol(), false);
	arma::mat e(er.begin(), er.nrow(), er.ncol(), false);


	for(t = 1; t < yr.ncol(); t++) {
		yHat.col(t) = wTranspose * x.col((t-1));
		e(0,t) = y(0, t) - yHat(0, t);
		x.col(t) = F * x.col((t-1)) + g * e(0,t);
	}

	return List::create(
			Named("y.hat") = yHat,
			Named("e") = e,
			Named("x") = x
	);

	END_RCPP
}

SEXP calcWTilda(SEXP wTildaTransposes, SEXP Ds) {
	BEGIN_RCPP

	NumericMatrix wTildaTransposer(wTildaTransposes);
	NumericMatrix Dr(Ds);

	int t;

	arma::mat wTildaTranspose(wTildaTransposer.begin(), wTildaTransposer.nrow(), wTildaTransposer.ncol(), false);
	arma::mat D(Dr.begin(), Dr.nrow(), Dr.ncol(), false);

	for(t = 1; t < wTildaTransposer.nrow(); t++) {
		wTildaTranspose.row(t) = wTildaTranspose.row((t-1)) * D;
	}

	return wTildaTransposer;

	END_RCPP
}


