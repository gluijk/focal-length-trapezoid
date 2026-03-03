#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector keystone_correct_cpp(NumericVector imgc,
								   NumericVector imgd,
								   NumericVector k) {
									   
	// This function writes on imgc the keystone corrected version of imgd
	IntegerVector dimd = imgd.attr("dim");
	IntegerVector dimc = imgc.attr("dim");

	int DIMYd = dimd[0];
	int DIMXd = dimd[1];
	int channels = dimd[2];

	int DIMYc = dimc[0];
	int DIMXc = dimc[1];

	for (int x = 0; x < DIMXc; x++) {
		for (int y = 0; y < DIMYc; y++) {

			double xd = x + 1.0;
			double yd = y + 1.0;

			// ---- undo.keystone inline ----
			double denom = k[6]*xd + k[7]*yd + 1.0;

			double xu = (k[0]*xd + k[1]*yd + k[2]) / denom;
			double yu = (k[3]*xd + k[4]*yd + k[5]) / denom;

			// ---- bilinear interpolation ----
			if (xu >= 1 && xu <= DIMXd-1 &&
				yu >= 1 && yu <= DIMYd-1) {

				int x1 = floor(xu) - 1;
				int y1 = floor(yu) - 1;
				int x2 = x1 + 1;
				int y2 = y1 + 1;

				double dx = xu - floor(xu);
				double dy = yu - floor(yu);

				for (int c = 0; c < channels; c++) {

					double Q11 = imgd[y1 + x1*DIMYd + c*DIMYd*DIMXd];
					double Q21 = imgd[y1 + x2*DIMYd + c*DIMYd*DIMXd];
					double Q12 = imgd[y2 + x1*DIMYd + c*DIMYd*DIMXd];
					double Q22 = imgd[y2 + x2*DIMYd + c*DIMYd*DIMXd];

					double val =
						Q11 * (1-dx)*(1-dy) +
						Q21 * dx*(1-dy) +
						Q12 * (1-dx)*dy +
						Q22 * dx*dy;

					imgc[y + x*DIMYc + c*DIMYc*DIMXc] = val;
				}
			}
		}
	}

	return imgc;
}



