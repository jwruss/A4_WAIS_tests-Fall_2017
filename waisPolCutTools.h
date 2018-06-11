#include <cmath>
#include "TMath.h"
#include "TString.h"
#include "FFTWComplex.h"
#include "FFTtools.h"


/*  TGraph::Integral() returns the area within a closed polygon of a TGraph after bin reshuffling,
 *  which is a different operation than integrating the function corresponding to TGraph::GetY() over
 *  the domain corresponding to TGraph::GetX(). Using trapezoidal rule and accounting for potential uneven
 *  sampling in the TGraph, we explicitly calculate the integral here.
 */
double getGraphIntegral(const TGraph * grInPtr) {

	double * XIn = grInPtr -> GetX();
	double * YIn = grInPtr -> GetY();

	double integral = 0;
	for (int i = 1; i < grInPtr -> GetN(); ++i) integral += (YIn[i] + YIn[i - 1]) * (XIn[i] - XIn[i - 1]) / 2;

	return integral;
}


/*  The normalization done in FFTtools::getNormalisedCorrelationGraph() appears to be wrong.
 *  We correct it here. Dropping "normalization" from the name of getCorrGraph() comes from the fact
 *  that what FFTtools::getCorrelationGraph() actually attempts to produce is the cross-covariance
 *  graph. Its normalization is what is usually referred to as the correlation, or cross-correlation.
 */
TGraph getCorrGraph(const TGraph * grIn1Ptr, const TGraph * grIn2Ptr) {

	TGraph * grPtr = FFTtools::getCorrelationGraph(grIn1Ptr, grIn2Ptr);
	TGraph grOut(* grPtr);

	delete grPtr;

	int NOut = grOut.GetN();
	double norm = 1 / (sqrt(grIn1Ptr -> GetN() * grIn2Ptr -> GetN()) * grIn1Ptr -> GetRMS(2) * grIn2Ptr -> GetRMS(2));
	norm *= NOut / 4;  //  Additional factor from Ryan's silly normalization convention.
	double * YOut = grOut.GetY();
	for (int i = 0; i < NOut; ++i) YOut[i] *= norm;

	return grOut;
}


/*  To use interpolation of TGraphs used in getCorrGraph() instead.
 */
TGraph getInterpolatedCorrGraph(const TGraph * grIn1Ptr, const TGraph * grIn2Ptr, double deltaT = 0.1) {

	TGraph * gr1Ptr = FFTtools::getInterpolatedGraph(grIn1Ptr, deltaT);
	TGraph * gr2Ptr = FFTtools::getInterpolatedGraph(grIn2Ptr, deltaT);
	TGraph grOut = getCorrGraph(gr1Ptr, gr2Ptr);

	delete gr1Ptr, delete gr2Ptr;

	return grOut;
}


/*  Create a truncation of input TGraph grInPtr to a ns-wide window centered about the pulse's zero crossing between maximum and minimum values.
 *  The input winFullWidth is the ns-width over the full window; it's automatically adjusted when outside the scope of grInPtr's range, such that
 *  the length of winFullWidth retains the same parity over the entire calculation.
 */
TGraph getWindowedGraph(const TGraph * grInPtr, double winFullWidth = 5) {

	int NIn = grInPtr -> GetN();
	double * XIn = grInPtr -> GetX();
	double * YIn = grInPtr -> GetY();

	//  Find index "zeroInd" corresponding to zero crossing.
	int maxInd = TMath::LocMax(NIn, YIn);
	int minInd = TMath::LocMin(NIn, YIn);
	int startInd = TMath::Min(minInd, maxInd);
	int endInd = TMath::Max(minInd, maxInd);
	int zeroInd = startInd;
	double zeroYSq = YIn[zeroInd] * YIn[zeroInd];
	for (int i = startInd + 1; i <= endInd; ++i) {

		double YSq = YIn[i] * YIn[i];
		if (YSq < zeroYSq) {

			zeroInd = i;
			zeroYSq = YSq;
		}
	}

	//  Find length of new TGraph from assuming same time step size.
	int NOut = (NIn - 1) * winFullWidth / (XIn[NIn - 1] - XIn[0]) + 1;
	//  What follows is for in case the window width is too wide.
	if (NOut > NIn) {
		
		int NOld = NOut;
		cout << "Window length is larger than input TGraph length? Let's pretend that wasn't requested." << endl;
		NOut = (NOld - NIn - 1) % 2 ? NIn : NIn - 1;
		cout << "Window length adjusted from " << NOld << " to " << NOut << " bins." << endl;
	}
	//  For setting the window's start and end indices, and also adjusting the window size
	//  when necessary, based on the assumptions
	//  zeroInd - NOut / 2 >= 0, zeroInd + NOut / 2 <= NIn - 1.
	int winStartInd, winEndInd;
	if (zeroInd <= NIn - 1 - zeroInd) {

		winStartInd = zeroInd - (NOut - 1) / 2;
		winEndInd = zeroInd + NOut / 2;
		if (winStartInd < 0) {

			cout << "Window length extends below smallest input TGraph bin." << endl;
			int NOld = NOut;
			winStartInd = 0;
			winEndInd = 2 * zeroInd + (NOld - 1) % 2;
			NOut = winEndInd - winStartInd + 1;
			cout << "Window length adjusted from " << NOld << " to " << NOut << " bins.\n" << endl;
		}
	} else {

		winStartInd = zeroInd - NOut / 2;
		winEndInd = zeroInd + (NOut - 1) / 2;
		if (winEndInd > NIn - 1) {

			cout << "Window length extends above largest input TGraph bin." << endl;
			int NOld = NOut;
			winStartInd = 2 * (NIn - 1) - zeroInd - (NOld - 1) % 2;
			winEndInd = NIn - 1;
			NOut = winEndInd - winStartInd + 1;
			cout << "Window length adjusted from " << NOld << " to " << NOut << " bins.\n" << endl;
		}
	}

	TGraph grOut(NOut);

	double * XOut = grOut.GetX();
	std::copy(XIn + winStartInd, XIn + winEndInd + 1, XOut);  //  The inputs are addresses, plus index (index 0 is droppable).

	double * YOut = grOut.GetY();
	std::copy(YIn + winStartInd, YIn + winEndInd + 1, YOut);

	grOut.SetTitle(grInPtr -> GetTitle());
	grOut.GetXaxis() -> SetTitle(grInPtr -> GetXaxis() -> GetTitle());
	grOut.GetYaxis() -> SetTitle(grInPtr -> GetYaxis() -> GetTitle());

	return grOut;
}


/*  To use interpolation of TGraph used in getWindowedGraph() instead.
 */
TGraph getWindowedInterpolatedGraph(const TGraph * grInPtr, double deltaT = 0.1, double winFullWidth = 5) {

	TGraph * grInterpPtr = FFTtools::getInterpolatedGraph(grInPtr, deltaT);
	TGraph grOut = getWindowedGraph(grInterpPtr, winFullWidth);

	delete grInterpPtr;

	return grOut;
}


/*  Get TGraph which is an affine transformation of input TGraph.
 *  Setting "scale = -1, offset = 0" flips polarity.
 */
TGraph getAffineTransformedGraph(const TGraph * grInPtr, double scale = 1, double shift = 0) {

	TGraph grOut(* grInPtr);
	grOut.SetTitle();
	grOut.GetYaxis() -> SetTitle();

	double * YOut = grOut.GetY();
	for (int i = 0; i < grOut.GetN(); ++i) YOut[i] = scale * YOut[i] + shift;

	return grOut;
}


/*  Take an input TGraph and generate another TGraph corresponding to its numeric central derivative.
 *  Technically, the input TGraph doesn't have to be uniformly sampled.
 */
TGraph getDerivGraph(const TGraph * grInPtr) {

	int NIn = grInPtr -> GetN();
	double * XIn = grInPtr -> GetX();
	double * YIn = grInPtr -> GetY();

	TGraph grOut(* grInPtr);
	grOut.SetTitle();
	grOut.GetYaxis() -> SetTitle();

	//  Construct array YOut corresponding to derivative of YIn.
	double * YOut = grOut.GetY();
	YOut[0] = (YIn[1] - YIn[0]) / (XIn[1] - XIn[0]);
	for (int i = 1; i < NIn - 1; ++i) YOut[i] = (YIn[i + 1] - YIn[i - 1]) / (XIn[i + 1] - XIn[i - 1]);
//	//  When assuming trapezoid rule antiderivative, with lower bound of 0.
//	for (int i = 1; i < NIn - 1; ++i) {
//
//		YOut[i] = (YIn[i + 1] - YIn[i]) / (XIn[i + 1] - XIn[i]);
//		YOut[i] += (YIn[i] - YIn[i - 1]) / (XIn[i] - XIn[i - 1]);
//		YOut[i] /= 2;
//	}
	YOut[NIn - 1] = (YIn[NIn - 1] - YIn[NIn - 2]) / (XIn[NIn - 1] - XIn[NIn - 2]);

//	TF1 fIn("fIn", [&](double * x, double *) {return grInPtr -> Eval(x[0]);}, XIn[0], XIn[NIn - 1], 0);
//	fIn.SetNpx(NIn);
//	fIn.SetLineColor(1);
//	fIn.SetLineWidth(1);
//
//	TGraph grOut(& fIn, "d");

	return grOut;
}


/*  Same as getDerivGraph(), but output TGraph's y-component has polarity reversed.
 */
TGraph getNegDerivGraph(const TGraph * grInPtr) {

	TGraph grDeriv = getDerivGraph(grInPtr);
	TGraph grOut = getAffineTransformedGraph(& grDeriv, -1);

	return grOut;
}


/*  Take an input TGraph and generate another TGraph corresponding to its numerical antiderivative.
 *  Technically, the input TGraph doesn't have to be uniformly sampled. Taking grInPtr -> GetY()
 *  and integrating up from the lower bound element, down from the upper bound element, and then
 *  averaging these results to find an accurate and invertible antiderivative, calling the method
 *  employed here a "superposition method" seems appropriate.
 */
TGraph getAntiderivGraph(const TGraph * grInPtr) {

	int NIn = grInPtr -> GetN();
	double * XIn = grInPtr -> GetX();
	double * YIn = grInPtr -> GetY();

	TGraph grOut(* grInPtr);
	grOut.SetTitle();
	grOut.GetYaxis() -> SetTitle();

	//  Construct array YOut corresponding to antiderivative of YIn.
	//  Integrating up from the lower bound. We know that the lower bound element should correspond to 0.
	double YLow[NIn];
	YLow[0] = 0;
	YLow[1] = YLow[0] + YIn[0] * (XIn[1] - XIn[0]);
	for (int i = 2; i < NIn; ++i) YLow[i] = YLow[i - 2] + YIn[i - 1] * (XIn[i] - XIn[i - 2]);
	YLow[NIn - 1] += YLow[NIn - 2] + YIn[NIn - 1] * (XIn[NIn - 1] - XIn[NIn - 2]);
	YLow[NIn - 1] /= 2;  //  Averaging the two upper bound relations produces the entire integral as approximated by the trapezoid rule.
	//  Integrating down from the upper bound. We know that the upper bound element should correspond to the entire integral.
	double YUp[NIn];
	YUp[NIn - 1] = YLow[NIn - 1];
	YUp[NIn - 2] = YUp[NIn - 1] - YIn[NIn - 1] * (XIn[NIn - 1] - XIn[NIn - 2]);
	for (int i = NIn - 3; i >= 0; --i) YUp[i] = YUp[i + 2] - YIn[i + 1] * (XIn[i + 2] - XIn[i]);
	YUp[0] += YUp[1] - YIn[0] * (XIn[1] - XIn[0]);
	YUp[0] /= 2;  //  Averaging the two lower bound relations produces 0 as determined by the trapezoid rule.
	//  Take the average of the two methods as the antiderivative.
	double * YOut = grOut.GetY();
	for (int i = 0; i < NIn; ++i) YOut[i] = (YLow[i] + YUp[i]) / 2;
//	//  Done with the trapezoidal rule result directly.
//	YOut[0] = 0;
//	for (int i = 1; i < grInPtr -> GetN(); ++i) YOut[i] = YOut[i - 1] + (YIn[i] + YIn[i - 1]) * (XIn[i] - XIn[i - 1]) / 2;

//	TF1 fIn("fIn", [&](double * x, double *) {return grInPtr -> Eval(x[0]);}, XIn[0], XIn[NIn - 1], 0);
//	fIn.SetNpx(NIn);
//	fIn.SetLineColor(1);
//	fIn.SetLineWidth(1);
//
//	TGraph grOut(& fIn, "i");

	return grOut;
}


/*  Modifying scipy.fftpack.diff from SciPy
 *  <https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.fftpack.diff.html>
 *  for usage with FFTtools. Assumes the waveform is uniformly sampled, zero-meaned when
 *  order is negative, isn't constant after differentiation, and the result is expected to be
 *  real over the input domain, so don't use it otherwise.
 */
TGraph getIFFTDiffGraph(const TGraph * grInPtr, double order = 1, int branchOrder = 0) {

	int NIn = grInPtr -> GetN();

	TGraph grOut(* grInPtr);
	grOut.SetTitle();
	grOut.GetYaxis() -> SetTitle();

	//  Construct array YOut corresponding to derivative of YIn.
	int NFFT = NIn / 2 + 1;  //  For a real array of odd length NIn, the largest positive component is indexed at (NIn + 1) / 2 - 1 and no Nyquist component. However, FFTW3 gets around this with unused padding.
	double dw = 2 * M_PI * (1 - 1 / NIn) / (grInPtr -> GetX()[NIn - 1] - grInPtr -> GetX()[0]);
	double cosOrder = cos(order * (0.5 + 2 * branchOrder) * M_PI);  //  i^x = cos(x * (0.5 + 2 * k) * pi) + i * sin(x * (0.5 + 2 * k) * pi), k = ..., -1, 0, 1, ...
	double sinOrder = sin(order * (0.5 + 2 * branchOrder) * M_PI);
	FFTWComplex * YFFT = FFTtools::doFFT(NIn, grInPtr -> GetY());
	FFTWComplex * YDiffFFT = new FFTWComplex[NFFT];
	YDiffFFT[0] = order ? 0 : YFFT[0];  //  Just in case the FFT itself (order == 0) is desired. But should be zero, anyway.
	for (int i = 1; i < NFFT; ++i) {

		double wOrder = pow(i * dw, order);
		YDiffFFT[i].re = wOrder * (cosOrder * YFFT[i].re - sinOrder * YFFT[i].im);
		YDiffFFT[i].im = wOrder * (sinOrder * YFFT[i].re + cosOrder * YFFT[i].im);
	}
	if (fmod(order, 2) && !(NIn % 2)) YDiffFFT[NFFT - 1] = 0;  //  "For odd order and even [NIn], Nyquist mode is taken zero."
	double * YDiff = FFTtools::doInvFFT(NIn, YDiffFFT);
	double * YOut = grOut.GetY();
	std::copy(YDiff, YDiff + NIn, YOut);

	delete [] YDiff, delete [] YDiffFFT, delete [] YFFT;

	return grOut;
}


/*  Same as getNegDerivGraph(), except using getIFFTDiffGraph() instead.
 */
TGraph getIFFTNegDerivGraph(const TGraph * grInPtr) {

	TGraph grDeriv = getIFFTDiffGraph(grInPtr);
	TGraph grOut = getAffineTransformedGraph(& grDeriv, -1);

	return grOut;
}
