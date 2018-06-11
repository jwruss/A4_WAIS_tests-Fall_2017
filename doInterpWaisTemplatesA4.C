/*  Based off of "doInterpWaisHPolTemplatesA3.C", but modifies accordingly for ANITA-4.
 *  Started 11/8/17.
 */


#include "FFTtools.h"


void doInterpWaisTemplatesA4(double deltaT = 0.1) {

	//  To point towards the correct WAIS event files.
	TString ANITA4_ROOT_DATA = gSystem -> Getenv("ANITA4_ROOT_DATA");
	gSystem -> Setenv("ANITA_ROOT_DATA", ANITA4_ROOT_DATA);

	//  Open the WAIS headFile so that we can iterate over selected runs and their event numbers.
	TFile waisFile("waisPolCutFiles/allWais_max_30005_sinsub_10_3_ad_2.root");
	TTree * waisTree = (TTree *) waisFile.Get("wais");
	RawAnitaHeader * headerWais = 0;
	waisTree -> SetBranchAddress("header", & headerWais);
	waisTree -> SetBranchStatus("*", 0);
	waisTree -> SetBranchStatus("run", 1);
	waisTree -> SetBranchStatus("eventNumber", 1);
	int & run = headerWais -> run;
	unsigned int & eventNumber = headerWais -> eventNumber;
	int numEntries = waisTree -> GetEntries();
	int numEntriesHV = waisTree -> GetEntries("run < 129 || run > 132");
	int numEntriesAD = waisTree -> GetEntries("run >= 129 && run <= 132");
	int unitStep = 1000;
	int numUnitEntriesHV = numEntriesHV / unitStep;
	int numUnitEntriesAD = numEntriesAD / unitStep;

	//  Create TFile of deconvolved WAIS template.
	TFile templateFile("interpWaisTemplatesA4.root", "recreate");

	//  UCorrelator stuff. Set up Analyzer objects.
	UCorrelator::AnalysisConfig cfg;  //  Default configuration (cfg.deconvolution_method) is with all-pass deconvolution.
	cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution;  //  Rather than use the band-limited default, use all-pass to account for lack of elephant trunk corrections.
	cfg.response_option = UCorrelator::AnalysisConfig::ResponseTUFF;  //  Deconvolving accounting for TUFF filter response.
//	cfg.response_option = UCorrelator::AnalysisConfig::ResponseSingleBRotter;  //  Deconvolving with the Ben Rotter's A3 impulse response.	UCorrelator::Analyzer analyzer(& cfg, true);  //  "true" is set so "analyzer" doesn't reset each time.
	UCorrelator::Analyzer analyzer(& cfg, true);  //  "true" is set so "analyzer" doesn't reset each time.
	FilterStrategy * fStrat = UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2");  //  Sine subtraction filter strategy used throughout UCorrelator.

	//  Create array, then fill with TGraphs.
	TGraph * cohWaisH[numUnitEntriesHV], * deconvWaisH[numUnitEntriesHV], * cohWaisV[numUnitEntriesHV], * deconvWaisV[numUnitEntriesHV];
	TGraph * cohWaisD[numUnitEntriesAD], * deconvWaisD[numUnitEntriesAD], * cohWaisA[numUnitEntriesAD], * deconvWaisA[numUnitEntriesAD];
	int idxHV = 0, idxAD = 0;  //  Indices for filling above arrays.
	for (int entryNum = 10; entryNum < numEntries; entryNum += unitStep) {

		waisTree -> GetEntry(entryNum);

		int runNum = run;

		AnitaDataset d(runNum);

		do {

			//  Use the analyzer declared above to acquire coherently summed waveform for given run and event number.
			d.getEvent(eventNumber);  //  In the run, pointing to the correct eventNumber.
			FilteredAnitaEvent filtEvent(d.useful(), fStrat, d.gps(), d.header());
			AnitaEventSummary eventSum;
			analyzer.analyze(& filtEvent, & eventSum);
//			AnitaEventSummary * eventSum = new AnitaEventSummary;
//			analyzer.analyze(& filtEvent, eventSum);

			//  Store the coherently summed waveform and its convolution into array elements.
			if (runNum < 129 || runNum > 132) {

				cohWaisH[idxHV] = new TGraph(* analyzer.getCoherent(AnitaPol::kHorizontal, 0) -> even());
				deconvWaisH[idxHV] = new TGraph(* analyzer.getDeconvolved(AnitaPol::kHorizontal, 0) -> even());
				cohWaisV[idxHV] = new TGraph(* analyzer.getCoherent(AnitaPol::kVertical, 0) -> even());
				deconvWaisV[idxHV] = new TGraph(* analyzer.getDeconvolved(AnitaPol::kVertical, 0) -> even());
				++idxHV;
			} else {

				cohWaisD[idxAD] = new TGraph(* analyzer.getCoherent(AnitaPol::kHorizontal, 0) -> even());
				deconvWaisD[idxAD] = new TGraph(* analyzer.getDeconvolved(AnitaPol::kHorizontal, 0) -> even());
				cohWaisA[idxAD] = new TGraph(* analyzer.getCoherent(AnitaPol::kVertical, 0) -> even());
				deconvWaisA[idxAD] = new TGraph(* analyzer.getDeconvolved(AnitaPol::kVertical, 0) -> even());
				++idxAD;
			}

			//  Further annoyance with the necessity of memory clearance.
			analyzer.clearInteractiveMemory();
//			delete eventSum;

			entryNum += unitStep;  //  Move on to next event in run.

			//  Handling end of waisTree.
			if (entryNum < numEntries) waisTree -> GetEntry(entryNum);
			else ++run;

		} while (runNum == run);

		entryNum -= unitStep;  //  So that no events are accidentally skipped over from last entryNum iteration inside the do-while loop.
	}

	//  Setting up Butterworth filters, and their descriptions to go into titles.
	double f3dBMult1 = 2 * deltaT * 0.65;  //  650 MHz 3dB point for first Butterworth low-pass filter in terms of Nyquist frequency multiplier.
	FFTtools::ButterworthFilter but1(FFTtools::LOWPASS, 4, f3dBMult1);
	TString but1Str = Form("{Butterworth low-pass, order 4, f_{3dB} = %g*f_{Nyq} (f_{3dB} = 650 MHz, #Deltat = %g ns)}", f3dBMult1, deltaT);

//	double f3dBMult2 = 2 * deltaT * 0.45;  //  450 MHz 3dB point for second Butterworth low-pass filter.
//	FFTtools::ButterworthFilter but2(FFTtools::LOWPASS, 4, f3dBMult2);
//	TString but2Str = Form("{Butterworth low-pass, order 4, f_{3dB} = %g*f_{Nyq} (f_{3dB} = 450 MHz, #Deltat = %g ns)}", f3dBMult2, deltaT);
//
//	double fCMult = 2 * deltaT * 0.7;  //  700 MHz central frequency of Butterworth band-pass filters in terms of Nyquist frequency multiplier.
//	double fHWHMMult = 2 * deltaT * 0.5;  //  500 MHz half width at half maximum frequency for Butterworth band-pass filters in terms of Nyquist frequency multiplier.
//	FFTtools::ButterworthFilter but3(FFTtools::BANDPASS, 2, fCMult, fHWHMMult);
//	TString but3Str = Form("{Butterworth band-pass, order 2, f_{C} = %g*f_{Nyq}, f_{HWHM} = %g*f_{Nyq} (f_{C} = 700 MHz, f_{HWHM} = 500 MHz, #Deltat = %g ns)}", fCMult, fHWHMMult, deltaT);
//	FFTtools::ButterworthFilter but4(FFTtools::BANDPASS, 4, fCMult, fHWHMMult);
//	TString but4Str = Form("{Butterworth band-pass, order 4, f_{C} = %g*f_{Nyq}, f_{HWHM} = %g*f_{Nyq} (f_{C} = 700 MHz, f_{HWHM} = 500 MHz, #Deltat = %g ns)}", fCMult, fHWHMMult, deltaT);

	//  Butterworth filter group delays.
	TString groupDelayStr = "#splitline{Corresponding group delay}";

	TGraph * groupDelay1 = but1.groupDelay();
	groupDelay1 -> SetTitle(groupDelayStr + but1Str);

//	TGraph * groupDelay2 = but2.groupDelay();
//	groupDelay2 -> SetTitle(groupDelayStr + but2Str);
//
//	TGraph * groupDelay3 = but3.groupDelay();
//	groupDelay3 -> SetTitle(groupDelayStr + but3Str);
//
//	TGraph * groupDelay4 = but4.groupDelay();
//	groupDelay4 -> SetTitle(groupDelayStr + but4Str);

	//  Coherently summed waveforms.
	TString cohHGraphStr = "#splitline{WAIS HPol template A4}";

	TGraph * cohHGraph0 = FFTtools::interpolateCorrelateAndAverage(deltaT, numUnitEntriesHV, cohWaisH);
	cohHGraph0 -> SetTitle(cohHGraphStr + "{Unfiltered}");
	cohHGraph0 -> GetXaxis() -> SetTitle("Time (ns)");
	cohHGraph0 -> GetYaxis() -> SetTitle("Voltage (mV)");

	TGraph * cohHGraph1 = (TGraph *) cohHGraph0 -> Clone();  //  We use "Clone()" instead of "new TGraph(* )"  because the latter doesn't seem to copy "SetObject()" commands correctly.
	but1.filterGraph(cohHGraph1);
	cohHGraph1 -> SetTitle(cohHGraphStr + but1Str);

	TString cohVGraphStr = "#splitline{WAIS VPol template A4}";

	TGraph * cohVGraph0 = FFTtools::interpolateCorrelateAndAverage(deltaT, numUnitEntriesHV, cohWaisV);
	cohVGraph0 -> SetTitle(cohVGraphStr + "{Unfiltered}");
	cohVGraph0 -> GetXaxis() -> SetTitle("Time (ns)");
	cohVGraph0 -> GetYaxis() -> SetTitle("Voltage (mV)");

	TGraph * cohVGraph1 = (TGraph *) cohVGraph0 -> Clone();  //  We use "Clone()" instead of "new TGraph(* )"  because the latter doesn't seem to copy "SetObject()" commands correctly.
	but1.filterGraph(cohVGraph1);
	cohVGraph1 -> SetTitle(cohVGraphStr + but1Str);

	TString cohDGraphStr = "#splitline{WAIS DPol template A4}";

	TGraph * cohDGraph0 = FFTtools::interpolateCorrelateAndAverage(deltaT, numUnitEntriesAD, cohWaisD);
	cohDGraph0 -> SetTitle(cohDGraphStr + "{Unfiltered}");
	cohDGraph0 -> GetXaxis() -> SetTitle("Time (ns)");
	cohDGraph0 -> GetYaxis() -> SetTitle("Voltage (mV)");

	TGraph * cohDGraph1 = (TGraph *) cohDGraph0 -> Clone();  //  We use "Clone()" instead of "new TGraph(* )"  because the latter doesn't seem to copy "SetObject()" commands correctly.
	but1.filterGraph(cohDGraph1);
	cohDGraph1 -> SetTitle(cohDGraphStr + but1Str);

	TString cohAGraphStr = "#splitline{WAIS APol template A4}";

	TGraph * cohAGraph0 = FFTtools::interpolateCorrelateAndAverage(deltaT, numUnitEntriesAD, cohWaisA);
	cohAGraph0 -> SetTitle(cohAGraphStr + "{Unfiltered}");
	cohAGraph0 -> GetXaxis() -> SetTitle("Time (ns)");
	cohAGraph0 -> GetYaxis() -> SetTitle("Voltage (mV)");

	TGraph * cohAGraph1 = (TGraph *) cohAGraph0 -> Clone();  //  We use "Clone()" instead of "new TGraph(* )"  because the latter doesn't seem to copy "SetObject()" commands correctly.
	but1.filterGraph(cohAGraph1);
	cohAGraph1 -> SetTitle(cohAGraphStr + but1Str);
//	TString cohGraphStr = "#splitline{WAIS HPol template A4}";
//
//	TGraph * cohGraph0 = FFTtools::interpolateCorrelateAndAverage(deltaT, numUnitEntries, cohWaisH);
//	cohGraph0 -> SetTitle(cohGraphStr + "{Unfiltered}");
//	cohGraph0 -> GetXaxis() -> SetTitle("Time (ns)");
//	cohGraph0 -> GetYaxis() -> SetTitle("Voltage (mV)");
//
//	TGraph * cohGraph1 = (TGraph *) cohGraph0 -> Clone();  //  We use "Clone()" instead of "new TGraph(* )"  because the latter doesn't seem to copy "SetObject()" commands correctly.
//	but1.filterGraph(cohGraph1);
//	cohGraph1 -> SetTitle(cohGraphStr + but1Str);
//
//	TGraph * cohGraph2 = (TGraph *) cohGraph0 -> Clone();
//	but2.filterGraph(cohGraph2);
//	cohGraph2 -> SetTitle(cohGraphStr + but2Str);
//
//	TGraph * cohGraph3 = (TGraph *) cohGraph0 -> Clone();
//	but3.filterGraph(cohGraph3);
//	cohGraph3 -> SetTitle(cohGraphStr + but3Str);
//
//	TGraph * cohGraph4 = (TGraph *) cohGraph0 -> Clone();
//	but4.filterGraph(cohGraph4);
//	cohGraph4 -> SetTitle(cohGraphStr + but4Str);

	//  Deconvolved coherently summed waveforms.
	TString deconvHGraphStr = "#splitline{Deconvolved WAIS HPol template A4}";

	TGraph * deconvHGraph0 = FFTtools::interpolateCorrelateAndAverage(deltaT, numUnitEntriesHV, deconvWaisH);
	deconvHGraph0 -> SetTitle(deconvHGraphStr + "{Unfiltered}");
	deconvHGraph0 -> GetXaxis() -> SetTitle("Time (ns)");
	deconvHGraph0 -> GetYaxis() -> SetTitle("Electric field strength (V/m)");

	TGraph * deconvHGraph1 = (TGraph *) deconvHGraph0 -> Clone();  //  We use "Clone()" instead of "new TGraph(* )"  because the latter doesn't seem to copy "SetObject()" commands correctly.
	but1.filterGraph(deconvHGraph1);
	deconvHGraph1 -> SetTitle(deconvHGraphStr + but1Str);

	TString deconvVGraphStr = "#splitline{Deconvolved WAIS VPol template A4}";

	TGraph * deconvVGraph0 = FFTtools::interpolateCorrelateAndAverage(deltaT, numUnitEntriesHV, deconvWaisV);
	deconvVGraph0 -> SetTitle(deconvVGraphStr + "{Unfiltered}");
	deconvVGraph0 -> GetXaxis() -> SetTitle("Time (ns)");
	deconvVGraph0 -> GetYaxis() -> SetTitle("Electric field strength (V/m)");

	TGraph * deconvVGraph1 = (TGraph *) deconvVGraph0 -> Clone();  //  We use "Clone()" instead of "new TGraph(* )"  because the latter doesn't seem to copy "SetObject()" commands correctly.
	but1.filterGraph(deconvVGraph1);
	deconvVGraph1 -> SetTitle(deconvVGraphStr + but1Str);

	TString deconvDGraphStr = "#splitline{Deconvolved WAIS DPol template A4}";

	TGraph * deconvDGraph0 = FFTtools::interpolateCorrelateAndAverage(deltaT, numUnitEntriesAD, deconvWaisD);
	deconvDGraph0 -> SetTitle(deconvDGraphStr + "{Unfiltered}");
	deconvDGraph0 -> GetXaxis() -> SetTitle("Time (ns)");
	deconvDGraph0 -> GetYaxis() -> SetTitle("Electric field strength (V/m)");

	TGraph * deconvDGraph1 = (TGraph *) deconvDGraph0 -> Clone();  //  We use "Clone()" instead of "new TGraph(* )"  because the latter doesn't seem to copy "SetObject()" commands correctly.
	but1.filterGraph(deconvDGraph1);
	deconvDGraph1 -> SetTitle(deconvDGraphStr + but1Str);

	TString deconvAGraphStr = "#splitline{Deconvolved WAIS APol template A4}";

	TGraph * deconvAGraph0 = FFTtools::interpolateCorrelateAndAverage(deltaT, numUnitEntriesAD, deconvWaisA);
	deconvAGraph0 -> SetTitle(cohAGraphStr + "{Unfiltered}");
	deconvAGraph0 -> GetXaxis() -> SetTitle("Time (ns)");
	deconvAGraph0 -> GetYaxis() -> SetTitle("Electric field strength (V/m)");

	TGraph * deconvAGraph1 = (TGraph *) deconvAGraph0 -> Clone();  //  We use "Clone()" instead of "new TGraph(* )"  because the latter doesn't seem to copy "SetObject()" commands correctly.
	but1.filterGraph(deconvAGraph1);
	deconvAGraph1 -> SetTitle(deconvAGraphStr + but1Str);
//	TString deconvGraphStr = "#splitline{Deconvolved WAIS HPol template A4}";
//
//	TGraph * deconvGraph0 = FFTtools::interpolateCorrelateAndAverage(deltaT, numUnitEntries, deconvWaisH);
//	deconvGraph0 -> SetTitle(deconvGraphStr + "{Unfiltered}");
//	deconvGraph0 -> GetXaxis() -> SetTitle("Time (ns)");
//	deconvGraph0 -> GetYaxis() -> SetTitle("Electric field (mV/m)");
//
//	TGraph * deconvGraph1 = (TGraph *) deconvGraph0 -> Clone();
//	but1.filterGraph(deconvGraph1);
//	deconvGraph1 -> SetTitle(deconvGraphStr + but1Str);
//
//	TGraph * deconvGraph2 = (TGraph *) deconvGraph0 -> Clone();
//	but2.filterGraph(deconvGraph2);
//	deconvGraph2 -> SetTitle(deconvGraphStr + but2Str);
//
//	TGraph * deconvGraph3 = (TGraph *) deconvGraph0 -> Clone();
//	but3.filterGraph(deconvGraph3);
//	deconvGraph3 -> SetTitle(deconvGraphStr + but3Str);
//
//	TGraph * deconvGraph4 = (TGraph *) deconvGraph0 -> Clone();
//	but4.filterGraph(deconvGraph4);
//	deconvGraph4 -> SetTitle(deconvGraphStr + but4Str);

	//  Writing TGraphs to templateFile.
	templateFile.cd();

	groupDelay1 -> Write("groupDelay_1");
//	groupDelay2 -> Write("groupDelay_2");
//	groupDelay3 -> Write("groupDelay_3");
//	groupDelay4 -> Write("groupDelay_4");

	cohHGraph0 -> Write("waisHPolTemplateA4_0");
	cohHGraph1 -> Write("waisHPolTemplateA4_1");
	cohVGraph0 -> Write("waisVPolTemplateA4_0");
	cohVGraph1 -> Write("waisVPolTemplateA4_1");
	cohDGraph0 -> Write("waisDPolTemplateA4_0");
	cohDGraph1 -> Write("waisDPolTemplateA4_1");
	cohAGraph0 -> Write("waisAPolTemplateA4_0");
	cohAGraph1 -> Write("waisAPolTemplateA4_1");
//	cohGraph0 -> Write("waisHPolTemplateA4_0");
//	cohGraph1 -> Write("waisHPolTemplateA4_1");
//	cohGraph2 -> Write("waisHPolTemplateA4_2");
//	cohGraph3 -> Write("waisHPolTemplateA4_3");
//	cohGraph4 -> Write("waisHPolTemplateA4_4");

	deconvHGraph0 -> Write("deconvolvedWaisHPolTemplateA4_0");
	deconvHGraph1 -> Write("deconvolvedWaisHPolTemplateA4_1");
	deconvVGraph0 -> Write("deconvolvedWaisVPolTemplateA4_0");
	deconvVGraph1 -> Write("deconvolvedWaisVPolTemplateA4_1");
	deconvDGraph0 -> Write("deconvolvedWaisDPolTemplateA4_0");
	deconvDGraph1 -> Write("deconvolvedWaisDPolTemplateA4_1");
	deconvAGraph0 -> Write("deconvolvedWaisAPolTemplateA4_0");
	deconvAGraph1 -> Write("deconvolvedWaisAPolTemplateA4_1");
//	deconvGraph0 -> Write("deconvolvedWaisHPolTemplateA4_0");
//	deconvGraph1 -> Write("deconvolvedWaisHPolTemplateA4_1");
//	deconvGraph2 -> Write("deconvolvedWaisHPolTemplateA4_2");
//	deconvGraph3 -> Write("deconvolvedWaisHPolTemplateA4_3");
//	deconvGraph4 -> Write("deconvolvedWaisHPolTemplateA4_4");

	//  Closing files.
	templateFile.Close();
	waisFile.Close();

	//  Freeing up memory.
	delete groupDelay1;
	delete cohHGraph0, delete cohHGraph1, delete cohVGraph0, delete cohVGraph1;
	delete cohDGraph0, delete cohDGraph1, delete cohAGraph0, delete cohAGraph1;
	delete deconvHGraph0, delete deconvHGraph1, delete deconvVGraph0, delete deconvVGraph1;
	delete deconvDGraph0, delete deconvDGraph1, delete deconvAGraph0, delete deconvAGraph1;
	for (int unitEntryNumHV = 0; unitEntryNumHV < numUnitEntriesHV; ++unitEntryNumHV) {
		delete cohWaisH[unitEntryNumHV], delete cohWaisV[unitEntryNumHV];
		delete deconvWaisH[unitEntryNumHV], delete deconvWaisV[unitEntryNumHV];
	}
	for (int unitEntryNumAD = 0; unitEntryNumAD < numUnitEntriesAD; ++unitEntryNumAD) {
		delete cohWaisD[unitEntryNumAD], delete cohWaisA[unitEntryNumAD];
		delete deconvWaisD[unitEntryNumAD], delete deconvWaisA[unitEntryNumAD];
	}
//	delete groupDelay1, delete groupDelay2, delete groupDelay3, delete groupDelay4;
//	delete cohGraph0, delete cohGraph1, delete cohGraph2, delete cohGraph3, delete cohGraph4;
//	delete deconvGraph0, delete deconvGraph1, delete deconvGraph2, delete deconvGraph3, delete deconvGraph4;
//	for (int unitEntryNum = 0; unitEntryNum < numUnitEntries; ++unitEntryNum) {
//		delete cohWaisH [unitEntryNum], delete deconvWaisH [unitEntryNum];
//	}
}
