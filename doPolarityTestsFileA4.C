/*  Modification of "doPolarityTestsFile.C" for use with ANITA-4 WAIS pulser data.
 *  Started 11/21/17 by John Russell.
 */


#include "FFTtools.h"
#include "waisPolCutTools.h"


/*  Structure corresponding to polarity data from event signals.
 *  Max and min values. Also their sum, difference, and ratio of these (p).
 */
struct polarity {

	double maxTime, maxVal, minTime, minVal, maxMinMeanVal, maxMinHalfWidthVal, p, RMS, maxMinMeanSNR, maxMinHalfWidthSNR, sgndExtrTime, sgndExtrVal, sgndExtrSNR;
};


/*  Shortcut function to fill in instances of polarity structure.
 */
void fillPolarityStruct(polarity & structIn, TGraph grIn) {

	int NIn = grIn.GetN();
	double * XIn = grIn.GetX();
	double * YIn = grIn.GetY();

	int maxIndIn = TMath::LocMax(NIn, YIn);
	structIn.maxTime = XIn[maxIndIn];
	structIn.maxVal = YIn[maxIndIn];  //  Maximum correlation.
	int minIndIn = TMath::LocMin(NIn, YIn);
	structIn.minTime = XIn[minIndIn];
	structIn.minVal = YIn[minIndIn];  //  Minimum correlation.
	structIn.maxMinMeanVal = (structIn.maxVal + structIn.minVal) / 2;
	structIn.maxMinHalfWidthVal = (structIn.maxVal - structIn.minVal) / 2;
	structIn.p = structIn.maxMinMeanVal / structIn.maxMinHalfWidthVal;
	structIn.RMS = grIn.GetRMS(2);
	structIn.maxMinMeanSNR = structIn.maxMinMeanVal / structIn.RMS;
	structIn.maxMinHalfWidthSNR = structIn.maxMinHalfWidthVal / structIn.RMS;
	int sgndExtrIndIn = structIn.maxVal * structIn.maxVal > structIn.minVal * structIn.minVal ? maxIndIn : minIndIn;
	structIn.sgndExtrTime = XIn[sgndExtrIndIn];
	structIn.sgndExtrVal = YIn[sgndExtrIndIn];
	structIn.sgndExtrSNR = structIn.sgndExtrVal / structIn.RMS;
}


void doPolarityTestsFileA4(int runNum = 123, double deltaT = 0.1, double winFullWidth = 5) {
	
	//  To point towards the correct WAIS event files.
	TString ANITA4_ROOT_DATA = gSystem -> Getenv("ANITA4_ROOT_DATA");
	gSystem -> Setenv("ANITA_ROOT_DATA", ANITA4_ROOT_DATA);
	
	// Butterworth low-pass filter in case interpolation introduces unphysical high frequency terms.
	FFTtools::ButterworthFilter but(FFTtools::LOWPASS, 4, 2 * deltaT * 0.65);  //  650 MHz 3dB point for first Butterworth low-pass filter in terms of Nyquist frequency multiplier.

	//  Get pointer addresses to WAIS template TGraph objects.
	TFile templateFile("interpWaisTemplatesA4.root");
	TGraph * templateGraph, * deconvTemplateGraph;
	if (deltaT == 0.1) {

		if (runNum < 129 || runNum > 132) {  //  Accounting for when the pulser is rotated.

			templateGraph = (TGraph *) templateFile.Get("waisHPolTemplateA4_1");
			deconvTemplateGraph = (TGraph *) templateFile.Get("deconvolvedWaisHPolTemplateA4_1");
		} else {

			templateGraph = (TGraph *) templateFile.Get("waisDPolTemplateA4_1");
			deconvTemplateGraph = (TGraph *) templateFile.Get("deconvolvedWaisDPolTemplateA4_1");
		}
	} else {

		if (runNum < 129 || runNum > 132) {

			templateGraph = (TGraph *) templateFile.Get("waisHPolTemplateA4_0");
			deconvTemplateGraph = (TGraph *) templateFile.Get("deconvolvedWaisHPolTemplateA4_0");
		} else {

			templateGraph = (TGraph *) templateFile.Get("waisDPolTemplateA4_0");
			deconvTemplateGraph = (TGraph *) templateFile.Get("deconvolvedWaisDPolTemplateA4_0");
		}
		but.filterGraph(templateGraph);
		but.filterGraph(deconvTemplateGraph);
	}

	//  Create TGraph objects with take negative derivatives of templates.
	TGraph negDerivTemplateGraph = getNegDerivGraph(templateGraph);
	TGraph negDerivDeconvTemplateGraph = getNegDerivGraph(deconvTemplateGraph);

	//  Create TGraph objects which take antiderivative of templates.
	TGraph antiderivTemplateGraph = getAntiderivGraph(templateGraph);
	TGraph antiderivDeconvTemplateGraph = getAntiderivGraph(deconvTemplateGraph);
	
	//  Create TGraph objects which window the template TGraph objects to 5-ns-width windows.
	TGraph windowedTemplateGraph = getWindowedGraph(templateGraph, winFullWidth);
	TGraph windowedDeconvTemplateGraph = getWindowedGraph(deconvTemplateGraph, winFullWidth);

	//  Create TGraph objects with take negative derivative of windowed templates.
	TGraph negDerivWindowedTemplateGraph = getNegDerivGraph(& windowedTemplateGraph);
	TGraph negDerivWindowedDeconvTemplateGraph = getNegDerivGraph(& windowedDeconvTemplateGraph);

	//  Create TGraph objects with take antiderivative of windowed templates.
	TGraph antiderivWindowedTemplateGraph = getAntiderivGraph(& windowedTemplateGraph);
	TGraph antiderivWindowedDeconvTemplateGraph = getAntiderivGraph(& windowedDeconvTemplateGraph);

	//  Open the WAIS headFile so that we can iterate over selected runs and their event numbers.
	TFile waisFile("waisPolCutFilesA4/allWais_max_30005_sinsub_10_3_ad_2.root");
	TTree * waisTree = (TTree *) waisFile.Get("wais");
	AnitaEventSummary * summaryWais = 0;
	waisTree -> SetBranchAddress("summary", & summaryWais);
	waisTree -> SetBranchStatus("*", 0);
	waisTree -> SetBranchStatus("run", 1);
	waisTree -> SetBranchStatus("eventNumber", 1);
	int & run = summaryWais -> run;
	unsigned int & eventNumber = summaryWais -> eventNumber;
	
	//  Create TFile of polarization cuts.
	const char * templateDir = deltaT == 0.1 ? "template_1" : "template_custom";
	const char * startFileName = deltaT == 0.1 ? "template1" : "templateCustom";
	TFile polarityTestsFile(Form("waisPolCutFilesA4/%s/%sPolarityTestsFileRun%d.root", templateDir, startFileName, runNum), "recreate");

	//  Tree which stores event number information common to other trees.
	TTree eventTree("eventTree", Form("TTree storing event numbers for run %d.", runNum));
	
	eventTree.Branch("run", & run);
	eventTree.Branch("eventNumber", & eventNumber);

	//  Phi and theta values for WAIS (waisAngle) and largest interferometric peak (peakAngle).
	TTree waisPointingTree("waisPointingTree",  Form("TTree storing WAIS pointing results from interferometric map of run %d.", runNum));

	struct {

			double waisAngle, peakAngle, peakDisp, peakDev;
	} waisPhiDev, waisThetaDev;

	TString angleDevBranchStruct = "waisAngle/D:peakAngle:peakDisp:peakDev";

	waisPointingTree.Branch("waisPhiDev", & waisPhiDev, angleDevBranchStruct);
	waisPointingTree.Branch("waisThetaDev", & waisThetaDev, angleDevBranchStruct);

	//  These trees and their leaves store the file's predominant data of polarity tests.
	TTree polarityTestsTree("polarityTestsTree", Form("TTree storing polarity test results of waveforms, run %d.", runNum));
	TTree deconvPolarityTestsTree("deconvPolarityTestsTree", Form("TTree storing polarity tests results of deconvolved waveforms, run %d.", runNum));
	TTree windowedPolarityTestsTree("windowedPolarityTestsTree", Form("TTree storing polarity test results of windowed waveforms, run %d.", runNum));
	TTree windowedDeconvPolarityTestsTree("windowedDeconvPolarityTestsTree", Form("TTree storing polarity tests results of windowed, deconvolved waveforms, run %d.", runNum));
	
	polarity waveformStruct, negDerivWaveformStruct, antiderivWaveformStruct;
	polarity corrStruct, negDerivCorrStruct, antiderivCorrStruct;

	polarity deconvWaveformStruct, negDerivDeconvWaveformStruct, antiderivDeconvWaveformStruct;
	polarity deconvCorrStruct, negDerivDeconvCorrStruct, antiderivDeconvCorrStruct;

	polarity windowedWaveformStruct, negDerivWindowedWaveformStruct, antiderivWindowedWaveformStruct;
	polarity windowedCorrStruct, negDerivWindowedCorrStruct, antiderivWindowedCorrStruct;

	polarity windowedDeconvWaveformStruct, negDerivWindowedDeconvWaveformStruct, antiderivWindowedDeconvWaveformStruct;
	polarity windowedDeconvCorrStruct, negDerivWindowedDeconvCorrStruct, antiderivWindowedDeconvCorrStruct;

	TString waveformBranchStruct = "maxTime/D:maxVal:minTime:minVal:maxMinMeanVal:maxMinHalfWidthVal:p:RMS:maxMinMeanSNR:maxMinHalfWidthSNR:sgndExtrTime:sgndExtrVal:sgndExtrSNR";
	TString corrBranchStruct = "maxTime/D:maxCorr:minTime:minCorr:maxMinMeanCorr:maxMinHalfWidthCorr:p:RMS:maxMinMeanSNR:maxMinHalfWidthSNR:sgndExtrTime:sgndExtrCorr:sgndExtrSNR";

	polarityTestsTree.Branch("waveform", & waveformStruct, waveformBranchStruct);
	polarityTestsTree.Branch("negDerivWaveform", & negDerivWaveformStruct, waveformBranchStruct);
	polarityTestsTree.Branch("antiderivWaveform", & antiderivWaveformStruct, waveformBranchStruct);
	polarityTestsTree.Branch("corr", & corrStruct, corrBranchStruct);
	polarityTestsTree.Branch("negDerivCorr", & negDerivCorrStruct, corrBranchStruct);
	polarityTestsTree.Branch("antiderivCorr", & antiderivCorrStruct, corrBranchStruct);

	deconvPolarityTestsTree.Branch("waveform", & deconvWaveformStruct, waveformBranchStruct);
	deconvPolarityTestsTree.Branch("negDerivWaveform", & negDerivDeconvWaveformStruct, waveformBranchStruct);
	deconvPolarityTestsTree.Branch("antiderivWaveform", & antiderivDeconvWaveformStruct, waveformBranchStruct);
	deconvPolarityTestsTree.Branch("corr", & deconvCorrStruct, corrBranchStruct);
	deconvPolarityTestsTree.Branch("negDerivCorr", & negDerivDeconvCorrStruct, corrBranchStruct);
	deconvPolarityTestsTree.Branch("antiderivCorr", & antiderivDeconvCorrStruct, corrBranchStruct);

	windowedPolarityTestsTree.Branch("waveform", & windowedWaveformStruct, waveformBranchStruct);
	windowedPolarityTestsTree.Branch("negDerivWaveform", & negDerivWindowedWaveformStruct, waveformBranchStruct);
	windowedPolarityTestsTree.Branch("antiderivWaveform", & antiderivWindowedWaveformStruct, waveformBranchStruct);
	windowedPolarityTestsTree.Branch("corr", & windowedCorrStruct, corrBranchStruct);
	windowedPolarityTestsTree.Branch("negDerivCorr", & negDerivWindowedCorrStruct, corrBranchStruct);
	windowedPolarityTestsTree.Branch("antiderivCorr", & antiderivWindowedCorrStruct, corrBranchStruct);

	windowedDeconvPolarityTestsTree.Branch("waveform", & windowedDeconvWaveformStruct, waveformBranchStruct);
	windowedDeconvPolarityTestsTree.Branch("negDerivWaveform", & negDerivWindowedDeconvWaveformStruct, waveformBranchStruct);
	windowedDeconvPolarityTestsTree.Branch("antiderivWaveform", & antiderivWindowedDeconvWaveformStruct, waveformBranchStruct);
	windowedDeconvPolarityTestsTree.Branch("corr", & windowedDeconvCorrStruct, corrBranchStruct);
	windowedDeconvPolarityTestsTree.Branch("negDerivCorr", & negDerivWindowedDeconvCorrStruct, corrBranchStruct);
	windowedDeconvPolarityTestsTree.Branch("antiderivCorr", & antiderivWindowedDeconvCorrStruct, corrBranchStruct);

	//  Set up AnitaDataset object for given runNum.
	AnitaDataset d(runNum);
	
	//  UCorrelator stuff. Set up Analyzer objects.
	UCorrelator::AnalysisConfig cfg;
	cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution;  //  Rather than use the band-limited default, use all-pass to account for lack of elephant trunk corrections.
	cfg.response_option = UCorrelator::AnalysisConfig::ResponseTUFF;  //  Deconvolving accounting for TUFF filter response.
//	cfg.response_option = UCorrelator::AnalysisConfig::ResponseSingleBRotter;  //  Deconvolving with the Ben Rotter's A3 impulse response.
	UCorrelator::Analyzer analyzer(& cfg, true);  //  "true" is set so "analyzer" doesn't reset each time.
	FilterStrategy * fStrat = UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2");  //  Basic sine subtraction filter strategy.

	//  Fill polarityTestsTree with data.
	for (int entryNum = 0; entryNum < waisTree -> GetEntries(); ++entryNum) {
		
		waisTree -> GetEntry(entryNum);
		
		if (run > runNum) break;
		if (run == runNum) {
			
			//  Use the analyzer declared above to acquire coherently summed waveform for given run and event number.
			d.getEvent(eventNumber);  //  In the run, pointing to the correct eventNumber.
			FilteredAnitaEvent filtEvent(d.useful(), fStrat, d.gps(), d.header());
			AnitaEventSummary eventSum;
			analyzer.analyze(& filtEvent, & eventSum);

			//  Fill branches corresponding to WAIS pointing values.
			waisPhiDev.waisAngle = eventSum.wais.phi;
			waisPhiDev.peakAngle = eventSum.peak[AnitaPol::kHorizontal][0].phi;
			waisPhiDev.peakDisp = FFTtools::wrap(waisPhiDev.peakAngle - waisPhiDev.waisAngle, 360, 0);  // fmod() doesn't behave as expected...
//			waisPhiDev.peakDisp = eventSum.peak[AnitaPol::kHorizontal][0].dPhiWais();
//			waisPhiDev.peakDisp = fmod(waisPhiDev.peakAngle - waisPhiDev.waisAngle, 360);
			waisPhiDev.peakDev = std::abs(waisPhiDev.peakDisp);

			waisThetaDev.waisAngle = eventSum.wais.theta;
			waisThetaDev.peakAngle = eventSum.peak[AnitaPol::kHorizontal][0].theta;
			waisThetaDev.peakDisp = waisThetaDev.peakAngle - waisThetaDev.waisAngle;
			waisThetaDev.peakDev = std::abs(waisThetaDev.peakDisp);

			//  The AnalysisWaveform class which the Analyzer class references makes constant objects.
			//  Annoyingly, we need to work around this.
			//  FFTtools creates NEW TGraph objects.
			//  Fill branches corresponding to entire event waveforms.
			TGraph * waveformWais = FFTtools::getInterpolatedGraph(analyzer.getCoherent(AnitaPol::kHorizontal, 0) -> even(), deltaT);
			but.filterGraph(waveformWais);
			fillPolarityStruct(waveformStruct, * waveformWais);

			TGraph * deconvWaveformWais = FFTtools::getInterpolatedGraph(analyzer.getDeconvolved(AnitaPol::kHorizontal, 0) -> even(), deltaT);
			but.filterGraph(deconvWaveformWais);
			fillPolarityStruct(deconvWaveformStruct, * deconvWaveformWais);

			//  Fill branches corresponding to negative derivative of entire event waveforms.
			TGraph negDerivWaveformWais = getNegDerivGraph(waveformWais);
			fillPolarityStruct(negDerivWaveformStruct, negDerivWaveformWais);

			TGraph negDerivDeconvWaveformWais = getNegDerivGraph(deconvWaveformWais);
			fillPolarityStruct(negDerivDeconvWaveformStruct, negDerivDeconvWaveformWais);

			//  Fill branches corresponding to antiderivative of entire event waveforms.
			TGraph antiderivWaveformWais = getAntiderivGraph(waveformWais);
			fillPolarityStruct(antiderivWaveformStruct, antiderivWaveformWais);

			TGraph antiderivDeconvWaveformWais = getAntiderivGraph(deconvWaveformWais);
			fillPolarityStruct(antiderivDeconvWaveformStruct, antiderivDeconvWaveformWais);

			//  Fill branches corresponding to correlation with entire event waveforms.
			TGraph corrWais = getCorrGraph(templateGraph, waveformWais);
			fillPolarityStruct(corrStruct, corrWais);

			TGraph deconvCorrWais = getCorrGraph(deconvTemplateGraph, deconvWaveformWais);
			fillPolarityStruct(deconvCorrStruct, deconvCorrWais);

			//  Fill branches corresponding to correlation with negative derivative of entire event waveforms.
			TGraph negDerivCorrWais = getCorrGraph(& negDerivTemplateGraph, & negDerivWaveformWais);
			fillPolarityStruct(negDerivCorrStruct, negDerivCorrWais);

			TGraph negDerivDeconvCorrWais = getCorrGraph(& negDerivDeconvTemplateGraph, & negDerivDeconvWaveformWais);
			fillPolarityStruct(negDerivDeconvCorrStruct, negDerivDeconvCorrWais);

			//  Fill branches corresponding to correlation with antiderivative of entire event waveforms.
			TGraph antiderivCorrWais = getCorrGraph(& antiderivTemplateGraph, & antiderivWaveformWais);
			fillPolarityStruct(antiderivCorrStruct, antiderivCorrWais);

			TGraph antiderivDeconvCorrWais = getCorrGraph(& antiderivDeconvTemplateGraph, & antiderivDeconvWaveformWais);
			fillPolarityStruct(antiderivDeconvCorrStruct, antiderivDeconvCorrWais);

			//  Create windowed graphs of event waveforms.
			//  Fill branches corresponding to windowed event waveform.
			TGraph windowedWaveformWais = getWindowedGraph(waveformWais, winFullWidth);
			fillPolarityStruct(windowedWaveformStruct, windowedWaveformWais);

			TGraph windowedDeconvWaveformWais = getWindowedGraph(deconvWaveformWais, winFullWidth);
			fillPolarityStruct(windowedDeconvWaveformStruct, windowedDeconvWaveformWais);

			//  Fill branches corresponding to negative derivative of windowed event waveforms.
			TGraph negDerivWindowedWaveformWais = getNegDerivGraph(& windowedWaveformWais);
			fillPolarityStruct(negDerivWindowedWaveformStruct, negDerivWindowedWaveformWais);

			TGraph negDerivWindowedDeconvWaveformWais = getNegDerivGraph(& windowedDeconvWaveformWais);
			fillPolarityStruct(negDerivWindowedDeconvWaveformStruct, negDerivWindowedDeconvWaveformWais);

			//  Fill branches corresponding to antiderivative of windowed event waveforms.
			TGraph antiderivWindowedWaveformWais = getAntiderivGraph(& windowedWaveformWais);
			fillPolarityStruct(antiderivWindowedWaveformStruct, antiderivWindowedWaveformWais);

			TGraph antiderivWindowedDeconvWaveformWais = getAntiderivGraph(& windowedDeconvWaveformWais);
			fillPolarityStruct(antiderivWindowedDeconvWaveformStruct, antiderivWindowedDeconvWaveformWais);

			//  Fill branches corresponding to correlation with windowed event waveforms.
			TGraph windowedCorrWais = getCorrGraph(& windowedTemplateGraph, & windowedWaveformWais);
			fillPolarityStruct(windowedCorrStruct, windowedCorrWais);
			
			TGraph deconvWindowedCorrWais = getCorrGraph(& windowedDeconvTemplateGraph, & windowedDeconvWaveformWais);
			fillPolarityStruct(windowedDeconvCorrStruct, deconvWindowedCorrWais);

			//  Fill branches corresponding to correlations with negative derivative of windowed event waveforms.
			TGraph negDerivWindowedCorrWais = getCorrGraph(& negDerivWindowedTemplateGraph, & negDerivWindowedWaveformWais);
			fillPolarityStruct(negDerivWindowedCorrStruct, negDerivWindowedCorrWais);

			TGraph negDerivWindowedDeconvCorrWais = getCorrGraph(& negDerivWindowedDeconvTemplateGraph, & negDerivWindowedDeconvWaveformWais);
			fillPolarityStruct(negDerivWindowedDeconvCorrStruct, negDerivWindowedDeconvCorrWais);

			//  Fill branches corresponding to correlations with antiderivative of windowed event waveforms.
			TGraph antiderivWindowedCorrWais = getCorrGraph(& antiderivWindowedTemplateGraph, & antiderivWindowedWaveformWais);
			fillPolarityStruct(antiderivWindowedCorrStruct, antiderivWindowedCorrWais);

			TGraph antiderivWindowedDeconvCorrWais = getCorrGraph(& antiderivWindowedDeconvTemplateGraph, & antiderivWindowedDeconvWaveformWais);
			fillPolarityStruct(antiderivWindowedDeconvCorrStruct, antiderivWindowedDeconvCorrWais);

			//  Filling the trees.
			eventTree.Fill();
			waisPointingTree.Fill();
			polarityTestsTree.Fill();
			deconvPolarityTestsTree.Fill();
			windowedPolarityTestsTree.Fill();
			windowedDeconvPolarityTestsTree.Fill();

			//  Further annoyance with the necessity of memory clearance.
			analyzer.clearInteractiveMemory();
			delete waveformWais, delete deconvWaveformWais;
		}
	}
	
	//  Writing trees to polarityTestsFile.
	polarityTestsFile.cd();
	eventTree.Write();
	waisPointingTree.Write();
	polarityTestsTree.Write();
	deconvPolarityTestsTree.Write();
	windowedPolarityTestsTree.Write();
	windowedDeconvPolarityTestsTree.Write();
	
	//  Closing files.
	polarityTestsFile.Close();
	waisFile.Close();
	templateFile.Close();

	//  Freeing up memory.
	if (deltaT != 0.1) delete templateGraph, delete deconvTemplateGraph;
}
