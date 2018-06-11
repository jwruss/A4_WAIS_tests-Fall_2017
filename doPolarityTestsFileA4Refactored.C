/*  Here, we create a refactor of the files produced by "doPolarityTestsA4.C", except now the
 *  branches from "polarityTestsTree", "deconvPolarityTestsTree", "windowedPolarityTestTree",
 *  and "windowedDeconvPolarityTestsTree" into a single tree under corresponding parent
 *  branches. Started 12/12/17 by John Russell.
 */


#include "waisPolCutTools.h"


/*  Class corresponding to polarity data from event signals.
 *  Max and min values. Also their sum, difference, and ratio of these (p).
 */
class polarity: public TObject {

	public:
		double maxTime, maxVal, minTime, minVal, maxMinMeanVal, maxMinHalfWidthVal;
		double p, RMS, maxMinMeanSNR, maxMinHalfWidthSNR;
		double sgndExtrTime, sgndExtrVal, sgndExtrSNR;
		polarity() {};
		ClassDef(polarity, 1);
};

/*  Class used to construct refactored tree branches and their subbranches.
 */
class parentBranch {

	public:
		polarity waveform, negDerivWaveform, antiderivWaveform, corr, negDerivCorr, antiderivCorr;
} waveformTests, deconvTests, windowedTests, windowedDeconvTests;


/*  Shortcut function to fill in instances of polarity class.
 */
void fillPolarityClass(polarity classIn, TBranch * brIn, TGraph grIn) {

	int NIn = grIn.GetN();
	double * XIn = grIn.GetX();
	double * YIn = grIn.GetY();

	int maxIndIn = TMath::LocMax(NIn, YIn);
	classIn.maxTime = XIn[maxIndIn];
	classIn.maxVal = YIn[maxIndIn];  //  Maximum correlation.
	int minIndIn = TMath::LocMin(NIn, YIn);
	classIn.minTime = XIn[minIndIn];
	classIn.minVal = YIn[minIndIn];  //  Minimum correlation.
	classIn.maxMinMeanVal = (classIn.maxVal + classIn.minVal) / 2;
	classIn.maxMinHalfWidthVal = (classIn.maxVal - classIn.minVal) / 2;
	classIn.p = classIn.maxMinMeanVal / classIn.maxMinHalfWidthVal;
	classIn.RMS = grIn.GetRMS(2);
	classIn.maxMinMeanSNR = classIn.maxMinMeanVal / classIn.RMS;
	classIn.maxMinHalfWidthSNR = classIn.maxMinHalfWidthVal / classIn.RMS;
	int sgndExtrIndIn = classIn.maxVal * classIn.maxVal > classIn.minVal * classIn.minVal ? maxIndIn : minIndIn;
	classIn.sgndExtrTime = XIn[sgndExtrIndIn];
	classIn.sgndExtrVal = YIn[sgndExtrIndIn];
	classIn.sgndExtrSNR = classIn.sgndExtrVal / classIn.RMS;

	brIn -> Fill();
}


//void fillPolarityClass(polarity classIn, TBranch * brIn, TGraph grIn) {
//
//	fillPolarityClass(classIn, grIn);
//	brIn -> Fill();
//}


void doPolarityTestsFileA4Refactored(int runNum = 123, double deltaT = 0.1, double winFullWidth = 5) {
	
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
	TFile waisFile("allWais_max_30005_sinsub_10_3_ad_2.root");
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
	TFile polarityTestsFile(Form("waisPolCutFilesA4/%s/%sPolarityTestsFileRun%dRefactored.root", templateDir, startFileName, runNum), "recreate");

	//  Phi and theta values for WAIS (waisAngle) and largest interferometric peak (peakAngle).
	TTree waisPointingTree("waisPointingTree",  Form("TTree storing WAIS pointing results from interferometric map of run %d.", runNum));

	class {

		public:

			double waisAngle, peakAngle, peakDisp, peakDev;
	} waisPhiDev, waisThetaDev;

	TString angleDevBranchClass = "waisAngle/D:peakAngle:peakDisp:peakDev";

	waisPointingTree.Branch("waisPhiDev", & waisPhiDev, angleDevBranchClass);
	waisPointingTree.Branch("waisThetaDev", & waisThetaDev, angleDevBranchClass);

	//  This trees and its parent branches to the file's predominant data of polarity tests.
	TTree polarityTestsTree("polarityTestsTree", Form("TTree storing polarity test results of run %d.", runNum));
	//  Add the branches to the tree.
	TBranch * waveformTestsBranch = polarityTestsTree.Branch("waveformTests", & waveformTests);
	TBranch * deconvTestsBranch = polarityTestsTree.Branch("deconvTests", & deconvTests);
	TBranch * windowedTestsBranch = polarityTestsTree.Branch("windowedTests", & windowedTests);
	TBranch * windowedDeconvTestsBranch = polarityTestsTree.Branch("windowedDeconvTests", & windowedDeconvTests);

	//  Set up AnitaDataset object for given runNum.
	AnitaDataset d(runNum);
	
	//  UCorrelator stuff. Set up Analyzer objects.
	UCorrelator::AnalysisConfig cfg;
	cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution;  //  Rather than use the band-limited default, use all-pass to account for lack of elephant trunk corrections.
	cfg.response_option = UCorrelator::AnalysisConfig::ResponseTUFF;  //  Deconvolving accounting for TUFF filter response.
	UCorrelator::Analyzer analyzer(& cfg, true);  //  "true" is set so "analyzer" doesn't reset each time.
	FilterStrategy * fStrat = UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2");  //  Basic sine subtraction filter strategy.

	//  Fill polarityTestsTree with data.
	int count = 0;
	for (int entryNum = 0; entryNum < waisTree -> GetEntries(); ++entryNum) {
		
		waisTree -> GetEntry(entryNum);
		
		if (count >= 10) break;
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
			waisPhiDev.peakDisp = fmod(waisPhiDev.peakAngle - waisPhiDev.waisAngle, 360);
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
			fillPolarityClass(waveformTests.waveform, waveformTestsBranch -> FindBranch("waveform"), * waveformWais);

			TGraph * deconvWaveformWais = FFTtools::getInterpolatedGraph(analyzer.getDeconvolved(AnitaPol::kHorizontal, 0) -> even(), deltaT);
			but.filterGraph(deconvWaveformWais);
			fillPolarityClass(deconvTests.waveform, deconvTestsBranch -> FindBranch("waveform"), * deconvWaveformWais);

			//  Fill branches corresponding to negative derivative of entire event waveforms.
			TGraph negDerivWaveformWais = getNegDerivGraph(waveformWais);
			fillPolarityClass(waveformTests.negDerivWaveform, waveformTestsBranch -> FindBranch("negDerivWaveform"), negDerivWaveformWais);

			TGraph negDerivDeconvWaveformWais = getNegDerivGraph(deconvWaveformWais);
			fillPolarityClass(deconvTests.negDerivWaveform, deconvTestsBranch -> FindBranch("negDerivWaveform"), negDerivDeconvWaveformWais);

			//  Fill branches corresponding to antiderivative of entire event waveforms.
			TGraph antiderivWaveformWais = getAntiderivGraph(waveformWais);
			fillPolarityClass(waveformTests.antiderivWaveform, waveformTestsBranch -> FindBranch("antiderivWaveform"), antiderivWaveformWais);

			TGraph antiderivDeconvWaveformWais = getAntiderivGraph(deconvWaveformWais);
			fillPolarityClass(deconvTests.antiderivWaveform, deconvTestsBranch -> FindBranch("antiderivWaveform"), antiderivDeconvWaveformWais);

			//  Fill branches corresponding to correlation with entire event waveforms.
			TGraph corrWais = getCorrGraph(templateGraph, waveformWais);
			fillPolarityClass(waveformTests.corr, waveformTestsBranch -> FindBranch("corr"), corrWais);

			TGraph deconvCorrWais = getCorrGraph(deconvTemplateGraph, deconvWaveformWais);
			fillPolarityClass(deconvTests.corr, deconvTestsBranch -> FindBranch("corr"), deconvCorrWais);

			//  Fill branches corresponding to correlation with negative derivative of entire event waveforms.
			TGraph negDerivCorrWais = getCorrGraph(& negDerivTemplateGraph, & negDerivWaveformWais);
			fillPolarityClass(waveformTests.negDerivCorr, waveformTestsBranch -> FindBranch("negDerivCorr"), negDerivCorrWais);

			TGraph negDerivDeconvCorrWais = getCorrGraph(& negDerivDeconvTemplateGraph, & negDerivDeconvWaveformWais);
			fillPolarityClass(deconvTests.negDerivCorr, deconvTestsBranch -> FindBranch("negDerivCorr"), negDerivDeconvCorrWais);

			//  Fill branches corresponding to correlation with antiderivative of entire event waveforms.
			TGraph antiderivCorrWais = getCorrGraph(& antiderivTemplateGraph, & antiderivWaveformWais);
			fillPolarityClass(waveformTests.antiderivCorr, waveformTestsBranch -> FindBranch("antiderivCorr"), antiderivCorrWais);

			TGraph antiderivDeconvCorrWais = getCorrGraph(& antiderivDeconvTemplateGraph, & antiderivDeconvWaveformWais);
			fillPolarityClass(deconvTests.antiderivCorr, deconvTestsBranch -> FindBranch("antiderivCorr"), antiderivDeconvCorrWais);

			//  Create windowed graphs of event waveforms.
			//  Fill branches corresponding to windowed event waveform.
			TGraph windowedWaveformWais = getWindowedGraph(waveformWais, winFullWidth);
			fillPolarityClass(windowedTests.waveform, windowedTestsBranch -> FindBranch("waveform"), windowedWaveformWais);

			TGraph windowedDeconvWaveformWais = getWindowedGraph(deconvWaveformWais, winFullWidth);
			fillPolarityClass(windowedDeconvTests.waveform, windowedDeconvTestsBranch -> FindBranch("waveform"), windowedDeconvWaveformWais);

			//  Fill branches corresponding to negative derivative of windowed event waveforms.
			TGraph negDerivWindowedWaveformWais = getNegDerivGraph(& windowedWaveformWais);
			fillPolarityClass(windowedTests.negDerivWaveform, windowedTestsBranch -> FindBranch("negDerivWaveform"), negDerivWindowedWaveformWais);

			TGraph negDerivWindowedDeconvWaveformWais = getNegDerivGraph(& windowedDeconvWaveformWais);
			fillPolarityClass(windowedDeconvTests.negDerivWaveform, windowedDeconvTestsBranch -> FindBranch("negDerivWaveform"), negDerivWindowedDeconvWaveformWais);

			//  Fill branches corresponding to antiderivative of windowed event waveforms.
			TGraph antiderivWindowedWaveformWais = getAntiderivGraph(& windowedWaveformWais);
			fillPolarityClass(windowedTests.antiderivWaveform, windowedTestsBranch -> FindBranch("antiderivWaveform"), antiderivWindowedWaveformWais);

			TGraph antiderivWindowedDeconvWaveformWais = getAntiderivGraph(& windowedDeconvWaveformWais);
			fillPolarityClass(windowedDeconvTests.antiderivWaveform, windowedDeconvTestsBranch -> FindBranch("antiderivWaveform"), antiderivWindowedDeconvWaveformWais);

			//  Fill branches corresponding to correlation with windowed event waveforms.
			TGraph windowedCorrWais = getCorrGraph(& windowedTemplateGraph, & windowedWaveformWais);
			fillPolarityClass(windowedTests.corr, windowedTestsBranch -> FindBranch("corr"), windowedCorrWais);
			
			TGraph deconvWindowedCorrWais = getCorrGraph(& windowedDeconvTemplateGraph, & windowedDeconvWaveformWais);
			fillPolarityClass(windowedDeconvTests.corr, windowedDeconvTestsBranch -> FindBranch("corr"), deconvWindowedCorrWais);

			//  Fill branches corresponding to correlations with negative derivative of windowed event waveforms.
			TGraph negDerivWindowedCorrWais = getCorrGraph(& negDerivWindowedTemplateGraph, & negDerivWindowedWaveformWais);
			fillPolarityClass(windowedTests.negDerivCorr, windowedTestsBranch -> FindBranch("negDerivCorr"), negDerivWindowedCorrWais);

			TGraph negDerivWindowedDeconvCorrWais = getCorrGraph(& negDerivWindowedDeconvTemplateGraph, & negDerivWindowedDeconvWaveformWais);
			fillPolarityClass(windowedDeconvTests.negDerivCorr, windowedDeconvTestsBranch -> FindBranch("negDerivCorr"), negDerivWindowedDeconvCorrWais);

			//  Fill branches corresponding to correlations with antiderivative of windowed event waveforms.
			TGraph antiderivWindowedCorrWais = getCorrGraph(& antiderivWindowedTemplateGraph, & antiderivWindowedWaveformWais);
			fillPolarityClass(windowedTests.antiderivCorr, windowedTestsBranch -> FindBranch("antiderivCorr"), antiderivWindowedCorrWais);

			TGraph antiderivWindowedDeconvCorrWais = getCorrGraph(& antiderivWindowedDeconvTemplateGraph, & antiderivWindowedDeconvWaveformWais);
			fillPolarityClass(windowedDeconvTests.antiderivCorr, windowedDeconvTestsBranch -> FindBranch("antiderivCorr"), antiderivWindowedDeconvCorrWais);

			//  Filling the trees.
			waisPointingTree.Fill();
			polarityTestsTree.GetBranch("waveformTests") -> Fill();
			polarityTestsTree.GetBranch("deconvTests") -> Fill();
			polarityTestsTree.GetBranch("windowedTests") -> Fill();
			polarityTestsTree.GetBranch("windowedDeconvTests") -> Fill();
			polarityTestsTree.Fill();

			//  Further annoyance with the necessity of memory clearance.
			analyzer.clearInteractiveMemory();
			delete waveformWais, delete deconvWaveformWais;

			++count;
		}
	}
	
	//  Writing trees to polarityTestsFile.
	polarityTestsFile.cd();
	waisPointingTree.Write();
	polarityTestsTree.Write();
	
	//  Closing files.
	polarityTestsFile.Close();
	waisFile.Close();
	templateFile.Close();

	//  Freeing up memory.
	if (deltaT != 0.1) delete templateGraph, delete deconvTemplateGraph;
}
