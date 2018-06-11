/*  With 4 of the trees in "polarityTestsFileA4.root" structured in the same way, here we refactor
 *  the file so that those 4 trees are placed under a larger tree as branches. Started 12/7/17
 *  by John Russell.
 */

#include "TTreeCloner.h"

/*  Structure corresponding to polarity data from event signals.
 *  Max and min values. Also their sum, difference, and ratio of these (p).
 */
//struct polarity {
//
//	double maxTime, maxVal, minTime, minVal, maxMinMeanVal, maxMinHalfWidthVal, p, RMS, maxMinMeanSNR, maxMinHalfWidthSNR, sgndExtrTime, sgndExtrVal, sgndExtrSNR;
//};


/*  Shortcut function which adds branches to parent branch.
 */
//void addSubbranches(TBranch & brParent, TBranch & brIn[6], polarity & structIn[6]) {
//
//	TString waveformBranchStruct = "maxTime/D:maxVal:minTime:minVal:maxMinMeanVal:maxMinHalfWidthVal:p:RMS:maxMinMeanSNR:maxMinHalfWidthSNR:sgndExtrTime:sgndExtrVal:sgndExtrSNR";
//	brIn[0] = TBranch(& brParent, "waveform", & structIn[0], waveformBranchStruct);
//	brIn[1] = TBranch(& brParent, "negDerivWaveform", & structIn[1], waveformBranchStruct);
//	brIn[2] = TBranch(& brParent, "antiderivWaveform", & structIn[2], waveformBranchStruct);
//
//	TString corrBranchStruct = "maxTime/D:maxCorr:minTime:minCorr:maxMinMeanCorr:maxMinHalfWidthCorr:p:RMS:maxMinMeanSNR:maxMinHalfWidthSNR:sgndExtrTime:sgndExtrCorr:sgndExtrSNR";
//	brIn[3] = TBranch(& brParent, "corr", & structIn[3], corrBranchStruct);
//	brIn[4] = TBranch(& brParent, "negDerivCorr", & structIn[4], corrBranchStruct);
//	brIn[5] = TBranch(& brParent, "antiderivCorr", & structIn[5], corrBranchStruct);
//}


//  Function to copy branches from tree to a parent branch.
void branchCopy(TTree * trParent, TBranch & brParent) {

	brParent.GetListOfBranches() -> Add((TBranch *) trParent -> GetBranch("waveform") -> Clone());
	brParent.GetListOfBranches() -> Add((TBranch *) trParent -> GetBranch("negDerivWaveform") -> Clone());
	brParent.GetListOfBranches() -> Add((TBranch *) trParent -> GetBranch("antiderivWaveform") -> Clone());
	brParent.GetListOfBranches() -> Add((TBranch *) trParent -> GetBranch("corr") -> Clone());
	brParent.GetListOfBranches() -> Add((TBranch *) trParent -> GetBranch("negDerivCorr") -> Clone());
	brParent.GetListOfBranches() -> Add((TBranch *) trParent -> GetBranch("antiderivCorr") -> Clone());
//	brParent -> waveform = trParent -> waveform;
//	brParent -> negDerivWaveform = trParent -> negDerivWaveform;
//	brParent -> antiderivWaveform = trParent -> antiderivWaveform;
//	brParent -> corr = trParent -> corr;
//	brParent -> negDerivCorr = trParent -> negDerivCorr;
//	brParent -> antiderivCorr = trParent -> antiderivCorr;
//	fillPolarityBranch(trParent -> waveform, brIn[0]);
//	fillPolarityBranch(trParent -> negDerivWaveform, brIn[1]);
//	fillPolarityBranch(trParent -> antiderivWaveform, brIn[2]);
//	fillPolarityBranch(trParent -> corr, brIn[3]);
//	fillPolarityBranch(trParent -> negDerivCorr, brIn[4]);
//	fillPolarityBranch(trParent -> antiderivCorr, brIn[5]);
}


///*  Shortcut function to fill in the branch corresponding to a tree.
// */
//void fillPolarityBranch(TBranch * brParent, TBranch & brIn) {
//
//	brIn.maxTime = brParent -> maxTime;
//	brIn.maxVal = brParent -> maxVal;
//	brIn.minTime = brParent -> minTime;
//	brIn.minVal = brParent -> minVal;
//	brIn.maxMinMeanVal = brParent -> maxMinMeanVal;
//	brIn.maxMinHalfWidthVal = brParent -> maxMinHalfWidthVal;
//	brIn.p = brParent -> p;
//	brIn.RMS = brParent -> RMS;
//	brIn.maxMinMeanSNR = brParent -> maxMinMeanSNR;
//	brIn.maxMinHalfWidthSNR = brParent -> maxMinHalfWidthSNR;
//	brIn.sgndExtrTime = brParent -> sgndExtrTime;
//	brIn.sgndExtrVal = brParent -> sgndExtrVal;
//	brIn.sgndExtrSNR = brParent -> sgndExtrSNR;
//}


void refactorPolarityTestsFileA4() {

	//  Open original format file and access its trees.
	TFile polarityTestsFile("waisPolCutFilesA4/polarityTestsFileA4.root");
	TTree * eventTree = (TTree *) polarityTestsFile.Get("eventTree");
	TTree * waisPointingTree = (TTree *) polarityTestsFile.Get("waisPointingTree");
	TTree * polarityTestsTree = (TTree *) polarityTestsFile.Get("polarityTestsTree");
	TTree * deconvPolarityTestsTree = (TTree *) polarityTestsFile.Get("deconvPolarityTestsTree");
	TTree * windowedPolarityTestsTree = (TTree *) polarityTestsFile.Get("windowedPolarityTestsTree");
	TTree * windowedDeconvPolarityTestsTree = (TTree *) polarityTestsFile.Get("windowedDeconvPolarityTestsTree");

	//  Create file to store trees in, then new tree to store 4 trees.
	TFile reformattedFile("waisPolCutFilesA4/polarityTestsFileA4reformatted.root", "recreate");
	TTree reformattedTree("polarityTestsTree", "TTree storing polarity test results.");
	TBranch waveformTests, deconvTests, windowedTests, windowedDeconvTests;
	reformattedTree.Branch("waveformTests.", & waveformTests);
	reformattedTree.Branch("deconvTests.", & deconvTests);
	reformattedTree.Branch("windowedTests.", & windowedTests);
	reformattedTree.Branch("windowedDeconvTests.", & windowedDeconvTests);

	//  Structure the branches of the reformatted tree.
//	TBranch waveformBranch[6], deconvBranch[6], windowedBranch[6], windowedDeconvBranch[6];
//	polarity waveformStruct[6], deconvStruct[6], windowedStruct[6], windowedDeconvStruct[6];
//	addSubbranches(waveformTests, waveformBranch, waveformStruct);
//	addSubbranches(deconvTests, deconvBranch, deconvStruct);
//	addSubbranches(windowedTests, windowedBranch, windowedStruct);
//	addSubbbranches(windowedDeconvTests, windowedDeconvBranch, windowedDeconvStruct);

	//  Now, copy the branches.
	CollectBranches(polarityTestsTree -> GetListOfBranches(), waveformTests.GetListOfBranches());
	reformattedTree.CollectBranches(deconvPolarityTestsTree -> GetListOfBranches(), deconvTests.GetListOfBranches());
	reformattedTree.CollectBranches(windowedPolarityTestsTree -> GetListOfBranches(), windowedTests.GetListOfBranches());
	reformattedTree.CollectBranches(windowedDeconvPolarityTestsTree -> GetListOfBranches(), windowedDeconvTests.GetListOfBranches());
//	waveformTests.GetListOfBranches() -> Add((TObjArray *) polarityTestsTree -> GetListOfBranches() -> Clone());
//	deconvTests.GetListOfBranches() -> Add((TObjArray *) deconvPolarityTestsTree -> GetListOfBranches() -> Clone());
//	windowedTests.GetListOfBranches() -> Add((TObjArray *) windowedPolarityTestsTree -> GetListOfBranches() -> Clone());
//	windowedDeconvTests.GetListOfBranches() -> Add((TObjArray *) windowedDeconvPolarityTestsTree -> GetListOfBranches() -> Clone());
//	branchCopy(polarityTestsTree, waveformTests);
//	branchCopy(deconvPolarityTestsTree, deconvTests);
//	branchCopy(windowedPolarityTestsTree, windowedTests);
//	branchCopy(windowedDeconvPolarityTestsTree, windowedDeconvTests);

	//  Now, to write to this new tree.
//	for (int entryNum = 0; entryNum < waisPointingTree -> GetEntries(); ++entryNum) {
//
//		polarityTestsTree -> GetEntry(entryNum);
//		deconvPolarityTestsTree -> GetEntry(entryNum);
//		windowedPolarityTestsTree -> GetEntry(entryNum);
//		windowedDeconvPolarityTestsTree -> GetEntry(entryNum);
//
////		reformattedTree.Fill();
//		waveformTests.Fill();
//		deconvTests.Fill();
//		windowedTests.Fill();
//		windowedDeconvTests.Fill();
//	}

	//  Write to the relevant trees to the reformatted file.
	reformattedFile.cd();
	eventTree -> CloneTree() -> Write();
	waisPointingTree -> CloneTree() -> Write();
	reformattedTree.Write();

	//  Closing files.
	reformattedFile.Close();
	polarityTestsFile.Close();
}
