/*  Modification of "collectPolarityTestsFile.C" for usage with ANITA-4 data. Started on
 *  11/21/17 by John Russell.
 */


void collectPolarityTestsFileA4() {

	//  Create file to store trees in.
	TFile polarityTestsFile("waisPolCutFilesA4/polarityTestsFileA4.root", "recreate");

	//  Chain up the trees.
//	TChain eventChain("eventTree");
	TChain waisPointingChain("waisPointingTree");
	TChain polarityTestsChain("polarityTestsTree");
	TChain deconvPolarityTestsChain("deconvPolarityTestsTree");
	TChain windowedPolarityTestsChain("windowedPolarityTestsTree");
	TChain windowedDeconvPolarityTestsChain("windowedDeconvPolarityTestsTree");

	for (int runNum = 123; runNum <= 152; ++runNum) {

		TString pathName = "waisPolCutFilesA4/template_1/template1PolarityTestsFileRun" + TString::Itoa(runNum, 10) + ".root";
//		eventChain.Add(pathName);
		waisPointingChain.Add(pathName);
		polarityTestsChain.Add(pathName);
		deconvPolarityTestsChain.Add(pathName);
		windowedPolarityTestsChain.Add(pathName);
		windowedDeconvPolarityTestsChain.Add(pathName);
	}

	//  Create trees to store chain entries in.
//	TTree * collectedEventTree = eventChain.CloneTree(-1, "fast");
//	collectedEventTree -> SetTitle("TTree storing event numbers.");

	TTree * collectedWaisPointingTree = waisPointingChain.CloneTree(-1, "fast");
	collectedWaisPointingTree -> SetTitle("TTree storing WAIS pointing results from interferometric map.");

	TTree * collectedPolarityTestsTree = polarityTestsChain.CloneTree(-1, "fast");
	collectedPolarityTestsTree -> SetTitle("TTree storing polarity test results of waveforms.");

	TTree * collectedDeconvPolarityTestsTree = deconvPolarityTestsChain.CloneTree(-1, "fast");
	collectedDeconvPolarityTestsTree -> SetTitle("TTree storing polarity test results of deconvolved waveforms.");

	TTree * collectedWindowedPolarityTestsTree = windowedPolarityTestsChain.CloneTree(-1, "fast");
	collectedWindowedPolarityTestsTree -> SetTitle("TTree storing polarity test results of windowed waveforms.");

	TTree * collectedWindowedDeconvPolarityTestsTree = windowedDeconvPolarityTestsChain.CloneTree(-1, "fast");
	collectedWindowedDeconvPolarityTestsTree -> SetTitle("TTree storing polarity test results of windowed, deconvolved waveforms.");

	//  Adding friends for convenience.
	collectedWaisPointingTree -> AddFriend("wais", "waisPolCutFilesA4/allWais_max_30005_sinsub_10_3_ad_2.root");
	collectedPolarityTestsTree -> AddFriend("waisPointingTree");
	collectedDeconvPolarityTestsTree -> AddFriend("waisPointingTree");
	collectedWindowedPolarityTestsTree -> AddFriend("waisPointingTree");
	collectedWindowedDeconvPolarityTestsTree -> AddFriend("waisPointingTree");

	//  Write to these trees, then delete them to disconnect the tree clones from the trees they were cloned from.
	polarityTestsFile.cd();
//	collectedEventTree -> Write();
	collectedWaisPointingTree -> Write();
	collectedPolarityTestsTree -> Write();
	collectedDeconvPolarityTestsTree -> Write();
	collectedWindowedPolarityTestsTree -> Write();
	collectedWindowedDeconvPolarityTestsTree -> Write();

//	delete collectedEventTree;
	delete collectedWaisPointingTree;
	delete collectedPolarityTestsTree;
	delete collectedDeconvPolarityTestsTree;
	delete collectedWindowedPolarityTestsTree;
	delete collectedWindowedDeconvPolarityTestsTree;

	//  Closing file.
	polarityTestsFile.Close();
}
