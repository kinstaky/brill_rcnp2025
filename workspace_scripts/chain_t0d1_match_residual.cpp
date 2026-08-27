void chain_t0d1_match_residual() {
	TH1F *hds = new TH1F("hds", "Distance of strips", 500, 0, 5000);
	TH1F *hdnns = new TH1F("hdnns", "Distance of normal-normal strips", 500, 0, 5000);
	TH1F *hdncs = new TH1F("hdncs", "Distance of normal-shared strips", 500, 0, 5000);
	TH1F *hdccs = new TH1F("hdccs", "Distance of shared-shared strips", 500, 0, 5000);
	TH1F *hdnls = new TH1F("hdnls", "Distance of normal-line strips", 500, 0, 5000);
	TH1F *hdnps = new TH1F("hdnps", "Distance of normal-plane strips", 500, 0, 5000);
	TH1F *hdcls = new TH1F("hdcls", "Distance of shared-line strips", 500, 0, 5000);
	TH1F *hdcps = new TH1F("hdcps", "Distance of shared-plane strips", 500, 0, 5000);
	TH1F *hdnscs = new TH1F("hdnscs", "Distance of normal-short_shared strips", 500, 0, 5000);
	TH1F *hdcscs = new TH1F("hdcscs", "Distance of shared-short_shared strips", 500, 0, 5000);

	TH1F *hfns = new TH1F("hfns", "Front residual hit distribution summary", 8, 0, 8);
	TH1F *hbns = new TH1F("hbns", "Back residual hit distribution summary", 8, 0, 8);
	TH1F *hrfss = new TH1F("hrfss", "Front residual strip distribution summary", 128, 0, 128);
	TH1F *hrbss = new TH1F("hrbss", "Back residual strip distribution summary", 128, 0, 128);
	TH2F *hrfsps = new TH2F("hrfsps", "Front residual strip pair distribution summary", 128, 0, 128, 128, 0, 128);
	TH2F *hrbsps = new TH2F("hrbsps", "Back residual strip pair distribution summary", 128, 0, 128, 128, 0, 128);

	for (int run = 1024; run <= 1127; ++run) {
		if (run == 1034 || run == 1041 || (run>=1052 && run<=1055) || (run>=1061 && run<=1071) || run == 1086 || run == 1087 || run == 1090 || run == 1124) continue;
		TFile ipf(TString::Format("match/t0d1_%04d.root", run));
		hds->Add((TH1F *)ipf.Get("hd"));
		hdnns->Add((TH1F *)ipf.Get("hdnn"));
		hdncs->Add((TH1F *)ipf.Get("hdnc"));
		hdccs->Add((TH1F *)ipf.Get("hdcc"));
		hdnls->Add((TH1F *)ipf.Get("hdnl"));
		hdnps->Add((TH1F *)ipf.Get("hdnp"));
		hdcls->Add((TH1F *)ipf.Get("hdcl"));
		hdcps->Add((TH1F *)ipf.Get("hdcp"));
		hdnscs->Add((TH1F *)ipf.Get("hdnsc"));
		hdcscs->Add((TH1F *)ipf.Get("hdcsc"));

		hfns->Add((TH1F *)ipf.Get("hrfn"));
		hbns->Add((TH1F *)ipf.Get("hrbn"));
		hrfss->Add((TH1F *)ipf.Get("hrfs"));
		hrbss->Add((TH1F *)ipf.Get("hrbs"));
		hrfsps->Add((TH2F *)ipf.Get("hrfsp"));
		hrbsps->Add((TH2F *)ipf.Get("hrbsp"));
		ipf.Close();
	}


	TChain *tree = new TChain("rtree");
	for (int run = 1024; run <= 1127; ++run) {
			if (run == 1034 || run == 1041 || (run>=1052 && run<=1055) || (run>=1061 && run<=1071) || run == 1086 || run == 1087 || run == 1090 || run == 1124) continue;
			tree->AddFile(TString::Format("match/t0d1_%04d.root", run));
	}
}
