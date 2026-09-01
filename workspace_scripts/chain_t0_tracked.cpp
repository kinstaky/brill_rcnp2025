void chain_t0_tracked() {
        TChain *tree = new TChain("tree");
        for (int run = 1024; run <= 1127; ++run) {
                if (run == 1034 || run == 1041 || (run>=1052 && run<=1055) || (run>=1061 && run<=1071) || run == 1086 || run == 1087 || run == 1090 || run == 1124) continue;
                tree->AddFile(TString::Format("telescope/t0_%04d.root", run));
        }
}
