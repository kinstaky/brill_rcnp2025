#include <iostream>
#include <set>

#include <TFile.h>
#include <TChain.h>
#include <TF1.h>
#include <TF2.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TH1F.h>

#include <Eigen/Dense>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/utils.h"
#include "include/t0/dssd.h"

#include "include/normalize_extra.h"

inline bool AreStrips(
	const int *strip,
	const int first,
	const int second
) {
	return (
		(strip[0] == first && strip[1] == second)
		|| (strip[0] == second && strip[1] == first)
	);
}

inline void GetEnergy(
	const int *strip,
	const double *energy,
	double &first_energy,
	double &second_energy
) {
	first_energy = strip[0] < strip[1] ? energy[0] : energy[1];
	second_energy = strip[0] < strip[1] ? energy[1] : energy[0];
}


struct PointD5 {
	double value[5];
};

int PrincipalComponentAnalysisDim5(
	const std::vector<PointD5> &points,
	Eigen::Matrix<double, 5, 1> &eigen_values,
	Eigen::Matrix<double, 5, 1> &mean,
	Eigen::Matrix<double, 5, 1> &direction
) {
	// PCA
	int n = points.size();
	if (n < 3) return -1;
	Eigen::MatrixXd matrix(n, 5);
	// fill matrix
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < 5; ++j) {
			matrix(i, j) = points[i].value[j];
		}
	}
	// center data
	Eigen::Matrix<double, 1, 5> center = matrix.colwise().mean();
	Eigen::MatrixXd centered = matrix.rowwise() - center;
	mean = center.transpose();
	// covariance matrix
	Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(n-1);
	// eigen decomposition
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 5, 5>> solver;
	solver.compute(cov);
	if (solver.info() != Eigen::Success) return -2;
	eigen_values = solver.eigenvalues();
	Eigen::Matrix<double, 5, 5> eigen_vectors = solver.eigenvectors();
	direction = eigen_vectors.col(4);
	return 0;
}


void LineNormalize(
	TChain &chain,
	brill::DssdEvent &raw_event,
	brill::DssdNormalizeParameters &parameters,
	brill::T0D2ExtraNormalizeParameters &extra_parameters
) {
	// for convenient
	const int &fn = raw_event.front_num;
	const int &bn = raw_event.back_num;
	const int *fs = raw_event.front_strip;
	const int *bs = raw_event.back_strip;
	const double *rfe = raw_event.front_integral;
	const double *rbe = raw_event.back_integral;

	// energy difference of first and second strip
	TH1F hbs32a("hbs32a", "energy difference of back strip 31 and 33", 100, -2000, 2000);
	TH1F hbs32b("hbs32b", "energy difference of back strip 31 and 34", 100, -2000, 2000);
	TH1F hbs32c("hbs32c", "energy difference of back strip 34 and 31", 100, -2000, 2000);
	TH1F hbs36("hbs36", "energy difference of back strip 35 and 37", 100, -2000, 2000);
	TH1F hbs48("hbs48", "energy difference of back strip 47 and 49", 100, -2000, 2000);
	TH1F hbs68("hbs68", "energy difference of back strip 67 and 69", 100, -2000, 2000);
	TH1F hfs70("hfs70", "energy difference of front strip 69 and 71", 100, -2000, 2000);
	TH1F hfs1819("hfs1819", "energy difference of front strip 17 and 20", 100, -2000, 2000);
	TH1F hfs2("hfs2", "energy difference of front strip 1 and 3", 100, -2000, 2000);
	// no energy cut
	TGraph2D gbs32n[4], gbs36n, gbs48n, gbs68n, gfs70n, gfs1819n, gfs2n;
	// with energy cut on ajdacent strips differences
	TGraph2D gbs32[3], gbs36, gbs48, gbs68, gfs70, gfs1819, gfs2;
	// the thrid strip
	TGraph grbs30, grbs35, grbs50, grbs70, grfs22, grfs15;
	// the third and fourth strip
	TGraph2D grbs3438, grbs5046, grbs7066, grfs7273, grfs34;
	// 2/3/4 hits 3D energy for comparing
	TGraph2D gbs68_hit[3], gfs70_hit[3];
	// normal strip for comparing
	TGraph2D gbs4041, gbs6566, gbs9293, gfs5253;
	// points for PCA
	std::vector<Eigen::Vector3d> pbs32[3], pbs36, pbs48, pbs68, pfs70, pfs1819, pfs2;
	// points of regular strip for comparing
	std::vector<Eigen::Vector3d> pbs4041, pbs6566, pbs9293, pfs5253;
	// estimate small signal's linearity
	std::vector<PointD5> fs70_points;

	long long total = chain.GetEntries();
	long long last_percentage = 0;
	printf("Normalizing line in t0d2   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		if (entry * 100 / total > last_percentage) {
			last_percentage = entry * 100 / total;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}
		chain.GetEntry(entry);
		if (fn == 1 && bn >= 2) {
			double fe = NormEnergy(parameters, 0, fs[0], rfe[0]);
			if (bn == 2 && AreStrips(bs, 31, 33)) {
				double be31, be33;
				GetEnergy(bs, rbe, be31, be33);
				gbs32n[0].AddPoint(be31, be33, fe);
				hbs32a.Fill(be31-be33);
				if (be31-be33 > -1000.0 && be31-be33 < 100.0) {
					gbs32[0].AddPoint(be31, be33, fe);
					pbs32[0].push_back(Eigen::Vector3d(be31, be33, fe));
				}
			} else if (bs[0] == 31 && bs[1] == 34) {
				gbs32n[1].AddPoint(rbe[0], rbe[1], fe);
				gbs32n[3].AddPoint(rbe[0], rbe[1], fe);
				hbs32b.Fill(rbe[0]-rbe[1]*2.3);
				if (rbe[0]-rbe[1]*2.3 > 200.0 && rbe[0]-rbe[1]*2.3 < 1200.0) {
					gbs32[1].AddPoint(rbe[0], rbe[1], fe);
					pbs32[1].push_back(Eigen::Vector3d(rbe[0], rbe[1], fe));
				}
				if (bn == 3) grbs30.AddPoint(rbe[0], rbe[2]);
			} else if (bs[0] == 34 && bs[1] == 31) {
				gbs32n[2].AddPoint(rbe[0], rbe[1], fe);
				gbs32n[3].AddPoint(rbe[1], rbe[0], fe);
				hbs32c.Fill(rbe[0]-rbe[1]*2.5);
				if (rbe[0]-rbe[1]*2.5 > -100.0 && rbe[0]-rbe[1]*2.5 < 600.0) {
					gbs32[2].AddPoint(rbe[0], rbe[1], fe);
					pbs32[2].push_back(Eigen::Vector3d(rbe[0], rbe[1], fe));
				}
				if (bn == 3) grbs35.AddPoint(rbe[0], rbe[2]);
			} else if (AreStrips(bs, 35, 37)) {
				double be35, be37;
				GetEnergy(bs, rbe, be35, be37);
				gbs36n.AddPoint(be35, be37, fe);
				hbs36.Fill(be35-be37);
				if (be35-be37 > -600.0 && be35-be37 < 0.0) {
					gbs36.AddPoint(be35, be37, fe);
					pbs36.push_back(Eigen::Vector3d(be35, be37, fe));
				}
				if (bn == 4) grbs3438.AddPoint(rbe[2], rbe[3], be35);
			} else if (AreStrips(bs, 47, 49)) {
				double be47, be49;
				GetEnergy(bs, rbe, be47, be49);
				gbs48n.AddPoint(be47, be49, fe);
				hbs48.Fill(be47-be49);
				if (be47-be49 > -1100.0 && be47-be49 < 0.0) {
					gbs48.AddPoint(be47, be49, fe);
					pbs48.push_back(Eigen::Vector3d(be47, be49, fe));
				}
				if (bn == 3 && bs[2] == 50) grbs50.AddPoint(be47, rbe[2]);
				else if (bn == 4) grbs5046.AddPoint(rbe[2], rbe[3], be47);
			} else if (AreStrips(bs, 67, 69)) {
				double be67, be69;
				GetEnergy(bs, rbe, be67, be69);
				gbs68n.AddPoint(be67, be69, fe);
				hbs68.Fill(be67-be69);
				if (be67-be69 > -500.0 && be67-be69 < 0.0) {
					gbs68.AddPoint(be67, be69, fe);
					pbs68.push_back(Eigen::Vector3d(be67, be69, fe));
				}
				if (bn == 2) {
					gbs68_hit[0].AddPoint(be67, be69, fe);
				} else if (bn == 3) {
					gbs68_hit[1].AddPoint(be67, be69, fe);
					if (bs[2] == 70) grbs70.AddPoint(be67, rbe[2]);
				} else if (bn == 4) {
					gbs68_hit[2].AddPoint(be67, be69, fe);
					grbs7066.AddPoint(rbe[2], rbe[3], be67);
				}
			} else if (bn == 2 && AreStrips(bs, 40, 41)) {
				double be40, be41;
				GetEnergy(bs, rbe, be40, be41);
				gbs4041.AddPoint(be40, be41, fe);
				pbs4041.push_back(Eigen::Vector3d(be40, be41, fe));
			} else if (bn == 2 && AreStrips(bs, 65, 66)) {
				double be65, be66;
				GetEnergy(bs, rbe, be65, be66);
				gbs6566.AddPoint(be65, be66, fe);
				pbs6566.push_back(Eigen::Vector3d(be65, be66, fe));
			} else if (bn == 2 && AreStrips(bs, 92, 93)) {
				double be92, be93;
				GetEnergy(bs, rbe, be92, be93);
				gbs9293.AddPoint(be92, be93, fe);
				pbs9293.push_back(Eigen::Vector3d(be92, be93, fe));
			}
		} else if (fn >= 2 && bn == 1) {
			double be = NormEnergy(parameters, 1, bs[0], rbe[0]);
			if (AreStrips(fs, 69, 71)) {
				double fe69, fe71;
				GetEnergy(fs, rfe, fe69, fe71);
				gfs70n.AddPoint(fe69, fe71, be);
				hfs70.Fill(fe69-fe71);
				if (fe69-fe71 > -500.0 && fe69-fe71 < 0.0) {
					gfs70.AddPoint(fe69, fe71, be);
					pfs70.push_back(Eigen::Vector3d(fe69, fe71, be));
				}
				if (fn == 2) {
					gfs70_hit[0].AddPoint(fe69, fe71, be);
				} else if (fn == 3) {
					gfs70_hit[1].AddPoint(fe69, fe71, be);
				} else if (fn == 4) {
					gfs70_hit[2].AddPoint(fe69, fe71, be);
					grfs7273.AddPoint(rfe[2], rfe[3], fe69);
					if (fe69-fe71 > -500.0 && fe69-fe71 < 0.0) {
						PointD5 point;
						point.value[0] = fe69;
						point.value[1] = fe71;
						point.value[2] = rfe[2];
						point.value[3] = rfe[3];
						point.value[4] = be;
						fs70_points.push_back(point);
					}
				}
			} else if (AreStrips(fs, 17, 20)) {
				double fe17, fe20;
				GetEnergy(fs, rfe, fe17, fe20);
				gfs1819n.AddPoint(fe17, fe20, be);
				hfs1819.Fill(fe17-fe20*3.9);
				if (fabs(fe17-fe20*3.9) < 500.0) {
					gfs1819.AddPoint(fe17, fe20, be);
					pfs1819.push_back(Eigen::Vector3d(fe17, fe20, be));
				}
				if (fn == 3) {
					if (fs[2] == 15) grfs15.AddPoint(fe17, rfe[2]);
					else if (fs[2] == 22) grfs22.AddPoint(fe20, rfe[2]);
				}
			} else if (AreStrips(fs, 1, 3)) {
				double fe1, fe3;
				GetEnergy(fs, rfe, fe1, fe3);
				gfs2n.AddPoint(fe1, fe3, be);
				hfs2.Fill(fe1-fe3*1.5);
				if (fabs(fe1-fe3*1.5) < 200.0) {
					gfs2.AddPoint(fe1, fe3, be);
					pfs2.push_back(Eigen::Vector3d(fe1, fe3, be));
				}
			} else if (fn == 2 && AreStrips(fs, 52, 53)) {
				double fe52, fe53;
				GetEnergy(fs, rfe, fe52, fe53);
				gfs5253.AddPoint(fe52, fe53, be);
				pfs5253.push_back(Eigen::Vector3d(fe52, fe53, be));
			}
		}
	}
	printf("\b\b\b\b100%%\n");

	// pca for special strips
	const std::vector<Eigen::Vector3d>* pca_points[] = {
		pbs32, pbs32+1, pbs32+2, &pbs36, &pbs48, &pbs68, &pfs70, &pfs1819, &pfs2
	};
	const std::string pca_descriptions[] = {
		"Back strip 31 and 33",
		"Back strip 31 and 34",
		"Back strip 34 and 31",
		"Back strip 35 and 37",
		"Back strip 47 and 49",
		"Back strip 67 and 69",
		"Front strip 69 and 71",
		"Front strip 17 and 20",
		"Front strip 1 and 3",
	};
	const std::string pca_names[] = {
		"bs32a", "bs32b", "bs32c", "bs36", "bs48", "bs68",
		"fs70", "fs1819", "fs2",
	};
	for (int i = 0; i < 9; ++i) {
		PCAPrint(
			*pca_points[i],
			pca_descriptions[i],
			pca_names[i],
			true,
			extra_parameters
		);
	}
	std::cout << "---------------------------------------------------------\n";

	// try 5-dim PCA for front strip 69, 71, 72, 73 + back strip
	Eigen::Matrix<double, 5, 1> fs70_eigen_values, fs70_mean, fs70_direction;
	int result = PrincipalComponentAnalysisDim5(
		fs70_points, fs70_eigen_values, fs70_mean, fs70_direction
	);
	if (result == -1) {
		std::cerr << "Error: Too few points for 5D PCA fs 70.\n";
	} else if (result == -2) {
		std::cerr << "Error: 5D PCA fs 70 failed.\n";
	} else {
		double linearity = fs70_eigen_values(4) / fs70_eigen_values.sum();
		std::cout << "5D front strip 69, 71, 72 and 73:\n"
			<< "  Eigen values: " << fs70_eigen_values.transpose() << "\n"
			<< "  Linearity: " << linearity << "\n"
			<< "  Mean: " << fs70_mean.transpose() << "\n"
			<< "  Direction: " << fs70_direction.transpose() << "\n";
	}
	std::cout << "---------------------------------------------------------\n";

	hbs32a.Write();
	hbs32b.Write();
	hbs32c.Write();
	hbs36.Write();
	hbs48.Write();
	hbs68.Write();
	hfs70.Write();
	hfs1819.Write();
	hfs2.Write();
	gbs32n[0].Write("gbs32na");
	gbs32n[1].Write("gbs32nb");
	gbs32n[2].Write("gbs32nc");
	gbs32n[3].Write("gbs32nbc");
	gbs36n.Write("gbs36n");
	gbs48n.Write("gbs48n");
	gbs68n.Write("gbs68n");
	gfs70n.Write("gfs70n");
	gfs1819n.Write("gfs1819n");
	gfs2n.Write("gfs2n");
	gbs32[0].Write("gbs32a");
	gbs32[1].Write("gbs32b");
	gbs32[2].Write("gbs32c");
	gbs36.Write("gbs36");
	gbs48.Write("gbs48");
	gbs68.Write("gbs68");
	gfs70.Write("gfs70");
	gfs1819.Write("gfs1819");
	gfs2.Write("gfs2");
	grbs30.Write("grbs30");
	grbs35.Write("grbs35");
	grbs50.Write("grbs50");
	grbs70.Write("grbs70");
	grfs22.Write("grfs22");
	grfs15.Write("grfs15");
	grbs3438.Write("grbs3438");
	grbs5046.Write("grbs5046");
	grbs7066.Write("grbs7066");
	grfs7273.Write("grfs7273");
	gbs6566.Write("gnbs6566");
	gbs4041.Write("gnbs4041");
	gbs9293.Write("gnbs9293");
	gfs5253.Write("gnfs5253");
	gbs68_hit[0].Write("gbs68h2");
	gbs68_hit[1].Write("gbs68h3");
	gbs68_hit[2].Write("gbs68h4");
	gfs70_hit[0].Write("gfs70h2");
	gfs70_hit[1].Write("gfs70h3");
	gfs70_hit[2].Write("gfs70h4");
}


void LineEstimate(
	TChain &chain,
	brill::DssdEvent &raw_event,
	brill::DssdNormalizeParameters &parameters,
	brill::T0D2ExtraNormalizeParameters &extra
) {
	// for convenient
	const int &fn = raw_event.front_num;
	const int &bn = raw_event.back_num;
	const int *fs = raw_event.front_strip;
	const int *bs = raw_event.back_strip;
	const double *rfe = raw_event.front_integral;
	const double *rbe = raw_event.back_integral;

	// point-line difference (energy) after PCA
	TH1F hdbs32a("hdbs32a", "Back strip 31 and 33 difference", 200, 0, 2000);
	TH1F hdbs32b("hdbs32b", "Back strip 31 and 34 difference", 200, 0, 2000);
	TH1F hdbs32c("hdbs32c", "Back strip 34 and 31 difference", 200, 0, 2000);
	TH1F hdbs36("hdbs36", "Back strip 35 and 37 difference", 200, 0, 2000);
	TH1F hdbs48("hdbs48", "Back strip 47 and 49 difference", 200, 0, 2000);
	TH1F hdbs68("hdbs68", "Back strip 67 and 69 difference", 200, 0, 2000);
	TH1F hdfs70("hdfs70", "Front strip 69 and 71 difference", 200, 0, 2000);
	TH1F hdfs1819("hdfs1819", "Front strip 17 and 20 difference", 200, 0, 2000);
	TH1F hdfs2("hdfs2", "Front strip 1 and 3 difference", 200, 0, 2000);

	long long total = chain.GetEntries();
	long long last_percentage = 0;
	printf("Estimating line in t0d2   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		if (entry * 100 / total > last_percentage) {
			last_percentage = entry * 100 / total;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}
		chain.GetEntry(entry);
		if (fn == 1 && bn >= 2) {
			double fe = NormEnergy(parameters, 0, fs[0], rfe[0]);
			if (bn == 2 && AreStrips(bs, 31, 33)) {
				double be31, be33;
				GetEnergy(bs, rbe, be31, be33);
				Eigen::Vector3d point(be31, be33, fe);
				hdbs32a.Fill(brill::PointLineDistance(
					point,
					extra.pca["bs32a"].mean,
					extra.pca["bs32a"].direction
				));
			} else if (bs[0] == 31 && bs[1] == 34) {
				double be31 = rbe[0];
				double be34 = rbe[1];
				Eigen::Vector3d point(be31, be34, fe);
				hdbs32b.Fill(brill::PointLineDistance(
					point,
					extra.pca["bs32b"].mean,
					extra.pca["bs32b"].direction
				));
			} else if (bs[0] == 34 && bs[1] == 31) {
				double be34 = rbe[0];
				double be31 = rbe[1];
				Eigen::Vector3d point(be34, be31, fe);
				hdbs32c.Fill(brill::PointLineDistance(
					point,
					extra.pca["bs32c"].mean,
					extra.pca["bs32c"].direction
				));
			} else if (AreStrips(bs, 35, 37)) {
				double be35, be37;
				GetEnergy(bs, rbe, be35, be37);
				Eigen::Vector3d point(be35, be37, fe);
				hdbs36.Fill(brill::PointLineDistance(
					point,
					extra.pca["bs36"].mean,
					extra.pca["bs36"].direction
				));
			} else if (AreStrips(bs, 47, 49)) {
				double be47, be49;
				GetEnergy(bs, rbe, be47, be49);
				Eigen::Vector3d point(be47, be49, fe);
				hdbs48.Fill(brill::PointLineDistance(
					point,
					extra.pca["bs48"].mean,
					extra.pca["bs48"].direction
				));
			} else if (AreStrips(bs, 67, 69)) {
				double be67, be69;
				GetEnergy(bs, rbe, be67, be69);
				Eigen::Vector3d point(be67, be69, fe);
				hdbs68.Fill(brill::PointLineDistance(
					point,
					extra.pca["bs68"].mean,
					extra.pca["bs68"].direction
				));
			}
		} else if (fn >= 2 && bn == 1) {
			double be = NormEnergy(parameters, 1, bs[0], rbe[0]);
			if (AreStrips(fs, 69, 71)) {
				double fe69, fe71;
				GetEnergy(fs, rfe, fe69, fe71);
				Eigen::Vector3d point(fe69, fe71, be);
				hdfs70.Fill(brill::PointLineDistance(
					point,
					extra.pca["fs70"].mean,
					extra.pca["fs70"].direction
				));
			} else if (AreStrips(fs, 17, 20)) {
				double fe17, fe20;
				GetEnergy(fs, rfe, fe17, fe20);
				Eigen::Vector3d point(fe17, fe20, be);
				hdfs1819.Fill(brill::PointLineDistance(
					point,
					extra.pca["fs1819"].mean,
					extra.pca["fs1819"].direction
				));
			} else if (AreStrips(fs, 1, 3)) {
				double fe1, fe3;
				GetEnergy(fs, rfe, fe1, fe3);
				Eigen::Vector3d point(fe1, fe3, be);
				hdfs2.Fill(brill::PointLineDistance(
					point,
					extra.pca["fs2"].mean,
					extra.pca["fs2"].direction
				));
			}
		}
	}
	printf("\b\b\b\b100%%\n");

	hdbs32a.Write();
	hdbs32b.Write();
	hdbs32c.Write();
	hdbs36.Write();
	hdbs48.Write();
	hdbs68.Write();
	hdfs70.Write();
	hdfs1819.Write();
	hdfs2.Write();
}


void PlaneNormalize(
	TChain &chain,
	brill::DssdEvent &raw_event,
	brill::DssdNormalizeParameters &parameters,
	brill::T0D2ExtraNormalizeParameters &extra
) {
	// for convenient
	const int &fn = raw_event.front_num;
	const int &bn = raw_event.back_num;
	const int *fs = raw_event.front_strip;
	const int *bs = raw_event.back_strip;
	const double *rfe = raw_event.front_integral;
	const double *rbe = raw_event.back_integral;

	TH1F habs32a("habs32a", "angle to line of back strip 31 and 33", 100, -1.0, 1.0);
	TH1F habs32b("habs32b", "angle to line of back strip 31 and 34", 100, -1.0, 1.0);
	TH1F habs32c("habs32c", "angle to line of back strip 34 and 31", 100, -1.0, 1.0);
	TH1F habs36("habs36", "angle to line of back strip 35 and 37", 100, -1.0, 1.0);
	TH1F habs48("habs48", "angle to line of back strip 47 and 49", 100, -1.0, 1.0);
	TH1F habs68("habs68", "angle to line of back strip 67 and 69", 100, -1.0, 1.0);
	TH1F hafs70("hafs70", "angle to line of front strip 69 and 71", 100, -1.0, 1.0);
	TH1F hafs1819("hafs1819", "angle to line of front strip 17 and 20", 100, -1.0, 1.0);
	TH1F hafs2("hafs2", "angle to line of front strip 1 and 3", 100, -1.0, 1.0);
	TGraph2D gbs32ag[3], gfs2g[3];
	std::vector<Eigen::Vector3d> pbs32a[2], pbs32b[2], pbs32c[2], pbs36[2];
	std::vector<Eigen::Vector3d> pbs48[2], pbs68[2], pfs70[2], pfs1819[2], pfs2[2];

	long long total = chain.GetEntries();
	long long last_percentage = 0;
	printf("Normalize plane in t0d2   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		if (entry * 100 / total > last_percentage) {
			last_percentage = entry * 100 / total;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}
		chain.GetEntry(entry);
		if (bn >= 2 && fn == 1) {
			double fe = NormEnergy(parameters, 0, fs[0], rfe[0]);
			if (bn == 2 && AreStrips(bs, 31, 33)) {
				double be31, be33;
				GetEnergy(bs, rbe, be31, be33);
				Eigen::Vector3d point(be31, be33, fe);
				const Eigen::Vector3d &centroid = extra.pca["bs32a"].mean;
				const Eigen::Vector3d &direction = extra.pca["bs32a"].direction;
				double distance = brill::PointLineDistance(point, centroid, direction);
				if (distance < 400.0) continue;
				double angle = brill::PointVerticalAngle(point, centroid, direction);
				habs32a.Fill(angle);
				if (angle > 0.1 && angle < 0.25) gbs32ag[0].AddPoint(be31, be33, fe);
				else if (angle > 0.25 && angle < 0.4) gbs32ag[1].AddPoint(be31, be33, fe);
				else if (fabs(angle) > 0.9) gbs32ag[2].AddPoint(be31, be33, fe);
				if (angle > 0.1 && angle < 0.25) pbs32a[0].push_back(point);
				else if (angle > 0.25 && angle < 0.4) pbs32a[1].push_back(point);
			} else if (bs[0] == 31 && bs[1] == 34) {
				double be31 = rbe[0];
				double be34 = rbe[1];
				Eigen::Vector3d point(be31, be34, fe);
				const Eigen::Vector3d &centroid = extra.pca["bs32b"].mean;
				const Eigen::Vector3d &direction = extra.pca["bs32b"].direction;
				double distance = brill::PointLineDistance(point, centroid, direction);
				if (distance < 500.0) continue;
				double angle = brill::PointVerticalAngle(point, centroid, direction);
				habs32b.Fill(angle);
				if (angle > -0.5 && angle < -0.3) pbs32b[0].push_back(point);
				else if (angle > 0.6 && angle < 0.9) pbs32b[1].push_back(point);
			} else if (bs[0] == 34 && bs[1] == 31) {
				double be34 = rbe[0];
				double be31 = rbe[1];
				Eigen::Vector3d point(be34, be31, fe);
				const Eigen::Vector3d &centroid = extra.pca["bs32c"].mean;
				const Eigen::Vector3d &direction = extra.pca["bs32c"].direction;
				double distance = brill::PointLineDistance(point, centroid, direction);
				if (distance < 500.0) continue;
				double angle = brill::PointVerticalAngle(point, centroid, direction);
				habs32c.Fill(angle);
				if (angle > -0.9 && angle < -0.6) pbs32c[0].push_back(point);
				else if (angle > 0.9) pbs32c[1].push_back(point);
			} else if (AreStrips(bs, 35, 37)) {
				double be35, be37;
				GetEnergy(bs, rbe, be35, be37);
				Eigen::Vector3d point(be35, be37, fe);
				const Eigen::Vector3d &centroid = extra.pca["bs36"].mean;
				const Eigen::Vector3d &direction = extra.pca["bs36"].direction;
				double distance = brill::PointLineDistance(point, centroid, direction);
				if (distance < 600.0) continue;
				double angle = brill::PointVerticalAngle(point, centroid, direction);
				habs36.Fill(angle);
				if (angle > -0.3 && angle < 0.0) pbs36[0].push_back(point);
				else if (angle > 0.5 && angle < 0.8) pbs36[1].push_back(point);
			} else if (AreStrips(bs, 47, 49)) {
				double be47, be49;
				GetEnergy(bs, rbe, be47, be49);
				Eigen::Vector3d point(be47, be49, fe);
				const Eigen::Vector3d &centroid = extra.pca["bs48"].mean;
				const Eigen::Vector3d &direction = extra.pca["bs48"].direction;
				double distance = brill::PointLineDistance(point, centroid, direction);
				if (distance < 600.0) continue;
				double angle = brill::PointVerticalAngle(point, centroid, direction);
				habs48.Fill(angle);
				if (angle > -0.1 && angle < 0.1) pbs48[0].push_back(point);
				else if (angle > 0.4 && angle < 0.65) pbs48[1].push_back(point);
			} else if (AreStrips(bs, 67, 69)) {
				double be67, be69;
				GetEnergy(bs, rbe, be67, be69);
				Eigen::Vector3d point(be67, be69, fe);
				const Eigen::Vector3d &centroid = extra.pca["bs68"].mean;
				const Eigen::Vector3d &direction = extra.pca["bs68"].direction;
				double distance = brill::PointLineDistance(point, centroid, direction);
				if (distance < 500.0) continue;
				double angle = brill::PointVerticalAngle(point, centroid, direction);
				habs68.Fill(angle);
				if (angle > -0.4 && angle < -0.1) pbs68[0].push_back(point);
				else if (angle > 0.6 && angle < 0.9) pbs68[1].push_back(point);
			}
		} else if (fn >= 2 && bn == 1) {
			double be = brill::NormEnergy(parameters, 1, bs[0], rbe[0]);
			if (AreStrips(fs, 69, 71)) {
				double fe69, fe71;
				GetEnergy(fs, rfe, fe69, fe71);
				Eigen::Vector3d point(fe69, fe71, be);
				const Eigen::Vector3d &centroid = extra.pca["fs70"].mean;
				const Eigen::Vector3d &direction = extra.pca["fs70"].direction;
				double distance = brill::PointLineDistance(point, centroid, direction);
				if (distance < 300.0) continue;
				double angle = brill::PointVerticalAngle(point, centroid, direction);
				hafs70.Fill(angle);
				if (angle > -0.6 && angle < -0.4) pfs70[0].push_back(point);
				else if (angle > 0.9) pfs70[1].push_back(point);
			} else if (AreStrips(fs, 17, 20)) {
				double fe17, fe20;
				GetEnergy(fs, rfe, fe17, fe20);
				Eigen::Vector3d point(fe17, fe20, be);
				const Eigen::Vector3d &centroid = extra.pca["fs1819"].mean;
				const Eigen::Vector3d &direction = extra.pca["fs1819"].direction;
				double distance = brill::PointLineDistance(point, centroid, direction);
				if (distance < 300.0) continue;
				double angle = brill::PointVerticalAngle(point, centroid, direction);
				hafs1819.Fill(angle);
				if (angle > -0.9 && angle < -0.6) pfs1819[0].push_back(point);
				else if (angle > 0.9) pfs1819[1].push_back(point);
			} else if (AreStrips(fs, 1, 3)) {
				double fe1, fe3;
				GetEnergy(fs, rfe, fe1, fe3);
				Eigen::Vector3d point(fe1, fe3, be);
				const Eigen::Vector3d &centroid = extra.pca["fs2"].mean;
				const Eigen::Vector3d &direction = extra.pca["fs2"].direction;
				double distance = brill::PointLineDistance(point, centroid, direction);
				if (distance < 1000.0) continue;
				double angle = brill::PointVerticalAngle(point, centroid, direction);
				hafs2.Fill(angle);
				if (angle > -0.9 && angle < -0.75) gfs2g[0].AddPoint(fe1, fe3, be);
				else if (angle > -0.6 && angle < -0.4) gfs2g[1].AddPoint(fe1, fe3, be);
				else if (fabs(angle) > 0.9) gfs2g[2].AddPoint(fe1, fe3, be);
				if (angle > -0.9 && angle < -0.75) pfs2[0].push_back(point);
				else if (angle > -0.6 && angle < -0.4) pfs2[1].push_back(point);
			}
		}
	}
	printf("\b\b\b\b100%%\n");

	// PCA for plane
	const std::vector<Eigen::Vector3d>* pca_points[] = {
		pbs32a, pbs32a+1, pbs32b, pbs32b+1, pbs32c, pbs32c+1,
		pbs36, pbs36+1, pbs48, pbs48+1, pbs68, pbs68+1,
		pfs70, pfs70+1, pfs1819, pfs1819+1, pfs2, pfs2+1
	};
	const std::string pca_descriptions[] = {
		"Back strip 31 and 33 plane 0",
		"Back strip 31 and 33 plane 1",
		"Back strip 31 and 34 plane 0",
		"Back strip 31 and 34 plane 1",
		"Back strip 34 and 31 plane 0",
		"Back strip 34 and 31 plane 1",
		"Back strip 35 and 37 plane 0",
		"Back strip 35 and 37 plane 1",
		"Back strip 47 and 49 plane 0",
		"Back strip 47 and 49 plane 1",
		"Back strip 67 and 69 plane 0",
		"Back strip 67 and 69 plane 1",
		"Front strip 69 and 71 plane 0",
		"Front strip 69 and 71 plane 1",
		"Front strip 17 and 20 plane 0",
		"Front strip 17 and 20 plane 1",
		"Front strip 1 and 3 plane 0",
		"Front strip 1 and 3 plane 1"
	};
	const std::string pca_names[] = {
		"bs32ap0", "bs32ap1", "bs32bp0", "bs32bp1", "bs32cp0", "bs32cp1",
		"bs36p0", "bs36p1", "bs48p0", "bs48p1", "bs68p0", "bs68p1",
		"fs70p0", "fs70p1", "fs1819p0", "fs1819p1", "fs2p0", "fs2p1"
	};
	for (int i = 0; i < 18; ++i) {
		PCAPrint(
			*pca_points[i],
			pca_descriptions[i],
			pca_names[i],
			false,
			extra
		);
	}
	std::cout << "---------------------------------------------------------\n";

	habs32a.Write();
	habs32b.Write();
	habs32c.Write();
	habs36.Write();
	habs48.Write();
	habs68.Write();
	hafs70.Write();
	hafs1819.Write();
	hafs2.Write();
	gbs32ag[0].Write("gbs32aga0");
	gbs32ag[1].Write("gbs32aga1");
	gbs32ag[2].Write("gbs32aga2");
	gfs2g[0].Write("gfs2ga0");
	gfs2g[1].Write("gfs2ga1");
	gfs2g[2].Write("gfs2ga2");
}


void PlaneEstimate(
	TChain &chain,
	brill::DssdEvent &raw_event,
	brill::DssdNormalizeParameters &parameters,
	brill::T0D2ExtraNormalizeParameters &extra
) {
	// for convenient
	const int &fn = raw_event.front_num;
	const int &bn = raw_event.back_num;
	const int *fs = raw_event.front_strip;
	const int *bs = raw_event.back_strip;
	const double *rfe = raw_event.front_integral;
	const double *rbe = raw_event.back_integral;

	TH1F hdbs32a_p0("hdp0bs32a", "distance to plane 0 for back strip 31 and 33", 100, 0, 5000);
	TH1F hdbs32a_p1("hdp1bs32a", "distance to plane 1 for back strip 31 and 33", 100, 0, 5000);
	TH1F hdbs32b_p0("hdp0bs32b", "distance to plane 0 for back strip 31 and 34", 100, 0, 5000);
	TH1F hdbs32b_p1("hdp1bs32b", "distance to plane 1 for back strip 31 and 34", 100, 0, 5000);
	TH1F hdbs32c_p0("hdp0bs32c", "distance to plane 0 for back strip 34 and 31", 100, 0, 5000);
	TH1F hdbs32c_p1("hdp1bs32c", "distance to plane 1 for back strip 34 and 31", 100, 0, 5000);
	TH1F hdbs36_p0("hdp0bs36", "distance to plane 0 for back strip 35 and 37", 100, 0, 5000);
	TH1F hdbs36_p1("hdp1bs36", "distance to plane 1 for back strip 35 and 37", 100, 0, 5000);
	TH1F hdbs48_p0("hdp0bs48", "distance to plane 0 for back strip 47 and 49", 100, 0, 5000);
	TH1F hdbs48_p1("hdp1bs48", "distance to plane 1 for back strip 47 and 49", 100, 0, 5000);
	TH1F hdbs68_p0("hdp0bs68", "distance to plane 0 for back strip 67 and 69", 100, 0, 5000);
	TH1F hdbs68_p1("hdp1bs68", "distance to plane 1 for back strip 67 and 69", 100, 0, 5000);
	TH1F hdfs70_p0("hdp0fs70", "distance to plane 0 for front strip 69 and 71", 100, 0, 5000);
	TH1F hdfs70_p1("hdp1fs70", "distance to plane 1 for front strip 69 and 71", 100, 0, 5000);
	TH1F hdfs1819_p0("hdp0fs1819", "distance to plane 0 for front strip 17 and 20", 100, 0, 5000);
	TH1F hdfs1819_p1("hdp1fs1819", "distance to plane 1 for front strip 17 and 20", 100, 0, 5000);
	TH1F hdfs2_p0("hdp0fs2", "distance to plane 0 for front strip 1 and 3", 100, 0, 5000);
	TH1F hdfs2_p1("hdp1fs2", "distance to plane 1 for front strip 1 and 3", 100, 0, 5000);
	std::vector<Eigen::Vector3d> pbs32a[2], pbs32b[2], pbs32c[2], pbs36[2], pbs48[2], pbs68[2];
	std::vector<Eigen::Vector3d> pfs70[2], pfs1819[2], pfs2[2];

	long long total = chain.GetEntries();
	long long last_percentage = 0;
	printf("Estimating plane in t0d1   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		if (entry * 100 / total > last_percentage) {
			last_percentage = entry * 100 / total;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}
		chain.GetEntry(entry);
		if (bn >= 2 && fn == 1) {
			double fe = NormEnergy(parameters, 0, fs[0], rfe[0]);
			if (bn == 2 && AreStrips(bs, 31, 33)) {
				double be31, be33;
				GetEnergy(bs, rbe, be31, be33);
				Eigen::Vector3d point(be31, be33, fe);
				double distance_p0 = brill::PointPlaneDistance(
					point, extra.pca["bs32ap0"].mean, extra.pca["bs32ap0"].direction
				);
				double distance_p1 = brill::PointPlaneDistance(
					point, extra.pca["bs32ap1"].mean, extra.pca["bs32ap1"].direction
				);
				hdbs32a_p0.Fill(distance_p0);
				hdbs32a_p1.Fill(distance_p1);
			} else if (bs[0] == 31 && bs[1] == 34) {
				double be31 = rbe[0];
				double be34 = rbe[1];
				Eigen::Vector3d point(be31, be34, fe);
				double distance_p0 = brill::PointPlaneDistance(
					point, extra.pca["bs32bp0"].mean, extra.pca["bs32bp0"].direction
				);
				double distance_p1 = brill::PointPlaneDistance(
					point, extra.pca["bs32bp1"].mean, extra.pca["bs32bp1"].direction
				);
				hdbs32b_p0.Fill(distance_p0);
				hdbs32b_p1.Fill(distance_p1);
			} else if (bs[0] == 34 && bs[1] == 31) {
				double be34 = rbe[0];
				double be31 = rbe[1];
				Eigen::Vector3d point(be34, be31, fe);
				double distance_p0 = brill::PointPlaneDistance(
					point, extra.pca["bs32cp0"].mean, extra.pca["bs32cp0"].direction
				);
				double distance_p1 = brill::PointPlaneDistance(
					point, extra.pca["bs32cp1"].mean, extra.pca["bs32cp1"].direction
				);
				hdbs32c_p0.Fill(distance_p0);
				hdbs32c_p1.Fill(distance_p1);
			} else if (AreStrips(bs, 35, 37)) {
				double be35, be37;
				GetEnergy(bs, rbe, be35, be37);
				Eigen::Vector3d point(be35, be37, fe);
				double distance_p0 = brill::PointPlaneDistance(
					point, extra.pca["bs36p0"].mean, extra.pca["bs36p0"].direction
				);
				double distance_p1 = brill::PointPlaneDistance(
					point, extra.pca["bs36p1"].mean, extra.pca["bs36p1"].direction
				);
				hdbs36_p0.Fill(distance_p0);
				hdbs36_p1.Fill(distance_p1);
			} else if (AreStrips(bs, 47, 49)) {
				double be47, be49;
				GetEnergy(bs, rbe, be47, be49);
				Eigen::Vector3d point(be47, be49, fe);
				double distance_p0 = brill::PointPlaneDistance(
					point, extra.pca["bs48p0"].mean, extra.pca["bs48p0"].direction
				);
				double distance_p1 = brill::PointPlaneDistance(
					point, extra.pca["bs48p1"].mean, extra.pca["bs48p1"].direction
				);
				hdbs48_p0.Fill(distance_p0);
				hdbs48_p1.Fill(distance_p1);
			} else if (AreStrips(bs, 67, 69)) {
				double be67, be69;
				GetEnergy(bs, rbe, be67, be69);
				Eigen::Vector3d point(be67, be69, fe);
				double distance_p0 = brill::PointPlaneDistance(
					point, extra.pca["bs68p0"].mean, extra.pca["bs68p0"].direction
				);
				double distance_p1 = brill::PointPlaneDistance(
					point, extra.pca["bs68p1"].mean, extra.pca["bs68p1"].direction
				);
				hdbs68_p0.Fill(distance_p0);
				hdbs68_p1.Fill(distance_p1);
			}
		} else if (fn >= 2 && bn == 1) {
			double be = brill::NormEnergy(parameters, 1, bs[0], rbe[0]);
			if (AreStrips(fs, 69, 71)) {
				double fe69, fe71;
				GetEnergy(fs, rfe, fe69, fe71);
				Eigen::Vector3d point(fe69, fe71, be);
				double distance_p0 = brill::PointPlaneDistance(
					point, extra.pca["fs70p0"].mean, extra.pca["fs70p0"].direction
				);
				double distance_p1 = brill::PointPlaneDistance(
					point, extra.pca["fs70p1"].mean, extra.pca["fs70p1"].direction
				);
				hdfs70_p0.Fill(distance_p0);
				hdfs70_p1.Fill(distance_p1);
			} else if (AreStrips(fs, 17, 20)) {
				double fe17, fe20;
				GetEnergy(fs, rfe, fe17, fe20);
				Eigen::Vector3d point(fe17, fe20, be);
				double distance_p0 = brill::PointPlaneDistance(
					point, extra.pca["fs1819p0"].mean, extra.pca["fs1819p0"].direction
				);
				double distance_p1 = brill::PointPlaneDistance(
					point, extra.pca["fs1819p1"].mean, extra.pca["fs1819p1"].direction
				);
				hdfs1819_p0.Fill(distance_p0);
				hdfs1819_p1.Fill(distance_p1);
			} else if (AreStrips(fs, 1, 3)) {
				double fe1, fe3;
				GetEnergy(fs, rfe, fe1, fe3);
				Eigen::Vector3d point(fe1, fe3, be);
				double distance_p0 = brill::PointPlaneDistance(
					point, extra.pca["fs2p0"].mean, extra.pca["fs2p0"].direction
				);
				double distance_p1 = brill::PointPlaneDistance(
					point, extra.pca["fs2p1"].mean, extra.pca["fs2p1"].direction
				);
				hdfs2_p0.Fill(distance_p0);
				hdfs2_p1.Fill(distance_p1);
			}
		}
	}
	printf("\b\b\b\b100%%\n");

	// save plots
	hdbs32a_p0.Write();
	hdbs32a_p1.Write();
	hdbs32b_p0.Write();
	hdbs32b_p1.Write();
	hdbs32c_p0.Write();
	hdbs32c_p1.Write();
	hdbs36_p0.Write();
	hdbs36_p1.Write();
	hdbs48_p0.Write();
	hdbs48_p1.Write();
	hdbs68_p0.Write();
	hdbs68_p1.Write();
	hdfs70_p0.Write();
	hdfs70_p1.Write();
	hdfs1819_p0.Write();
	hdfs1819_p1.Write();
	hdfs2_p0.Write();
	hdfs2_p1.Write();
}


void ClassifyEvent(
	const Eigen::Vector3d &point,
	const brill::T0D2ExtraNormalizeParameters &extra,
	const std::string &name,
	const double line_threshold,
	const double plane0_threshold,
	const double plane1_threshold,
	double &line_distance,
	double &plane0_distance,
	double &plane1_distance,
	int &type,
	int &flag
) {
	auto line_parameter = extra.pca.at(name);
	line_distance = brill::PointLineDistance(
		point, line_parameter.mean, line_parameter.direction
	);
	if (line_distance < line_threshold) {
		type = 1;
		flag |= 1;
	} else {
		auto plane0_parameter = extra.pca.at(name + "p0");
		auto plane1_parameter = extra.pca.at(name + "p1");
		plane0_distance = brill::PointPlaneDistance(
			point, plane0_parameter.mean, plane0_parameter.direction
		);
		plane1_distance = brill::PointPlaneDistance(
			point, plane1_parameter.mean, plane1_parameter.direction
		);
		if (plane0_distance < plane0_threshold && plane1_distance < plane1_threshold) {
			type = plane0_distance < plane1_distance ? 2 : 3;
			flag |= 6;
		} else if (plane0_distance < plane0_threshold) {
			type = 2;
			flag |= 2;
		} else if (plane1_distance < plane1_threshold) {
			type = 3;
			flag |= 4;
		}
	}
}


void ReportClassificationResult(
	int *num,
	const std::string &description
) {
	std::cout << "Classified " << description << ": \n"
		<< "  line   : " << num[1] << " / " << num[0] << " = "
		<< std::fixed << std::setprecision(2) << num[1] * 100.0 / num[0] << "%\n"
		<< "  plane 0: " << num[2] << " / " << num[0] << " = "
		<< std::fixed << std::setprecision(2) << num[2] * 100.0 / num[0] << "%\n"
		<< "  plane 1: " << num[3] << " / " << num[0] << " = "
		<< std::fixed << std::setprecision(2) << num[3] * 100.0 / num[0] << "%\n"
		<< "  total  : " << num[1]+num[2]+num[3] << " / " << num[0] << " = "
		<< std::fixed << std::setprecision(2)
		<< (num[1]+num[2]+num[3]) * 100.0 / num[0] << "%\n";
}


void Classify(
	TChain &chain,
	brill::DssdEvent &raw_event,
	brill::DssdNormalizeParameters &parameters,
	brill::T0D2ExtraNormalizeParameters &extra
) {
	// for convenient
	const int &fn = raw_event.front_num;
	const int &bn = raw_event.back_num;
	const int *fs = raw_event.front_strip;
	const int *bs = raw_event.back_strip;
	const double *rfe = raw_event.front_integral;
	const double *rbe = raw_event.back_integral;

	TGraph2D gbs68[4];

	TTree opt("tree", "tree for T0D1 extra normalize classification");
	int side;
	int strip[2], opposite_strip;
	double energy[2], opposite_energy;
	double line_distance, plane_distance[2];
	// 0: not classified, 1: line, 2: plane 0, 3: plane 1
	int type;
	// bit0: line, bit1: plane 0, bit2: plane 1
	int flag;
	// setup branches
	opt.Branch("side", &side, "side/I");
	opt.Branch("strip", strip, "s[2]/I");
	opt.Branch("opposite_strip", &opposite_strip, "os/I");
	opt.Branch("energy", energy, "e[2]/D");
	opt.Branch("opposite_energy", &opposite_energy, "oe/D");
	opt.Branch("line_distance", &line_distance, "ld/D");
	opt.Branch("plane_distance", plane_distance, "pd[2]/D");
	opt.Branch("type", &type, "type/I");
	opt.Branch("flag", &flag, "flag/I");

	int num_bs32a[4], num_bs32b[4], num_bs32c[4];
	int num_bs36[4], num_bs48[4], num_bs68[4];
	int num_fs70[4], num_fs1819[4], num_fs2[4];
	for (int i = 0; i < 4; ++i) {
		num_bs32a[i] = num_bs32b[i] = num_bs32c[i] = 0;
		num_bs36[i] = num_bs48[i] = num_bs68[i] = 0;
		num_fs70[i] = num_fs1819[i] = num_fs2[i] = 0;
	}

	long long total = chain.GetEntries();
	long long last_percentage = 0;
	printf("Classifying t0d1   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		if (entry * 100 / total > last_percentage) {
			last_percentage = entry * 100 / total;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}
		chain.GetEntry(entry);
		flag = type = 0;
		bool valid = false;
		if (bn >= 2 && fn == 1) {
			opposite_energy = NormEnergy(parameters, 0, fs[0], rfe[0]);
			opposite_strip = fs[0];
			side = 1;
			if (bn == 2 && AreStrips(bs, 31, 33)) {
				strip[0] = 31;
				strip[1] = 33;
				GetEnergy(bs, rbe, energy[0], energy[1]);
				valid = true;
				ClassifyEvent(
					Eigen::Vector3d(energy[0], energy[1], opposite_energy),
					extra, "bs32a",
					400.0, 500.0, 600.0,
					line_distance, plane_distance[0], plane_distance[1],
					type, flag
				);
				++num_bs32a[0];
				if (type) ++num_bs32a[type];
			} else if (bs[0] == 31 && bs[1] == 34) {
				strip[0] = 31;
				strip[1] = 33;
				energy[0] = rbe[0];
				energy[1] = rbe[1];
				valid = true;
				ClassifyEvent(
					Eigen::Vector3d(energy[0], energy[1], opposite_energy),
					extra, "bs32b",
					500.0, 500.0, 500.0,
					line_distance, plane_distance[0], plane_distance[1],
					type, flag
				);
				++num_bs32b[0];
				if (type) ++num_bs32b[type];
			} else if (bs[0] == 34 && bs[1] == 31) {
				strip[0] = 34;
				strip[1] = 31;
				energy[0] = rbe[0];
				energy[1] = rbe[1];
				valid = true;
				ClassifyEvent(
					Eigen::Vector3d(energy[0], energy[1], opposite_energy),
					extra, "bs32c",
					500.0, 500.0, 600.0,
					line_distance, plane_distance[0], plane_distance[1],
					type, flag
				);
				++num_bs32c[0];
				if (type) ++num_bs32c[type];
			} else if (AreStrips(bs, 35, 37)) {
				strip[0] = 35;
				strip[1] = 37;
				GetEnergy(bs, rbe, energy[0], energy[1]);
				valid = true;
				ClassifyEvent(
					Eigen::Vector3d(energy[0], energy[1], opposite_energy),
					extra, "bs36",
					600.0, 600.0, 600.0,
					line_distance, plane_distance[0], plane_distance[1],
					type, flag
				);
				++num_bs36[0];
				if (type) ++num_bs36[type];
			} else if (AreStrips(bs, 47, 49)) {
				strip[0] = 47;
				strip[1] = 49;
				GetEnergy(bs, rbe, energy[0], energy[1]);
				valid = true;
				ClassifyEvent(
					Eigen::Vector3d(energy[0], energy[1], opposite_energy),
					extra, "bs48",
					600.0, 700.0, 700.0,
					line_distance, plane_distance[0], plane_distance[1],
					type, flag
				);
				++num_bs48[0];
				if (type) ++num_bs48[type];
			} else if (AreStrips(bs, 67, 69)) {
				strip[0] = 67;
				strip[1] = 69;
				GetEnergy(bs, rbe, energy[0], energy[1]);
				valid = true;
				ClassifyEvent(
					Eigen::Vector3d(energy[0], energy[1], opposite_energy),
					extra, "bs68",
					500.0, 600.0, 600.0,
					line_distance, plane_distance[0], plane_distance[1],
					type, flag
				);
				++num_bs68[0];
				if (type) ++num_bs68[type];
				gbs68[0].AddPoint(energy[0], energy[1], opposite_energy);
				if (type) gbs68[type].AddPoint(energy[0], energy[1], opposite_energy);
			}
		} else if (fn >= 2 && bn == 1) {
			opposite_energy = brill::NormEnergy(parameters, 1, bs[0], rbe[0]);
			opposite_strip = bs[0];
			side = 0;
			if (AreStrips(fs, 69, 71)) {
				strip[0] = 69;
				strip[1] = 71;
				GetEnergy(fs, rfe, energy[0], energy[1]);
				valid = true;
				ClassifyEvent(
					Eigen::Vector3d(energy[0], energy[1], opposite_energy),
					extra, "fs70",
					300.0, 300.0, 300.0,
					line_distance, plane_distance[0], plane_distance[1],
					type, flag
				);
				++num_fs70[0];
				if (type) ++num_fs70[type];
			} else if (AreStrips(fs, 17, 20)) {
				strip[0] = 17;
				strip[1] = 20;
				GetEnergy(fs, rfe, energy[0], energy[1]);
				valid = true;
				ClassifyEvent(
					Eigen::Vector3d(energy[0], energy[1], opposite_energy),
					extra, "fs1819",
					300.0, 300.0, 200.0,
					line_distance, plane_distance[0], plane_distance[1],
					type, flag
				);
				++num_fs1819[0];
				if (type) ++num_fs1819[type];
			} else if (AreStrips(fs, 1, 3)) {
				strip[0] = 1;
				strip[1] = 3;
				GetEnergy(fs, rfe, energy[0], energy[1]);
				valid = true;
				ClassifyEvent(
					Eigen::Vector3d(energy[0], energy[1], opposite_energy),
					extra, "fs2",
					200.0, 600.0, 600.0,
					line_distance, plane_distance[0], plane_distance[1],
					type, flag
				);
				++num_fs2[0];
				if (type) ++num_fs2[type];
			}
		}
		if (valid) opt.Fill();
	}
	printf("\b\b\b\b100%%\n");
	gbs68[0].Write("gcbs68a");
	gbs68[1].Write("gcbs68b");
	gbs68[2].Write("gcbs68c");
	gbs68[3].Write("gcbs68d");
	opt.Write();

	ReportClassificationResult(num_bs32a, "front strip 31 and 33");
	ReportClassificationResult(num_bs32b, "front strip 31 and 34");
	ReportClassificationResult(num_bs32c, "front strip 34 and 31");
	ReportClassificationResult(num_bs36, "front strip 35 and 37");
	ReportClassificationResult(num_bs48, "front strip 47 and 49");
	ReportClassificationResult(num_bs68, "front strip 67 and 69");
	ReportClassificationResult(num_fs70, "back strip 69 and 71");
	ReportClassificationResult(num_fs1819, "back strip 17 and 20");
	ReportClassificationResult(num_fs2, "back strip 1 and 3");
	std::cout << "---------------------------------------------------------\n";
}

int main(int argc, char **argv) {
	cxxopts::Options options(
		"normalize_extra_t0d2",
		"Normalize special strips in T0D2 only in RCNP2025 experiment."
	);
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Run number.", cxxopts::value<int>(), "run")
		("e,end-run", "End run number.", cxxopts::value<int>(), "run")
		("t,trigger", "Trigger type.", cxxopts::value<std::string>(), "trigger")
		(
			"c,config",
			"Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"),
			"file"
		);

	auto result = options.parse(argc, argv);
	if (result.count("help")) {
		std::cout << options.help() << std::endl;
		return 0;
	}
	if (!result.count("run")) {
		std::cerr << "Error: Missing required option --run.\n";
		std::cout << options.help() << std::endl;
		return 1;
	}

	brill::AppConfig config;
	if (brill::LoadConfig(result["config"].as<std::string>(), config)) {
		return 1;
	}
	if (result.count("trigger")) {
		config.trigger = result["trigger"].as<std::string>();
	}
	const int run = result["run"].as<int>();
	const int end_run = result.count("end-run") ? result["end-run"].as<int>() : run;
	if (end_run < run) {
		std::cerr << "Error: end run " << end_run << " is smaller than run " << run << ".\n";
		return -1;
	}

	const brill::SiliconDetectorConfig *detector =
		brill::FindDetectorConfig(config, "t0d2");
	if (!detector) {
		std::cerr << "Error: Detector t0d2 is not found in config.\n";
		return 1;
	}

	TChain chain("tree");
	int added_runs = 0;
	for (int current_run = run; current_run <= end_run; ++current_run) {
		if (brill::IsJumpRun(config, current_run)) continue;
		++added_runs;
		chain.Add(TString::Format(
			"%s/t0d2_%s%04d.root",
			brill::JoinPath(config.workspace, config.paths.ingot).c_str(),
			brill::TriggerInfix(config.trigger).c_str(),
			current_run
		));
	}
	if (added_runs == 0) {
		std::cout << "No runs to process after jumping runs.\n";
		return 0;
	}
	brill::DssdEvent raw_event;
	brill::SetupInput(&chain, raw_event);

	brill::DssdNormalizeParameters parameters;
	parameters.front_strips = detector->front_strips;
	parameters.back_strips = detector->back_strips;

	std::string normalize_dir = brill::JoinPath(config.workspace, config.paths.normalize);
	TString front_path = TString::Format(
		"%s/t0d2_front_%04d.txt",
		normalize_dir.c_str(),
		run
	);
	TString back_path = TString::Format(
		"%s/t0d2_back_%04d.txt",
		normalize_dir.c_str(),
		run
	);
	if (brill::ReadDssdNormalizeParameters(
		front_path.Data(), back_path.Data(), parameters
	)) {
		std::cerr << "Error: Read normalize parameters failed.\n";
		return 1;
	}

	TString output_path = TString::Format(
		"%s/t0d2_%sextra_%04d_%04d.root",
		normalize_dir.c_str(),
		brill::TriggerInfix(config.trigger).c_str(),
		run,
		end_run
	);
	TFile opf(output_path, "recreate");
	brill::T0D2ExtraNormalizeParameters extra_parameters;
	LineNormalize(chain, raw_event, parameters, extra_parameters);
	LineEstimate(chain, raw_event, parameters, extra_parameters);
	PlaneNormalize(chain, raw_event, parameters, extra_parameters);
	PlaneEstimate(chain, raw_event, parameters, extra_parameters);
	Classify(chain, raw_event, parameters, extra_parameters);
	opf.Close();

	// write to file
	TString t0d2_extra_path = TString::Format(
		"%s/t0d2_extra_%04d.txt",
		normalize_dir.c_str(),
		run
	);
	if (brill::WriteExtraNormalizeParameters(
		t0d2_extra_path.Data(),
		extra_parameters
	)) {
		std::cerr << "Error: write t0d2_extra parameter file "
			<< t0d2_extra_path.Data() << " failed." << std::endl;
	}

	return 0;
}