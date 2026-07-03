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

void FitNonLinear(
	TGraph *gfs101,
	brill::T0D1ExtraNormalizeParameters &extra_parameters
) {
	// piecewise fit fs101
	TF1 ffs101a("ffs101a", "pol1", 0, 6000);
	gfs101[0].Fit(&ffs101a, "R+ ROB=0.9");
	TF1 ffs101b("ffs101b", "pol1", 4600, 5200);
	gfs101[1].Fit(&ffs101b, "R+ ROB=0.9");
	TF1 ffs101c0("ffs101c0", "pol2", 0, 5900);
	gfs101[2].Fit(&ffs101c0, "R+ ROB=0.9");
	TF1 ffs101c1("ffs101c1", "pol2", 5900, 50000);
	gfs101[2].Fit(&ffs101c1, "R+");

	const double range[4][4] = {
		-1.0, -1.0, 0.0, 3900.0,
		-1.0, -1.0, 0.0, 5900.0,
		0.0, 5900.0, 5900.0, -1.0,
		5900.0, 50000.0, 5900.0, -1.0
	};
	for (int i = 0; i < 4; ++i) {
		extra_parameters.non_linear.xmin[i] = range[i][0];
		extra_parameters.non_linear.xmax[i] = range[i][1];
		extra_parameters.non_linear.ymin[i] = range[i][2];
		extra_parameters.non_linear.ymax[i] = range[i][3];
	}
	extra_parameters.non_linear.p0[0] = ffs101a.GetParameter(0);
	extra_parameters.non_linear.p1[0] = ffs101a.GetParameter(1);
	extra_parameters.non_linear.p2[0] = 0.0;
	extra_parameters.non_linear.p0[1] = ffs101b.GetParameter(0);
	extra_parameters.non_linear.p1[1] = ffs101b.GetParameter(1);
	extra_parameters.non_linear.p2[1] = 0.0;
	extra_parameters.non_linear.p0[2] = ffs101c0.GetParameter(0);
	extra_parameters.non_linear.p1[2] = ffs101c0.GetParameter(1);
	extra_parameters.non_linear.p2[2] = ffs101c0.GetParameter(2);
	extra_parameters.non_linear.p0[3] = ffs101c1.GetParameter(0);
	extra_parameters.non_linear.p1[3] = ffs101c1.GetParameter(1);
	extra_parameters.non_linear.p2[3] = ffs101c1.GetParameter(2);
}

void FillStripLine(
	const int *strip,
	const double *energy,
	const double other_energy,
	const double min,
	const double max,
	TGraph2D &gn,
	TH1F &hd,
	TGraph2D &g,
	std::vector<Eigen::Vector3d> &points
) {
	double e0, e1;
	GetEnergy(strip, energy, e0, e1);
	gn.AddPoint(e0, e1, other_energy);
	hd.Fill(e0 - e1);
	if (e0-e1 > min && e0-e1 < max) {
		g.AddPoint(e0, e1, other_energy);
		Eigen::Vector3d point(e0, e1, other_energy);
		points.push_back(point);
	}
}

void LineNormalize(
	TChain &chain,
	brill::DssdEvent &raw_event,
	brill::DssdNormalizeParameters &parameters,
	brill::T0D1ExtraNormalizeParameters &extra_parameters
) {
	// for convenient
	const int &fn = raw_event.front_num;
	const int &bn = raw_event.back_num;
	const int *fs = raw_event.front_strip;
	const int *bs = raw_event.back_strip;
	const double *rfe = raw_event.front_integral;
	const double *rbe = raw_event.back_integral;

	// energy difference of ajacent strips
	TH1F hfs50("hfs50", "energy difference of front strip 49 and 51", 100, -2000, 2000);
	TH1F hfs94("hfs94", "energy difference of front strip 93 and 95", 100, -2000, 2000);
	TH1F hfs120("hfs120", "energy difference of front strip 119 and 121", 100, -2000, 2000);
	TH1F hbs104("hbs104", "energy difference of back strip 104 and 109", 100, -2000, 2000);
	// non-linear strip
	TGraph gfs101[3];
	// no energy cut
	TGraph2D gfs50n, gfs94n, gfs120n, gbs104n;
	// cut by threshold get from energy difference histogram
	TGraph2D gfs50, gfs94, gfs120, gbs104;
	// points for PCA
	std::vector<Eigen::Vector3d> pfs50, pfs94, pfs120, pbs104;

	long long total = chain.GetEntries();
	long long last_percentage = 0;
	printf("Normalizing extra t0d1   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		if (entry * 100 / total > last_percentage) {
			last_percentage = entry * 100 / total;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}
		chain.GetEntry(entry);
		if (fn == 1 && bn == 1 && fs[0] == 101) {
			double be = NormEnergy(parameters, 1, bs[0], rbe[0]);
			if (be < 3900) gfs101[0].AddPoint(rfe[0], be);
			else if (be < 5900) gfs101[1].AddPoint(rfe[0], be);
			else gfs101[2].AddPoint(rfe[0], be);
		} else if (fn == 1 && bn == 2 && fs[0] != 101 && AreStrips(bs, 104, 109)) {
			double be104, be109;
			GetEnergy(bs, rbe, be104, be109);
			double fe = NormEnergy(parameters, 0, fs[0], rfe[0]);
			hbs104.Fill(be104 - be109);
			gbs104n.AddPoint(be104, be109, fe);
			if (fabs(be104-be109) < 500.0) {
				gbs104.AddPoint(be104, be109, fe);
				Eigen::Vector3d point(be104, be109, fe);
				pbs104.push_back(point);
			}
		} else if (fn == 2 && bn == 1 && bs[0] != 104 && bs[0] != 109) {
			double be = NormEnergy(parameters, 1, bs[0], rbe[0]);
			if (AreStrips(fs, 49, 51)) {
				double fe49, fe51;
				GetEnergy(fs, rfe, fe49, fe51);
				gfs50n.AddPoint(fe49, fe51, be);
				hfs50.Fill(fe49 - fe51);
				if (fe49-fe51 > -600.0 && fe49-fe51 < 0.0) {
					gfs50.AddPoint(fe49, fe51, be);
					Eigen::Vector3d point(fe49, fe51, be);
					pfs50.push_back(point);
				}
			} else if (AreStrips(fs, 93, 95)) {
				double fe93, fe95;
				GetEnergy(fs, rfe, fe93, fe95);
				hfs94.Fill(fe93 - fe95);
				gfs94n.AddPoint(fe93, fe95, be);
				if (fe93-fe95 > 0.0 && fe93-fe95 < 1000.0) {
					gfs94.AddPoint(fe93, fe95, be);
					Eigen::Vector3d point(fe93, fe95, be);
					pfs94.push_back(point);
				}
			} else if (AreStrips(fs, 119, 121)) {
				double fe119, fe121;
				GetEnergy(fs, rfe, fe119, fe121);
				hfs120.Fill(fe119 - fe121);
				if (fe119-fe121 > -1000.0 && fe119-fe121 < 0.0) {
					gfs120.AddPoint(fe119, fe121, be);
					Eigen::Vector3d point(fe119, fe121, be);
					pfs120.push_back(point);
				}
				gfs120n.AddPoint(fe119, fe121, be);
			}
		}
	}
	printf("\b\b\b\b100%%\n");

	// fit non-linear fs101
	FitNonLinear(gfs101, extra_parameters);
	std::cout << "---------------------------------------------------------\n";

	// PCA
	const std::vector<Eigen::Vector3d>* pca_points[] = {
		&pfs50, &pfs94, &pfs120, &pbs104
	};
	const std::string pca_desriptions[] = {
		"Front strip 50",
		"Front strip 94",
		"Front strip 120",
		"Back strip 104 and 109"
	};
	const std::string pca_names[] = {
		"fs50",
		"fs94",
		"fs120",
		"bs104"
	};
	for (int i = 0; i < 4; ++i) {
		PCAPrint(
			*pca_points[i],
			pca_desriptions[i],
			pca_names[i],
			true,
			extra_parameters
		);
	}
	std::cout << "---------------------------------------------------------\n";

	hfs50.Write();
	hfs94.Write();
	hfs120.Write();
	hbs104.Write();
	gfs101[0].Write("gfs101a");
	gfs101[1].Write("gfs101b");
	gfs101[2].Write("gfs101c");
	gfs50.Write("gfs50");
	gfs94.Write("gfs94");
	gfs120.Write("gfs120");
	gbs104n.Write("gbs104n");
	gfs50n.Write("gfs50n");
	gfs94n.Write("gfs94n");
	gfs120n.Write("gfs120n");
	gbs104.Write("gbs104");
}


void LineEstimate(
	TChain &chain,
	brill::DssdEvent &raw_event,
	brill::DssdNormalizeParameters &parameters,
	brill::T0D1ExtraNormalizeParameters &extra
) {
	// for convenient
	const int &fn = raw_event.front_num;
	const int &bn = raw_event.back_num;
	const int *fs = raw_event.front_strip;
	const int *bs = raw_event.back_strip;
	const double *rfe = raw_event.front_integral;
	const double *rbe = raw_event.back_integral;

	TGraph gfs101;
	TH1F hfs101("hdfs101", "Front strip 101 difference", 1000, -5000, 5000);
	TH1F hfs50("hdfs50", "Front strip 50 distance", 200, 0, 2000);
	TH1F hfs94("hdfs94", "Front strip 94 distance", 200, 0, 2000);
	TH1F hfs120("hdfs120", "Front strip 120 distance", 200, 0, 2000);
	TH1F hbs104("hdbs104", "Back strip 104 distance", 200, 0, 2000);

	long long total = chain.GetEntries();
	long long last_percentage = 0;
	printf("Estimating extra t0d1   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		if (entry * 100 / total > last_percentage) {
			last_percentage = entry * 100 / total;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}
		chain.GetEntry(entry);
		if (fn == 1 && bn == 1 && fs[0] == 101) {
			double be = NormEnergy(parameters, 1, bs[0], rbe[0]);
			double fe = 0.0;
			if (be < 3900) {
				fe = extra.non_linear.p0[0] + extra.non_linear.p1[0] * rfe[0];
			} else if (be < 5900) {
				fe = extra.non_linear.p0[1] + extra.non_linear.p1[1] * rfe[0];
			} else if (raw_event.front_integral[0] < 5900) {
				fe = extra.non_linear.p0[2]
					+ extra.non_linear.p1[2] * rfe[0]
					+ extra.non_linear.p2[2] * rfe[0] * rfe[0];
			} else {
				fe = extra.non_linear.p0[3]
					+ extra.non_linear.p1[3] * rfe[0]
					+ extra.non_linear.p2[3] * rfe[0] * rfe[0];
			}
			gfs101.AddPoint(fe, be);
			hfs101.Fill(fe-be);
		} else if (fn == 1 && bn == 2 && fs[0] != 101 && AreStrips(bs, 104, 109)) {
			double be104, be109;
			GetEnergy(bs, rbe, be104, be109);
			double fe = NormEnergy(parameters, 0, fs[0], rfe[0]);
			Eigen::Vector3d point(be104, be109, fe);
			hbs104.Fill(brill::PointLineDistance(
				point,
				extra.pca["bs104"].mean,
				extra.pca["bs104"].direction
			));
		} else if (fn == 2 && bn == 1 && bs[0] != 104 && bs[0] != 109) {
			double be = NormEnergy(parameters, 1, bs[0], rbe[0]);
			if (AreStrips(fs, 49, 51)) {
				double fe49, fe51;
				GetEnergy(fs, rfe, fe49, fe51);
				Eigen::Vector3d point(fe49, fe51, be);
				hfs50.Fill(brill::PointLineDistance(
					point,
					extra.pca["fs50"].mean,
					extra.pca["fs50"].direction
				));
			} else if (AreStrips(fs, 93, 95)) {
				double fe93, fe95;
				GetEnergy(fs, rfe, fe93, fe95);
				Eigen::Vector3d point(fe93, fe95, be);
				hfs94.Fill(brill::PointLineDistance(
					point,
					extra.pca["fs94"].mean,
					extra.pca["fs94"].direction
				));
			} else if (AreStrips(fs, 119, 121)) {
				double fe119, fe121;
				GetEnergy(fs, rfe, fe119, fe121);
				Eigen::Vector3d point(fe119, fe121, be);
				hfs120.Fill(brill::PointLineDistance(
					point,
					extra.pca["fs120"].mean,
					extra.pca["fs120"].direction
				));
			}
		}
	}
	printf("\b\b\b\b100%%\n");

	gfs101.Write("gfs101n");
	hfs101.Write();
	hfs50.Write();
	hfs94.Write();
	hfs120.Write();
	hbs104.Write();
}


void PlaneNormalize(
	TChain &chain,
	brill::DssdEvent &raw_event,
	brill::DssdNormalizeParameters &parameters,
	brill::T0D1ExtraNormalizeParameters &extra
) {
	// for convenient
	const int &fn = raw_event.front_num;
	const int &bn = raw_event.back_num;
	const int *fs = raw_event.front_strip;
	const int *bs = raw_event.back_strip;
	const double *rfe = raw_event.front_integral;
	const double *rbe = raw_event.back_integral;

	TH1F hafs50("hafs50", "angle to line of front strip 49 and 51", 100, -1.0, 1.0);
	TH1F hafs94("hafs94", "angle to line of front strip 93 and 95", 100, -1.0, 1.0);
	TH1F hafs120("hafs120", "angle to line of front strip 119 and 121", 100, -1.0, 1.0);
	std::vector<Eigen::Vector3d> pfs50[2], pfs94[2], pfs120[2];

	long long total = chain.GetEntries();
	long long last_percentage = 0;
	printf("Normalizing plane in t0d1   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		if (entry * 100 / total > last_percentage) {
			last_percentage = entry * 100 / total;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}
		chain.GetEntry(entry);
		if (fn == 2 && bn == 1) {
			double be = NormEnergy(parameters, 1, bs[0], rbe[0]);
			if (AreStrips(fs, 49, 51)) {
				double fe49, fe51;
				GetEnergy(fs, rfe, fe49, fe51);
				Eigen::Vector3d point(fe49, fe51, be);
				const Eigen::Vector3d &centroid = extra.pca["fs50"].mean;
				const Eigen::Vector3d &direction = extra.pca["fs50"].direction;
				double distance = brill::PointLineDistance(point, centroid, direction);
				if (distance < 400.0) continue;
				double angle = brill::PointVerticalAngle(point, centroid, direction);
				hafs50.Fill(angle);
				if (angle > -0.8 && angle < -0.4) pfs50[0].push_back(point);
				else if (angle > 0.8) pfs50[1].push_back(point);
			} else if (AreStrips(fs, 93, 95)) {
				double fe93, fe95;
				GetEnergy(fs, rfe, fe93, fe95);
				Eigen::Vector3d point(fe93, fe95, be);
				const Eigen::Vector3d &centroid = extra.pca["fs94"].mean;
				const Eigen::Vector3d &direction = extra.pca["fs94"].direction;
				double distance = brill::PointLineDistance(point, centroid, direction);
				if (distance < 600.0) continue;
				double angle = brill::PointVerticalAngle(point, centroid, direction);
				hafs94.Fill(angle);
				if (angle > -0.4 && angle < 0.0) pfs94[0].push_back(point);
				else if (angle > 0.6) pfs94[1].push_back(point);
			} else if (AreStrips(fs, 119, 121)) {
				double fe119, fe121;
				GetEnergy(fs, rfe, fe119, fe121);
				Eigen::Vector3d point(fe119, fe121, be);
				const Eigen::Vector3d &centroid = extra.pca["fs120"].mean;
				const Eigen::Vector3d &direction = extra.pca["fs120"].direction;
				double distance = brill::PointLineDistance(point, centroid, direction);
				if (distance < 400.0) continue;
				double angle = brill::PointVerticalAngle(point, centroid, direction);
				hafs120.Fill(angle);
				if (angle < -0.8) pfs120[0].push_back(point);
				else if (angle > 0.5 && angle < 0.9) pfs120[1].push_back(point);
			}
		}
	}
	printf("\b\b\b\b100%%\n");

	// PCA for plane
	const std::vector<Eigen::Vector3d>* pca_points[6] = {
		pfs50, pfs50+1, pfs94, pfs94+1, pfs120, pfs120+1
	};
	const std::string pca_descriptions[6] = {
		"Front strip 50 plane 0",
		"Front strip 50 plane 1",
		"Front strip 94 plane 0",
		"Front strip 94 plane 1",
		"Front strip 120 plane 0",
		"Front strip 120 plane 1"
	};
	const std::string pca_names[6] = {
		"fs50p0", "fs50p1", "fs94p0", "fs94p1", "fs120p0", "fs120p1"
	};
	for (int i = 0; i < 6; ++i) {
		PCAPrint(
			*pca_points[i],
			pca_descriptions[i],
			pca_names[i],
			false,
			extra
		);
	}
	std::cout << "---------------------------------------------------------\n";

	hafs50.Write();
	hafs94.Write();
	hafs120.Write();
}


void PlaneEstimate(
	TChain &chain,
	brill::DssdEvent &raw_event,
	brill::DssdNormalizeParameters &parameters,
	brill::T0D1ExtraNormalizeParameters &extra
) {
	// for convenient
	const int &fn = raw_event.front_num;
	const int &bn = raw_event.back_num;
	const int *fs = raw_event.front_strip;
	const int *bs = raw_event.back_strip;
	const double *rfe = raw_event.front_integral;
	const double *rbe = raw_event.back_integral;

	TH1F hdfs50_p0("hdp0fs50", "distance to plane 0 for front strip 49 and 51", 100, 0, 5000);
	TH1F hdfs50_p1("hdp1fs50", "distance to plane 1 for front strip 49 and 51", 100, 0, 5000);
	TH1F hdfs94_p0("hdp0fs94", "distance to plane 0 for front strip 93 and 95", 100, 0, 5000);
	TH1F hdfs94_p1("hdp1fs94", "distance to plane 1 for front strip 93 and 95", 100, 0, 5000);
	TH1F hdfs120_p0("hdp0fs120", "distance to plane 0 for front strip 119 and 121", 100, 0, 5000);
	TH1F hdfs120_p1("hdp1fs120", "distance to plane 1 for front strip 119 and 121", 100, 0, 5000);
	TH1F hdfs50_p0g("hdp0fs50g", "distance to plane 0 for front strip 49 and 51 gated", 100, 0, 5000);
	TH1F hdfs50_p1g("hdp1fs50g", "distance to plane 1 for front strip 49 and 51 gated", 100, 0, 5000);
	std::vector<Eigen::Vector3d> pfs50[2], pfs94[2], pfs120[2];

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
		if (fn == 2 && bn == 1) {
			double be = NormEnergy(parameters, 1, bs[0], rbe[0]);
			if (AreStrips(fs, 49, 51)) {
				double fe49, fe51;
				GetEnergy(fs, rfe, fe49, fe51);
				Eigen::Vector3d point(fe49, fe51, be);
				double distance_p0 = brill::PointPlaneDistance(
					point, extra.pca["fs50p0"].mean, extra.pca["fs50p0"].direction
				);
				double distance_p1 = brill::PointPlaneDistance(
					point, extra.pca["fs50p1"].mean, extra.pca["fs50p1"].direction
				);
				hdfs50_p0.Fill(distance_p0);
				hdfs50_p1.Fill(distance_p1);
				const Eigen::Vector3d &centroid = extra.pca["fs50"].mean;
				const Eigen::Vector3d &direction = extra.pca["fs50"].direction;
				double distance = brill::PointLineDistance(point, centroid, direction);
				if (distance > 400.0) {
					double angle = brill::PointVerticalAngle(point, centroid, direction);
					if (angle > -0.8 && angle < -0.4) hdfs50_p0g.Fill(distance_p0);
					else if (angle > 0.8) hdfs50_p1g.Fill(distance_p1);
				}
			} else if (AreStrips(fs, 93, 95)) {
				double fe93, fe95;
				GetEnergy(fs, rfe, fe93, fe95);
				Eigen::Vector3d point(fe93, fe95, be);
				double distance_p0 = brill::PointPlaneDistance(
					point, extra.pca["fs94p0"].mean, extra.pca["fs94p0"].direction
				);
				double distance_p1 = brill::PointPlaneDistance(
					point, extra.pca["fs94p1"].mean, extra.pca["fs94p1"].direction
				);
				hdfs94_p0.Fill(distance_p0);
				hdfs94_p1.Fill(distance_p1);
			} else if (AreStrips(fs, 119, 121)) {
				double fe119, fe121;
				GetEnergy(fs, rfe, fe119, fe121);
				Eigen::Vector3d point(fe119, fe121, be);
				double distance_p0 = brill::PointPlaneDistance(
					point, extra.pca["fs120p0"].mean, extra.pca["fs120p0"].direction
				);
				double distance_p1 = brill::PointPlaneDistance(
					point, extra.pca["fs120p1"].mean, extra.pca["fs120p1"].direction
				);
				hdfs120_p0.Fill(distance_p0);
				hdfs120_p1.Fill(distance_p1);
			}
		}
	}
	printf("\b\b\b\b100%%\n");

	hdfs50_p0.Write();
	hdfs50_p1.Write();
	hdfs94_p0.Write();
	hdfs94_p1.Write();
	hdfs120_p0.Write();
	hdfs120_p1.Write();
	hdfs50_p0g.Write();
	hdfs50_p1g.Write();
}


void Classify(
	TChain &chain,
	brill::DssdEvent &raw_event,
	brill::DssdNormalizeParameters &parameters,
	brill::T0D1ExtraNormalizeParameters &extra
) {
	// for convenient
	const int &fn = raw_event.front_num;
	const int &bn = raw_event.back_num;
	const int *fs = raw_event.front_strip;
	const int *bs = raw_event.back_strip;
	const double *rfe = raw_event.front_integral;
	const double *rbe = raw_event.back_integral;

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

	int num_fs50[4], num_fs94[4], num_fs120[4], num_bs104[2];
	for (int i = 0; i < 4; ++i) {
		num_fs50[i] = 0;
		num_fs94[i] = 0;
		num_fs120[i] = 0;
	}
	for (int i = 0; i < 2; ++i) num_bs104[i] = 0;

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
		flag = 0;
		type = 0;
		bool valid = false;
		if (fn == 2 && bn == 1) {
			side = 0;
			opposite_energy = NormEnergy(parameters, 1, bs[0], rbe[0]);
			opposite_strip = bs[0];
			if (AreStrips(fs, 49, 51)) {
				strip[0] = 49;
				strip[1] = 51;
				GetEnergy(fs, rfe, energy[0], energy[1]);

				valid = true;
				Eigen::Vector3d point(energy[0], energy[1], opposite_energy);
				line_distance = brill::PointLineDistance(
					point, extra.pca["fs50"].mean, extra.pca["fs50"].direction
				);
				if (line_distance < 400.0) {
					type = 1;
					flag |= 1;
				} else {
					plane_distance[0] = brill::PointPlaneDistance(
						point, extra.pca["fs50p0"].mean, extra.pca["fs50p0"].direction
					);
					plane_distance[1] = brill::PointPlaneDistance(
						point, extra.pca["fs50p1"].mean, extra.pca["fs50p1"].direction
					);
					if (plane_distance[0] < 400.0 && plane_distance[1] < 400.0) {
						type = plane_distance[0] < plane_distance[1] ? 2 : 3;
						flag |= 6;
					} else if (plane_distance[0] < 400.0) {
						type = 2;
						flag |= 2;
					} else if (plane_distance[1] < 400.0) {
						type = 3;
						flag |= 4;
					}
				}
				++num_fs50[0];
				if (type) ++num_fs50[type];
			} else if (AreStrips(fs, 93, 95)) {
				strip[0] = 93;
				strip[1] = 95;
				GetEnergy(fs, rfe, energy[0], energy[1]);

				valid = true;
				Eigen::Vector3d point(energy[0], energy[1], opposite_energy);
				line_distance = brill::PointLineDistance(
					point, extra.pca["fs94"].mean, extra.pca["fs94"].direction
				);
				if (line_distance < 600.0) {
					type = 1;
					flag |= 1;
					++num_fs94[1];
				Eigen::Vector3d point(energy[0], energy[1], opposite_energy);
				} else {
					plane_distance[0] = brill::PointPlaneDistance(
						point, extra.pca["fs94p0"].mean, extra.pca["fs94p0"].direction
					);
					plane_distance[1] = brill::PointPlaneDistance(
						point, extra.pca["fs94p1"].mean, extra.pca["fs94p1"].direction
					);
					if (plane_distance[0] < 400.0 && plane_distance[1] < 400.0) {
						type = plane_distance[0] < plane_distance[1] ? 2 : 3;
						flag |= 6;
					} else if (plane_distance[0] < 400.0) {
						type = 2;
						flag |= 2;
					} else if (plane_distance[1] < 400.0) {
						type = 3;
						flag |= 4;
					}
				}
				++num_fs94[0];
				if (type) ++num_fs94[type];
			} else if (AreStrips(fs, 119, 121)) {
				strip[0] = 119;
				strip[1] = 121;
				GetEnergy(fs, rfe, energy[0], energy[1]);

				valid = true;
				Eigen::Vector3d point(energy[0], energy[1], opposite_energy);
				line_distance = brill::PointLineDistance(
					point, extra.pca["fs120"].mean, extra.pca["fs120"].direction
				);
				if (line_distance < 400.0) {
					type = 1;
					flag |= 1;
					++num_fs120[1];
				} else {
					plane_distance[0] = brill::PointPlaneDistance(
						point, extra.pca["fs120p0"].mean, extra.pca["fs120p0"].direction
					);
					plane_distance[1] = brill::PointPlaneDistance(
						point, extra.pca["fs120p1"].mean, extra.pca["fs120p1"].direction
					);
					if (plane_distance[0] < 400.0 && plane_distance[1] < 400.0) {
						type = plane_distance[0] < plane_distance[1] ? 2 : 3;
						flag |= 6;
						++num_fs120[2];
					} else if (plane_distance[0] < 500.0) {
						type = 2;
					} else if (plane_distance[1] < 500.0) {
						type = 3;
					}
				}
				++num_fs120[0];
				if (type) ++num_fs120[type];
			}
		} else if (bn == 2 && fn == 1 && AreStrips(bs, 104, 109)) {
			side = 1;
			strip[0] = 104;
			strip[1] = 109;
			opposite_strip = bs[0];
			opposite_energy = NormEnergy(parameters, 0, fs[0], rfe[0]);
			GetEnergy(bs, rbe, energy[0], energy[1]);

			++num_bs104[0];
			valid = true;
			Eigen::Vector3d point(energy[0], energy[1], opposite_energy);
			line_distance = brill::PointLineDistance(
				point, extra.pca["bs104"].mean, extra.pca["bs104"].direction
			);
			if (line_distance < 800.0) {
				type = 1;
				flag |= 1;
				++num_bs104[1];
			}
			plane_distance[0] = plane_distance[1] = 0.0;
		}
		if (valid) opt.Fill();
	}
	printf("\b\b\b\b100%%\n");

	// report
	std::cout
		<< "Classified front strip 50: \n"
		<< "  line   : " << num_fs50[1] << " / " << num_fs50[0] << " = "
		<< std::fixed << std::setprecision(2) << num_fs50[1] * 100.0 / num_fs50[0] << "%\n"
		<< "  plane 0: " << num_fs50[2] << " / " << num_fs50[0] << " = "
		<< std::fixed << std::setprecision(2) << num_fs50[2] * 100.0 / num_fs50[0] << "%\n"
		<< "  plane 1: " << num_fs50[3] << " / " << num_fs50[0] << " = "
		<< std::fixed << std::setprecision(2) << num_fs50[3] * 100.0 / num_fs50[0] << "%\n"
		<< "  total  : " << num_fs50[1]+num_fs50[2]+num_fs50[3] << " / " << num_fs50[0] << " = "
		<< std::fixed << std::setprecision(2)
		<< (num_fs50[1]+num_fs50[2]+num_fs50[3]) * 100.0 / num_fs50[0] << "%\n"
		<< "Classified front strip 94: \n"
		<< "  line   : " << num_fs94[1] << " / " << num_fs94[0] << " = "
		<< std::fixed << std::setprecision(2) << num_fs94[1] * 100.0 / num_fs94[0] << "%\n"
		<< "  plane 0: " << num_fs94[2] << " / " << num_fs94[0] << " = "
		<< std::fixed << std::setprecision(2) << num_fs94[2] * 100.0 / num_fs94[0] << "%\n"
		<< "  plane 1: " << num_fs94[3] << " / " << num_fs94[0] << " = "
		<< std::fixed << std::setprecision(2) << num_fs94[3] * 100.0 / num_fs94[0] << "%\n"
		<< "  total  : " << num_fs94[1]+num_fs94[2]+num_fs94[3] << " / " << num_fs94[0] << " = "
		<< std::fixed << std::setprecision(2)
		<< (num_fs94[1]+num_fs94[2]+num_fs94[3]) * 100.0 / num_fs94[0] << "%\n"
		<< "Classified front strip 120: \n"
		<< "  line   : " << num_fs120[1] << " / " << num_fs120[0] << " = "
		<< std::fixed << std::setprecision(2) << num_fs120[1] * 100.0 / num_fs120[0] << "%\n"
		<< "  plane 0: " << num_fs120[2] << " / " << num_fs120[0] << " = "
		<< std::fixed << std::setprecision(2) << num_fs120[2] * 100.0 / num_fs120[0] << "%\n"
		<< "  plane 1: " << num_fs120[3] << " / " << num_fs120[0] << " = "
		<< std::fixed << std::setprecision(2) << num_fs120[3] * 100.0 / num_fs120[0] << "%\n"
		<< "  total  : " << num_fs120[1]+num_fs120[2]+num_fs120[3] << " / " << num_fs120[0] << " = "
		<< std::fixed << std::setprecision(2)
		<< (num_fs120[1]+num_fs120[2]+num_fs120[3]) * 100.0 / num_fs120[0] << "%\n"
		<< "Classified back strip 104: \n"
		<< "  line   : " << num_bs104[1] << " / " << num_bs104[0] << " = "
		<< std::fixed << std::setprecision(2) << num_bs104[1] * 100.0 / num_bs104[0] << "%\n"
		<< "---------------------------------------------------------\n";

	opt.Write();
}


int main(int argc, char **argv) {
	cxxopts::Options options(
		"normalize_extra_t0d1",
		"Normalize special strips in T0D1 only in RCNP2025 experiment."
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
		brill::FindDetectorConfig(config, "t0d1");
	if (!detector) {
		std::cerr << "Error: Detector t0d1 is not found in config.\n";
		return 1;
	}

	TChain chain("tree");
	int added_runs = 0;
	for (int current_run = run; current_run <= end_run; ++current_run) {
		if (brill::IsJumpRun(config, current_run)) continue;
		++added_runs;
		chain.Add(TString::Format(
			"%s/t0d1_%s%04d.root",
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
		"%s/t0d1_front_%04d.txt",
		normalize_dir.c_str(),
		run
	);
	TString back_path = TString::Format(
		"%s/t0d1_back_%04d.txt",
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
		"%s/t0d1_%sextra_%04d_%04d.root",
		normalize_dir.c_str(),
		brill::TriggerInfix(config.trigger).c_str(),
		run,
		end_run
	);
	TFile opf(output_path, "recreate");
	brill::T0D1ExtraNormalizeParameters extra_parameters;
	LineNormalize(chain, raw_event, parameters, extra_parameters);
	LineEstimate(chain, raw_event, parameters, extra_parameters);
	PlaneNormalize(chain, raw_event, parameters, extra_parameters);
	PlaneEstimate(chain, raw_event, parameters, extra_parameters);
	Classify(chain, raw_event, parameters, extra_parameters);
	opf.Close();

	// write to file
	TString t0d1_extra_path = TString::Format(
		"%s/t0d1_extra_%04d.txt",
		normalize_dir.c_str(),
		run
	);
	if (brill::WriteExtraNormalizeParameters(
		t0d1_extra_path.Data(),
		extra_parameters
	)) {
		std::cerr << "Error: write t0d1_extra parameter file "
			<< t0d1_extra_path.Data() << " failed." << std::endl;
	}

	return 0;
}