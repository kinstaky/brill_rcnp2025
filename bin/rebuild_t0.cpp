#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include "carquet/carquet.h"

#include "external/cxxopts.hpp"
#include "include/config.h"
// #include "include/energy_calculator/lost_energy_calculator.h"
#include "include/event/particle_event.h"
#include "include/event/t0/t0_event.h"
#include "include/t0/dssd.h"
#include "include/utils.h"

constexpr int kCalibrationLayers = 2;
constexpr int kFirstStopLayer = 2;
constexpr int kLastStopLayer = 2;
constexpr int kD2Bit = 0x2;

struct ParticleEntry {
	int charge = 0;
	int mass = 0;
	double energy = 0.0;
	double time = 0.0;
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	bool stop = false;
	int last = -1;
};

// using StopEnergyCache = std::map<std::tuple<int, int, int>, double>;

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

double CalibrateEnergy(
	const brill::CalibrationParameters &parameters,
	int layer,
	double raw_energy
) {
	return parameters.p0[layer] + parameters.p1[layer] * raw_energy;
}

// double RebuildEnergy(
// 	const brill::T0Event &input,
// 	int particle_index,
// 	const brill::CalibrationParameters &parameters,
// 	const double max_stop_energy
// ) {
// 	double energy = 0.0;
// 	int actual_layer = input.layer[particle_index] == 6 ? 5 : input.layer[particle_index];
// 	for (int i = 1; i < input.layer[particle_index]; ++i) {
// 		double layer_energy = CalibrateEnergy(parameters, i, input.energy[particle_index][i]);
// 		if (i == input.layer[particle_index] - 1 && layer_energy > max_stop_energy) {
// 			return 0.0;
// 		}
// 		energy += layer_energy;
// 	}
// 	return energy;
// }

// std::string RangeCachePath(
// 	const brill::AppConfig &config,
// 	int charge,
// 	int mass
// ) {
// 	return TString::Format(
// 		"%s/si_z%d_a%d.root",
// 		brill::JoinPath(config.root.workspace, config.paths.energy_calculator).c_str(),
// 		charge,
// 		mass
// 	).Data();
// }

// double MaximumStopEnergy(
// 	const brill::AppConfig &config,
// 	int charge,
// 	int mass,
// 	int stop_layer,
// 	StopEnergyCache &cache
// ) {
// 	if (stop_layer > 5) return 0.0;
// 	const auto key = std::make_tuple(charge, mass, stop_layer);
// 	auto iterator = cache.find(key);
// 	if (iterator != cache.end()) return iterator->second;

// 	if (
// 		stop_layer < kFirstStopLayer
// 		|| stop_layer > kLastStopLayer
// 		|| size_t(stop_layer - 1) >= config.t0.silicon.size()
// 	) {
// 		return 0.0;
// 	}

// 	const std::string &detector_name = config.t0.silicon[stop_layer - 1];
// 	const auto *detector = config.FindDetector(detector_name);
// 	if (!detector) {
// 		std::cerr
// 			<< "Error: Detector config for "
// 			<< detector_name
// 			<< " not found when rebuilding T0.\n";
// 		return 0.0;
// 	}

// 	brill::LostEnergyCalculator calculator(
// 		charge,
// 		mass,
// 		brill::SiliconMaterial(),
// 		detector->thickness_um,
// 		RangeCachePath(config, charge, mass)
// 	);
// 	double maximum_stop_energy = calculator.IncidentEnergy(0.0);
// 	cache[key] = maximum_stop_energy;
// 	return maximum_stop_energy;
// }


bool CompareParticleEntry(const ParticleEntry &left, const ParticleEntry &right) {
	if (left.charge != right.charge) return left.charge < right.charge;
	if (left.mass != right.mass) return left.mass < right.mass;
	return left.last < right.last;
}

bool CheckCarquetStatus(carquet_status_t status, const char *operation) {
	if (status == CARQUET_OK) return true;
	std::cerr
		<< "Error: " << operation << " failed: "
		<< carquet_status_string(status) << ".\n";
	return false;
}

carquet_schema_t *CreateParquetSchema() {
	carquet_error_t error{};
	carquet_error_init(&error);
	carquet_schema_t *schema = carquet_schema_create(&error);
	if (!schema) {
		std::cerr << "Error: Create parquet schema failed: "
			<< error.message << ".\n";
		return nullptr;
	}

	const std::array<std::pair<const char *, carquet_physical_type_t>, 11> columns = {{
		{"run", CARQUET_PHYSICAL_INT32},
		{"event_id", CARQUET_PHYSICAL_INT32},
		{"multiplicity", CARQUET_PHYSICAL_INT32},
		{"index", CARQUET_PHYSICAL_INT32},
		{"charge", CARQUET_PHYSICAL_INT32},
		{"mass", CARQUET_PHYSICAL_INT32},
		{"energy", CARQUET_PHYSICAL_DOUBLE},
		{"x", CARQUET_PHYSICAL_DOUBLE},
		{"y", CARQUET_PHYSICAL_DOUBLE},
		{"z", CARQUET_PHYSICAL_DOUBLE},
		{"stop", CARQUET_PHYSICAL_BOOLEAN},
	}};
	for (const auto &[name, type] : columns) {
		if (!CheckCarquetStatus(
			carquet_schema_add_column(
				schema,
				name,
				type,
				nullptr,
				CARQUET_REPETITION_REQUIRED,
				0,
				0
			),
			"Add parquet column"
		)) {
			carquet_schema_free(schema);
			return nullptr;
		}
	}
	return schema;
}

bool WriteParquetRows(
	carquet_writer_t *writer,
	const brill::T0Event &input_event,
	const brill::ParticleEvent &output_event
) {
	const int count = output_event.num;
	if (count == 0) return true;

	std::array<int32_t, 8> run;
	std::array<int32_t, 8> event_id;
	std::array<int32_t, 8> multiplicity;
	std::array<int32_t, 8> index;
	std::array<int32_t, 8> charge;
	std::array<int32_t, 8> mass;
	for (int i = 0; i < count; ++i) {
		run[i] = input_event.run;
		event_id[i] = input_event.entry;
		multiplicity[i] = count;
		index[i] = i;
		charge[i] = output_event.charge[i];
		mass[i] = output_event.mass[i];
	}

	const std::array<const void *, 11> values = {{
		event_id.data(), multiplicity.data(), run.data(), index.data(),
		charge.data(), mass.data(), output_event.energy, output_event.x,
		output_event.y, output_event.z, output_event.stop,
	}};
	for (int column = 0; column < int(values.size()); ++column) {
		if (!CheckCarquetStatus(
			carquet_writer_write_batch(writer, column, values[column], count, nullptr, nullptr),
			"Write parquet batch"
		)) return false;
	}
	return true;
}

int main(int argc, char **argv) {
	cxxopts::Options options("rebuild_t0", "Rebuild T0 particles from tracked events.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Run number.", cxxopts::value<int>(), "run")
		("t,trigger", "Trigger type.", cxxopts::value<std::string>(), "trigger")
		(
			"parquet",
			"Also save to parquet file.",
			cxxopts::value<bool>()->default_value("false")->implicit_value("true")
		)
		(
			"c,config",
			"Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"),
			"file"
		);

	auto result = options.parse(argc, argv);
	if (result.count("help")) {
		PrintUsage(options);
		return 0;
	}
	if (!result.count("run")) {
		std::cerr << "Error: Missing required option --run.\n";
		PrintUsage(options);
		return 1;
	}

	brill::AppConfig config;
	if (config.Load(result["config"].as<std::string>())) {
		return 1;
	}
	if (result.count("trigger")) {
		config.root.trigger = result["trigger"].as<std::string>();
	}

	const int run = result["run"].as<int>();
	if (config.IsJumpRun(run)) {
		std::cout << "Skipping jump run " << run << ".\n";
		return 0;
	}

	const std::string trigger_infix = brill::TriggerInfix(config.root.trigger);
	const std::string tele_dir = brill::JoinPath(config.root.workspace, config.paths.telescope);

	brill::CalibrationParameters calibration(2);
	const std::string calibration_path = TString::Format(
		"%s/t0.txt",
		brill::JoinPath(config.root.workspace, config.paths.calibration).c_str()
	).Data();
	if (calibration.Read(calibration_path)) {
		return 1;
	}

	TString input_path = TString::Format(
		"%s/t0_%s%04d.root",
		tele_dir.c_str(),
		trigger_infix.c_str(),
		run
	);
	TFile ipf(input_path, "read");
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from " << input_path << " failed.\n";
		return 1;
	}
	brill::T0Event input_event;
	brill::SetupInput(ipt, input_event);

	const std::string particle_dir = brill::JoinPath(config.root.workspace, config.paths.particle);
	std::filesystem::path particle_path(particle_dir);
	if (!particle_path.empty()) {
		std::filesystem::create_directories(particle_path);
	}
	const std::string output_path = TString::Format(
		"%s/t0_%s%04d.root",
		particle_dir.c_str(),
		trigger_infix.c_str(),
		run
	).Data();
	TFile opf(output_path.c_str(), "recreate");
	TTree opt("tree", "rebuild t0 particles");
	brill::ParticleEvent output_event;
	brill::Reset(output_event);
	brill::SetupOutput(&opt, output_event);

	carquet_schema_t *parquet_schema = nullptr;
	carquet_writer_t *parquet_writer = nullptr;
	if (result["parquet"].as<bool>()) {
		const std::string parquet_path = TString::Format(
			"%s/t0_%s%04d.par",
			particle_dir.c_str(),
			trigger_infix.c_str(),
			run
		).Data();
		parquet_schema = CreateParquetSchema();
		if (!parquet_schema) return 1;

		carquet_error_t error{};
		carquet_error_init(&error);
		parquet_writer = carquet_writer_create(
			parquet_path.c_str(), parquet_schema, nullptr, &error
		);
		if (!parquet_writer) {
			std::cerr << "Error: Create parquet file " << parquet_path
				<< " failed: " << error.message << ".\n";
			carquet_schema_free(parquet_schema);
			return 1;
		}
	}

	const long long total = ipt->GetEntries();
	long long last_percentage = -1;
	std::printf("Rebuilding T0   0%%");
	std::fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}

		ipt->GetEntry(entry);
		brill::Reset(output_event);
		int &num = output_event.num;
		for (int i = 0; i < input_event.num; ++i) {
			output_event.charge[num] = input_event.charge[i];
			output_event.mass[num] = input_event.mass[i];
			output_event.energy[num] =
				CalibrateEnergy(calibration, 0, input_event.energy[i][0])
				+ CalibrateEnergy(calibration, 1, input_event.energy[i][1]);
			output_event.time[num] = input_event.time[i][0];
			output_event.x[num] = input_event.x[i][0];
			output_event.y[num] = input_event.y[i][0];
			output_event.z[num] = input_event.z[i][0];
			output_event.stop[num] = true;
			++num;
		}

		opt.Fill();
		if (parquet_writer && !WriteParquetRows(parquet_writer, input_event, output_event)) {
			carquet_writer_abort(parquet_writer);
			carquet_schema_free(parquet_schema);
			opf.Close();
			ipf.Close();
			return 1;
		}
	}
	std::printf("\b\b\b\b100%%\n");

	opf.cd();
	opt.Write();
	opf.Close();
	ipf.Close();
	if (parquet_writer) {
		if (!CheckCarquetStatus(carquet_writer_close(parquet_writer), "Close parquet file")) {
			carquet_schema_free(parquet_schema);
			return 1;
		}
		carquet_schema_free(parquet_schema);
	}
	return 0;
}
