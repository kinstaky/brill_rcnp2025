#include "include/strip_combination.h"

#include <cmath>
#include <vector>

namespace brill {

StripCombination::StripCombination(const int type, const int flag)
: type_(type), flag_(flag) {}


MatchResult StripCombination::Match(
	const StripCombination &other,
	const double threshold,
	const bool print
) const {
	MatchResult result;
	result.valid = false;
	std::vector<EnergyGuess> this_guesses = Guesses();
	std::vector<EnergyGuess> other_guesses = other.Guesses();
	for (size_t i = 0; i < this_guesses.size(); ++i) {
		for (size_t j = 0; j < other_guesses.size(); ++j) {
			double guess = WeightedGuessValue(this_guesses[i], other_guesses[j]);

			// for debug
			if (print) {
				std::cout << "Combination A: Guess " << i
					<< ", guess " << this_guesses[i].value << ", type " << this_guesses[i].type
					<< ", flag " << this_guesses[i].flag << std::endl;
				std::cout << "Combination B: Guess " << j
					<< ", guess " << other_guesses[j].value << ", type " << other_guesses[j].type
					<< ", flag " << other_guesses[j].flag << std::endl;
				std::cout << "Weighted guess value: " << guess << std::endl;
			}

			double distance = fabs(this_guesses[i].value - other_guesses[j].value);
			if (distance > threshold) continue;
			if (result.valid && distance >= result.distance) continue;
			result.valid = true;
			result.flag = this_guesses[i].flag | other_guesses[j].flag;
			result.type[0] = this_guesses[i].type;
			result.type[1] = other_guesses[j].type;
			result.strip[0] = Strip();
			result.strip[1] = other.Strip();
			result.energy = guess;
			result.time = this_guesses[i].time;
			result.distance = distance;
		}
	}
	return result;
}

/*
 * NormalStripCombination
 */
NormalStripCombination::NormalStripCombination(
	const int flag,
	const StripInfo &strip,
	const double threshold
)
: StripCombination(StripType::Normal, flag)
, energy_(
	strip.p0 + strip.p1*strip.raw_energy
	+ strip.p2*strip.raw_energy*strip.raw_energy
	+ strip.p3*strip.raw_energy*strip.raw_energy*strip.raw_energy
)
, strip_(strip.strip)
, time_(strip.time)
, threshold_(threshold) {
}

std::vector<EnergyGuess> NormalStripCombination::Guesses() const {
	std::vector<EnergyGuess> result;
	result.push_back(EnergyGuess{energy_, 1.0, time_, type_, flag_});
	return result;
}

double NormalStripCombination::Distance(const double value, const int) const {
	return fabs(value - energy_);
}

int NormalStripCombination::Strip() const {
	return strip_;
}

std::ostream &operator<<(std::ostream &os, const NormalStripCombination &strip) {
	os << strip.energy_ << " " << strip.strip_
		<< " " << strip.time_ << " " << strip.threshold_;
	return os;
}


/*
 * NormalSharedStripCombination
 */
NormalSharedStripCombination::NormalSharedStripCombination(
	const int flag,
	const StripInfo &strip1,
	const StripInfo &strip2,
	const double threshold
)
: StripCombination(StripType::Normal | StripType::ChargeShared, flag)
, time_(strip1.time)
, threshold_(threshold) {
	energy_[0] = strip1.p0
		+ strip1.p1*strip1.raw_energy
		+ strip1.p2*strip1.raw_energy*strip1.raw_energy
		+ strip1.p3*strip1.raw_energy*strip1.raw_energy*strip1.raw_energy;
	energy_[1] = strip2.p0
		+ strip2.p1*strip2.raw_energy
		+ strip2.p2*strip2.raw_energy*strip2.raw_energy
		+ strip2.p3*strip2.raw_energy*strip2.raw_energy*strip2.raw_energy;
	strip_ = energy_[0] > energy_[1] ? strip1.strip : strip2.strip;
}

std::vector<EnergyGuess> NormalSharedStripCombination::Guesses() const {
	std::vector<EnergyGuess> result;
	result.push_back(EnergyGuess{energy_[0]+energy_[1], 2.0, time_, type_, flag_});
	return result;
}

double NormalSharedStripCombination::Distance(const double value, const int) const {
	return fabs(value - energy_[0] - energy_[1]);
}

int NormalSharedStripCombination::Strip() const {
	return strip_;
}


/*
 * BrokenAdjacentStripCombination
 */
BrokenAdjacentStripCombination::BrokenAdjacentStripCombination(
	const int flag,
	const StripInfo &strip1,
	const StripInfo &strip2,
	const int broken_strip,
	const Line3D &line,
	const Plane &plane1,
	const Plane &plane2,
	const double line_threshold,
	const double plane1_threshold,
	const double plane2_threshold
)
: StripCombination(StripType::BrokenAdjacent, flag)
, broken_strip_(broken_strip)
, line_(line)
, line_threshold_(line_threshold)
, time_(strip1.time) {
	raw_energy_[0] = strip1.raw_energy;
	raw_energy_[1] = strip2.raw_energy;
	plane_[0] = plane1;
	plane_[1] = plane2;
	plane_threshold_[0] = plane1_threshold;
	plane_threshold_[1] = plane2_threshold;
}

std::vector<EnergyGuess> BrokenAdjacentStripCombination::Guesses() const {
	std::vector<EnergyGuess> result;
	const double &xm = raw_energy_[0];
	const double &ym = raw_energy_[1];

	// line
	const double &lcx = line_.centroid(0);
	const double &lcy = line_.centroid(1);
	const double &lcz = line_.centroid(2);
	const double &ex = line_.direction(0);
	const double &ey = line_.direction(1);
	const double &ez = line_.direction(2);
	double line_value = lcz + (ex*ez*(xm-lcx) + ey*ez*(ym-lcy)) / (1.0 - ez*ez);
	result.push_back(EnergyGuess{line_value, ez*ez/(1.0 - ez*ez), time_, type_, flag_});

	// plane
	for (int i = 0; i < 2; ++i) {
		const double &cx = plane_[i].centroid(0);
		const double &cy = plane_[i].centroid(1);
		const double &cz = plane_[i].centroid(2);
		const double &nx = plane_[i].normal(0);
		const double &ny = plane_[i].normal(1);
		const double &nz = plane_[i].normal(2);
		double plane_value = nx/nz*(cx-xm) + ny/nz*(cy-ym) + cz;
		result.push_back(EnergyGuess{
			plane_value,
			(1.0-nz*nz)/(nz*nz),
			time_,
			type_ | StripType::ChargeShared,
			flag_
		});
	}

	return result;
}

double BrokenAdjacentStripCombination::Distance(
	const double value,
	const int piece
) const {
	if (piece == 0) return LineDistance(value);
	else if (piece == 1) return PlaneDistance(value, 0);
	else return PlaneDistance(value, 1);
}

double BrokenAdjacentStripCombination::LineDistance(const double value) const {
	return line_.DistanceFromPoint(
		Eigen::Vector3d(raw_energy_[0], raw_energy_[1], value)
	);
}

double BrokenAdjacentStripCombination::PlaneDistance(
	const double value,
	const int index
) const {
	return plane_[index].DistanceFromPoint(
		Eigen::Vector3d(raw_energy_[0], raw_energy_[1], value)
	);
}

int BrokenAdjacentStripCombination::Strip() const {
	return broken_strip_;
}

// bool BrokenAdjacentStripCombination::VerifyGuess(
// 	const double value,
// 	const int piece
// ) const {
// 	if (piece == 0) return LineDistance(value) < line_threshold_;
// 	else if (piece == 1) return PlaneDistance(value, 0) < plane_threshold_[0];
// 	else return PlaneDistance(value, 1) < plane_threshold_[1];
// }


/*
 * ShortStripCombination
 */
ShortStripCombination::ShortStripCombination(
	const int flag,
	const StripInfo &strip1,
	const StripInfo &strip2,
	const Line3D &line,
	const double threshold
)
: StripCombination(StripType::Short, flag)
, strip_(strip1.strip)
, time_(strip1.time)
, line_(line)
, threshold_(threshold) {
	raw_energy_[0] = strip1.raw_energy;
	raw_energy_[1] = strip2.raw_energy;
}

std::vector<EnergyGuess> ShortStripCombination::Guesses() const {
	std::vector<EnergyGuess> result;
	const double &xm = raw_energy_[0];
	const double &ym = raw_energy_[1];
	const double &cx = line_.centroid(0);
	const double &cy = line_.centroid(1);
	const double &cz = line_.centroid(2);
	const double &ex = line_.direction(0);
	const double &ey = line_.direction(1);
	const double &ez = line_.direction(2);
	result.push_back(EnergyGuess{
		cz + (ex*ez*(xm-cx) + ey*ez*(ym-cy)) / (1.0 - ez*ez),
		ez*ez/(1.0 - ez*ez),
		time_,
		type_,
		flag_
	});
	return result;
}

double ShortStripCombination::Distance(const double value, const int) const {
	return line_.DistanceFromPoint(Eigen::Vector3d(raw_energy_[0], raw_energy_[1], value));
}

int ShortStripCombination::Strip() const {
	return strip_;
}

// bool ShortStripCombination::VerifyGuess(const double value, const int) const {
// 	return Distance(value) < threshold_;
// }


/*
 * ShortStripSharedCombination
 */
ShortSharedStripCombination::ShortSharedStripCombination(
	const int flag,
	const StripInfo &strip1,
	const StripInfo &strip2,
	const StripInfo &strip3,
	const Line3D &line,
	const double threshold
)
: StripCombination(StripType::Short | StripType::ChargeShared, flag)
, shared_energy_(
	strip3.p0
	+ strip3.p1 * strip3.raw_energy
	+ strip3.p2 * strip3.raw_energy * strip3.raw_energy
	+ strip3.p3 * strip3.raw_energy * strip3.raw_energy * strip3.raw_energy
)
, shared_strip_(strip3.strip)
, time_(strip1.time)
, line_(line)
, threshold_(threshold) {
	short_energy_[0] = strip1.raw_energy;
	short_energy_[1] = strip2.raw_energy;
}

std::vector<EnergyGuess> ShortSharedStripCombination::Guesses() const {
	std::vector<EnergyGuess> result;
	const double &xm = short_energy_[0];
	const double &ym = short_energy_[1];
	const double &cx = line_.centroid(0);
	const double &cy = line_.centroid(1);
	const double &cz = line_.centroid(2);
	const double &ex = line_.direction(0);
	const double &ey = line_.direction(1);
	const double &ez = line_.direction(2);
	double line_value = cz + (ex*ez*(xm-cx) + ey*ez*(ym-cy)) / (1.0 - ez*ez);
	result.push_back(EnergyGuess{
		line_value + shared_energy_,
		ez*ez/(1.0 - ez*ez) + 1.0,
		time_,
		type_,
		flag_
	});
	return result;
}

double ShortSharedStripCombination::Distance(const double value, const int) const {
	return line_.DistanceFromPoint(
		Eigen::Vector3d(short_energy_[0], short_energy_[1], value-shared_energy_)
	);
}

int ShortSharedStripCombination::Strip() const {
	return shared_strip_;
}

// bool ShortSharedStripCombination::VerifyGuess(const double value, const int) const {
// 	return Distance(value) < threshold_;
// }


/*
 * PiecewiseStripCombination
 */
PiecewiseStripCombination::PiecewiseStripCombination(
	const int flag,
	const StripInfo &strip,
	const PiecewiseParameter &piecewise,
	const double threshold
)
: StripCombination(StripType::Piecewise, flag)
, raw_energy_(strip.raw_energy)
, strip_(strip.strip)
, time_(strip.time)
, piecewise_(piecewise)
, threshold_(threshold) {
}

std::vector<EnergyGuess> PiecewiseStripCombination::Guesses() const {
	std::vector<EnergyGuess> result;
	for (size_t i = 0; i < piecewise_.range.size(); ++i) {
		const Range &rg = piecewise_.range[i];
		if (rg.xmin > 0 && raw_energy_ < rg.xmin) continue;
		if (rg.xmax > 0 && raw_energy_ >= rg.xmax) continue;
		double guess = piecewise_.pol2[i].p0
			+ piecewise_.pol2[i].p1 * raw_energy_
			+ piecewise_.pol2[i].p2 * raw_energy_ * raw_energy_;
		if (rg.ymin > 0 && guess < rg.ymin) continue;
		if (rg.ymax > 0 && guess > rg.ymax) continue;
		result.push_back(EnergyGuess{guess, 1.0, time_, type_, flag_});
	}
	return result;
}

double PiecewiseStripCombination::Distance(
	const double value,
	const int piece
) const {
	double compare_value = piecewise_.pol2[piece].p0
		+ piecewise_.pol2[piece].p1 * raw_energy_
		+ piecewise_.pol2[piece].p2 * raw_energy_ * raw_energy_;
	return fabs(value - compare_value);
}

int PiecewiseStripCombination::Strip() const {
	return strip_;
}

// bool PiecewiseStripCombination::VerifyGuess(
// 	const double value,
// 	const int piece
// ) const {
// 	return Distance(value, piece) < threshold_;
// }


/*
 * PiecewiseStripSharedCombination
 */
PiecewiseSharedStripCombination::PiecewiseSharedStripCombination(
	const int flag,
	const StripInfo &piecewise_strip,
	const StripInfo &shared_strip,
	const PiecewiseParameter &piecewise,
	const double threshold
)
: StripCombination(StripType::Piecewise | StripType::ChargeShared, flag)
, raw_energy_(piecewise_strip.raw_energy)
, shared_energy_(
	shared_strip.p0
	+ shared_strip.p1 * shared_strip.raw_energy
	+ shared_strip.p2 * shared_strip.raw_energy * shared_strip.raw_energy
	+ shared_strip.p3 * shared_strip.raw_energy * shared_strip.raw_energy * shared_strip.raw_energy
)
, strip_(piecewise_strip.strip)
, time_(piecewise_strip.time)
, piecewise_(piecewise)
, threshold_(threshold) {
}

std::vector<EnergyGuess> PiecewiseSharedStripCombination::Guesses() const {
	std::vector<EnergyGuess> result;
	for (size_t i = 0; i < piecewise_.range.size(); ++i) {
		const Range &rg = piecewise_.range[i];
		if (rg.xmin > 0 && raw_energy_ < rg.xmin) continue;
		if (rg.xmax > 0 && raw_energy_ >= rg.xmax) continue;
		double guess = piecewise_.pol2[i].p0
			+ piecewise_.pol2[i].p1 * raw_energy_
			+ piecewise_.pol2[i].p2 * raw_energy_ * raw_energy_;
		if (rg.ymin > 0 && guess < rg.ymin) continue;
		if (rg.ymax > 0 && guess > rg.ymax) continue;
		result.push_back(EnergyGuess{guess + shared_energy_, 2.0, time_, type_, flag_});
	}
	return result;
}

double PiecewiseSharedStripCombination::Distance(
	const double value,
	const int piece
) const {
	double compare_value = piecewise_.pol2[piece].p0
		+ piecewise_.pol2[piece].p1 * raw_energy_
		+ piecewise_.pol2[piece].p2 * raw_energy_ * raw_energy_;
	return fabs(value - compare_value - shared_energy_);
}

int PiecewiseSharedStripCombination::Strip() const {
	return strip_;
}

// bool PiecewiseSharedStripCombination::VerifyGuess(
// 	const double value,
// 	const int piece
// ) const {
// 	return Distance(value, piece) < threshold_;
// }

}