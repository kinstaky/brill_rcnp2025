#pragma once

#include <vector>
#include <iostream>

#include "include/geometry3d.h"
#include "include/normalize_extra.h"

namespace brill {

namespace StripType {
	int Normal = 0;
	int ChargeShared = 1;
	int BrokenAdjacent = 2;
	int Short = 4;
	int Piecewise = 8;
}

struct MatchResult {
	bool valid = false;
	int flag = 0;
	int type[2] = {0, 0};
	int strip[2] = {-1, -1};
	double energy = 0.0;
	double distance = 0.0;
	double time = 0.0;
};

struct EnergyGuess {
	double value = 0.0;
	double weight = 0.0;
	double time = 0.0;
	int type = 0;
	int flag = 0;
};

struct StripInfo {
	int strip;
	double raw_energy = 0.0;
	double time = 0.0;
	double p0 = 0.0;
	double p1 = 1.0;
	double p2 = 0.0;
	double p3 = 0.0;
};

inline double WeightedGuessValue(
	const EnergyGuess &guess1,
	const EnergyGuess &guess2
) {
	const double &value1 = guess1.value;
	const double &value2 = guess2.value;
	const double &weight1 = guess1.weight;
	const double &weight2 = guess2.weight;
	return (weight2*value1 + weight1*value2) / (weight1 + weight2);
}

class StripCombination {
public:
	StripCombination(const int type, const int flag);
	virtual ~StripCombination() = default;

	virtual std::vector<EnergyGuess> Guesses() const = 0;

	virtual double Distance(const double value, const int piece) const = 0;

	virtual int Strip() const = 0;

	// virtual bool VerifyGuess(const double value) const = 0;

	virtual MatchResult Match(
		const StripCombination &other,
		const double threshold,
		const bool print = false
	) const;

protected:
	int type_;
	int flag_;
};


class NormalStripCombination: public StripCombination {
public:
	NormalStripCombination(
		const int flag,
		const StripInfo &strip,
		const double threshold
	);

	virtual ~NormalStripCombination() = default;

	virtual std::vector<EnergyGuess> Guesses() const override;
	virtual double Distance(const double value, const int piece) const override;
	virtual int Strip() const override;

	friend std::ostream &operator<<(
		std::ostream &os,
		const NormalStripCombination &strip
	);
private:
	double energy_;
	int strip_;
	double time_;
	double threshold_;
};


class NormalSharedStripCombination: public StripCombination {
public:
	NormalSharedStripCombination(
		const int flag,
		const StripInfo &strip1,
		const StripInfo &strip2,
		const double threshold
	);

	virtual ~NormalSharedStripCombination() = default;

	virtual std::vector<EnergyGuess> Guesses() const override;
	virtual double Distance(const double value, const int piece) const override;
	virtual int Strip() const override;

private:
	double energy_[2];
	int strip_;
	double time_;
	double threshold_;
};


class BrokenAdjacentStripCombination: public StripCombination {
public:
	BrokenAdjacentStripCombination(
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
	);

	virtual ~BrokenAdjacentStripCombination() = default;

	virtual std::vector<EnergyGuess> Guesses() const override;
	virtual double Distance(const double value, const int piece) const override;
	double LineDistance(const double value) const;
	double PlaneDistance(const double value, const int index) const;
	virtual int Strip() const override;

private:
	double raw_energy_[2];
	int broken_strip_;
	Line3D line_;
	double line_threshold_;
	Plane plane_[2];
	double plane_threshold_[2];
	double time_;
};

class ShortStripCombination: public StripCombination {
public:
	ShortStripCombination(
		const int flag,
		const StripInfo &strip1,
		const StripInfo &strip2,
		const Line3D &line,
		const double threshold
	);

	virtual ~ShortStripCombination() = default;

	virtual std::vector<EnergyGuess> Guesses() const override;
	virtual double Distance(const double value, const int piece) const override;
	virtual int Strip() const override;

private:
	double raw_energy_[2];
	int strip_;
	double time_;
	Line3D line_;
	double threshold_;
};

class ShortSharedStripCombination: public StripCombination {
public:
	ShortSharedStripCombination(
		const int flag,
		const StripInfo &strip1,
		const StripInfo &strip2,
		const StripInfo &strip3,
		const Line3D &line,
		const double threshold
	);

	virtual ~ShortSharedStripCombination() = default;

	virtual std::vector<EnergyGuess> Guesses() const override;
	virtual double Distance(const double value, const int piece) const override;
	virtual int Strip() const override;

private:
	double short_energy_[2];
	double shared_energy_;
	int shared_strip_;
	double time_;
	Line3D line_;
	double threshold_;
};

class PiecewiseStripCombination: public StripCombination {
public:
	PiecewiseStripCombination(
		const int flag,
		const StripInfo &strip,
		const PiecewiseParameter &piecewise,
		const double threshold
	);

	~PiecewiseStripCombination() = default;

	virtual std::vector<EnergyGuess> Guesses() const override;
	virtual double Distance(const double value, const int piece) const override;
	virtual int Strip() const override;

private:
	double raw_energy_;
	int strip_;
	double time_;
	PiecewiseParameter piecewise_;
	double threshold_;
};

class PiecewiseSharedStripCombination: public StripCombination {
public:
	PiecewiseSharedStripCombination(
		const int flag,
		const StripInfo &piecewise_strip,
		const StripInfo &shared_strip,
		const PiecewiseParameter &piecewise,
		const double threshold
	);

	virtual ~PiecewiseSharedStripCombination() = default;

	virtual std::vector<EnergyGuess> Guesses() const override;
	virtual double Distance(const double value, const int piece) const override;
	virtual int Strip() const override;

private:
	double raw_energy_, shared_energy_;
	int strip_;
	double time_;
	PiecewiseParameter piecewise_;
	double threshold_;
};

}