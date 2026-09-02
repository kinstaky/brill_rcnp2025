add_subdirectory(estimate)

add_executable(normalize_extra_t0d1 normalize_extra_t0d1.cpp)
target_link_libraries(
	normalize_extra_t0d1
	PRIVATE
	brill::core brill::t0 ROOT::Graf Eigen3::Eigen normalize_extra
)

add_executable(normalize_extra_t0d2 normalize_extra_t0d2.cpp)
target_link_libraries(
	normalize_extra_t0d2
	PRIVATE
	brill::core brill::t0 ROOT::Graf Eigen3::Eigen normalize_extra
)

add_executable(match_dssd match_dssd.cpp)
target_link_libraries(
	match_dssd PRIVATE
	brill::core brill::t0 t0_utils ROOT::Graf Eigen3::Eigen normalize_extra strip_combination
)

add_executable(track_t0 track_t0.cpp)
target_link_libraries(
	track_t0 PRIVATE
	brill::core brill::t0 t0_utils ROOT::Graf
)

add_executable(calibrate_t0 calibrate_t0.cpp)
target_link_libraries(
	calibrate_t0 PRIVATE
	brill::core brill::t0 brill::energy_calculator ROOT::Graf t0_utils
)

add_executable(rebuild_t0 rebuild_t0.cpp)
target_link_libraries(
	rebuild_t0 PRIVATE
	brill::core brill::t0 brill::energy_calculator ROOT::Graf t0_utils brill::particle carquet
)
