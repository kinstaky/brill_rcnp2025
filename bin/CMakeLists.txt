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

add_executable(match match.cpp)
target_link_libraries(
	match PRIVATE
	brill::core brill::t0 ROOT::Graf Eigen3::Eigen normalize_extra strip_combination
)