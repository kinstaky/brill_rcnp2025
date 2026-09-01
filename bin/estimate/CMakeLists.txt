add_executable(estimate_t0_pid estimate_t0_pid.cpp)
target_link_libraries(
	estimate_t0_pid
	PRIVATE
	brill::core brill::t0 ROOT::Graf t0_utils
)