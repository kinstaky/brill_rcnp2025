add_executable(estimate_t0_pid estimate_t0_pid.cpp)
target_link_libraries(
	estimate_t0_pid
	PRIVATE
	brill::core brill::t0 ROOT::Graf t0_utils
)

add_executable(estimate_gagg_trace estimate_gagg_trace.cpp)
target_link_libraries(
	estimate_gagg_trace
	PRIVATE
	brill::core ROOT::Graf ROOT::Tree
)