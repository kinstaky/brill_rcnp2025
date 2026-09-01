add_library(pca pca.cpp)
target_include_directories(pca PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(pca PUBLIC Eigen3::Eigen)

add_library(normalize_extra normalize_extra.cpp)
target_include_directories(
	normalize_extra PUBLIC
	${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS}
)
target_link_libraries(normalize_extra PUBLIC pca)

add_library(strip_combination strip_combination.cpp)
target_include_directories(
	strip_combination PUBLIC
	${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS}
)
target_link_libraries(strip_combination PUBLIC normalize_extra)

add_library(t0_utils t0_utils.cpp)
target_include_directories(t0_utils PUBLIC ${PROJECT_SOURCE_DIR})