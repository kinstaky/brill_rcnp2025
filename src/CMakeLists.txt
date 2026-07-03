add_library(pca pca.cpp)
target_include_directories(pca PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(pca PUBLIC Eigen3::Eigen)

add_library(normalize_extra normalize_extra.cpp)
target_include_directories(
	normalize_extra PUBLIC
	${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS}
)
target_link_libraries(normalize_extra PUBLIC pca)