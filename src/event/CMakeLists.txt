add_library(gagg_event gagg_event.cpp)
target_include_directories(gagg_event PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(gagg_event PUBLIC ROOT::Core ROOT::Tree)