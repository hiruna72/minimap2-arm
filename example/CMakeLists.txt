# Get all the directories with needed .hpp files
include_directories ("${PROJECT_SOURCE_DIR}/src")

set_source_files_properties(minimap2_example.c PROPERTIES LANGUAGE CXX)

# Create the executable
add_executable(minimap2_example minimap2_example.c)

# Link the required libraries
target_link_libraries(minimap2_example minimapshared)