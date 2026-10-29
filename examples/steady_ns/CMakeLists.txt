# Copyright 2026 Lawrence Livermore National Security, LLC.
# SPDX-License-Identifier: MIT

set(STEADY_NS_EXAMPLE_FILES
    README.md
    generate_configs.py
    prepare_meshes.sh
    report_results.py
    steady_ns.yml
    submit_workflow.sh)

foreach(example_file IN LISTS STEADY_NS_EXAMPLE_FILES)
  configure_file(${example_file} ${example_file} COPYONLY)
endforeach()

file(COPY jobs DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

