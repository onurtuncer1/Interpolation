#-------------------------------------------------------------------------------------------
#  Copyright 2022-2024 MILTEKSAN
#
#  Use of this software is restricted to MILTEKSAN
#
#--------------------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.22)

project(interpolation LANGUAGES CXX)

set(SRC_DIR "${PROJECT_SOURCE_DIR}/src")
set(3RD_PARTY_DIR "${PROJECT_SOURCE_DIR}/vendor")

message(STATUS "Third party directory:  " ${3RD_PARTY_DIR})



# ----------------Looking for git and updating submodules---------------------------------------
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")

# Git version
execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
                WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                OUTPUT_VARIABLE PACKAGE_GIT_VERSION
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE)
# Update submodules as needed
if(GIT_SUBMODULE)
    message(STATUS "Submodule update")
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
        message(WARNING "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}")
        endif()
    endif()
endif()

#--------------------------------------------------------------------------------------------------

# Remove ruckig tests, examples and cloud client which are on by default in its own CMakeLists.txt
set(BUILD_EXAMPLES OFF CACHE BOOL ""      FORCE)
set(BUILD_TESTS OFF CACHE BOOL ""         FORCE)
set(BUILD_CLOUD_CLIENT  OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS  OFF CACHE BOOL ""  FORCE)

add_subdirectory(${3RD_PARTY_DIR}/ruckig)

FILE(GLOB SOURCE_FILES ${SRC_DIR}/*.cpp)

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES})

set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX ""
                                                  OUTPUT_NAME ${PROJECT_NAME}
                                                  CXX_STANDARD 20
                                                  CXX_STANDARD_REQUIRED ON
                                                  CXX_EXTENSIONS OFF)

target_link_libraries(${PROJECT_NAME} PRIVATE ruckig)

target_compile_options(${PROJECT_NAME} PUBLIC -DRELEASE
                                               -Wall
                                               -Werror
                                               -Wpedantic)

target_include_directories(${PROJECT_NAME} PUBLIC  ${3RD_PARTY_DIR}/glm)
                                               #    ${3RD_PARTY_DIR}/ruckig/include)                                                  
