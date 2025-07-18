# Copyright (c), ETH Zurich and UNC Chapel Hill.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#
#     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
#       its contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.


# Find package module for Metis library.
#
# The following variables are set by this module:
#
#   METIS_FOUND: TRUE if Metis is found.
#   metis: Imported target to link against.
#
# The following variables control the behavior of this module:
#
# METIS_INCLUDE_DIR_HINTS: List of additional directories in which to
#                              search for Metis includes.
# METIS_LIBRARY_DIR_HINTS: List of additional directories in which to
#                              search for Metis libraries.

set(METIS_INCLUDE_DIR_HINTS "" CACHE PATH "Metis include directory")
set(METIS_LIBRARY_DIR_HINTS "" CACHE PATH "Metis library directory")

unset(METIS_FOUND)

find_package(metis CONFIG QUIET)
if(TARGET metis)
    set(METIS_FOUND TRUE)
    message(STATUS "Found Metis")
    message(STATUS "  Target : metis")
else()
    list(APPEND METIS_CHECK_INCLUDE_DIRS
        ${METIS_INCLUDE_DIR_HINTS}
        /usr/include
        /usr/local/include
        /opt/include
        /opt/local/include
    )

    list(APPEND METIS_CHECK_LIBRARY_DIRS
        ${METIS_LIBRARY_DIR_HINTS}
        /usr/lib
        /usr/local/lib
        /opt/lib
        /opt/local/lib
    )

    find_path(METIS_INCLUDE_DIRS
        NAMES
        metis.h
        PATHS
        ${METIS_CHECK_INCLUDE_DIRS})
    find_library(METIS_LIBRARIES
        NAMES
        metis
        PATHS
        ${METIS_CHECK_LIBRARY_DIRS})
    find_library(GK_LIBRARIES
        NAMES
        GKlib
        PATHS
        ${METIS_CHECK_LIBRARY_DIRS})

    if(GK_LIBRARIES)
        set(METIS_LIBRARIES ${METIS_LIBRARIES} ${GK_LIBRARIES})
        message(STATUS "Found GKlib")
    endif()

    if(METIS_INCLUDE_DIRS AND METIS_LIBRARIES)
        set(METIS_FOUND TRUE)
        message(STATUS "Found Metis")
        message(STATUS "  Includes : ${METIS_INCLUDE_DIRS}")
        message(STATUS "  Libraries : ${METIS_LIBRARIES}")
    endif()

    add_library(metis INTERFACE IMPORTED)
    target_include_directories(
        metis INTERFACE ${METIS_INCLUDE_DIRS})
    target_link_libraries(
        metis INTERFACE ${METIS_LIBRARIES})
endif()

if(NOT METIS_FOUND AND METIS_FIND_REQUIRED)
    message(FATAL_ERROR "Could not find Metis")
endif()
