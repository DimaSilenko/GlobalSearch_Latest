cmake_minimum_required(VERSION 3.1) #�������� ������ CMake.

set(PROJECT_NAME 000_Lebedev_AllLab) #������� ������� ��������� ���������� � ������ �������
set(CMAKE_CXX_STANDARD 11)
project(${PROJECT_NAME}) # �������� �������

set(PROBLEM_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/_bin CACHE STRING "Output directory")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${PROBLEM_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PROBLEM_OUTPUT_DIRECTORY})

set(GlobalSearch GlobalSearch) #���������� � ������ ����������
set(pugixml pugixml) #���������� � ������ ����������
set(src src)
set(Drawing Drawing)
option(USE_MPI "Set OFF to disable MPI" OFF)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Configs" FORCE) #��������� ���������� CMAKE_CONFIGURATION_TYPES �� ��������� "Debug;Release" � ������/�������� Configs, ���������� ����������� � ����������� ����� ��������� �� ���� ������,
if(NOT CMAKE_BUILD_TYPE) # �������� ��������� � ����������� ��������� ������
  set(CMAKE_BUILD_TYPE Release) #���� ��� �� ������������� Release
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/pugixml)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Drawing)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/SRC/)

if(WIN32)
    option(gtest_force_shared_crt "" TRUE)
endif(WIN32)

#include( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Configure.cmake" )


################################ MPI ################################

if( ${USE_MPI} )

    find_package( MPI )
    if( MPI_FOUND )
	message ("MPI FOUND!!!")
        include_directories( ${MPI_INCLUDE_PATH} )
    else( MPI_FOUND )
		message ("MPI Not Found")
        set( USE_MPI OFF )
    endif( MPI_FOUND )
endif( USE_MPI )

add_subdirectory(pugixml)
add_subdirectory(Drawing)
add_subdirectory(src)
add_subdirectory(GlobalSearch)
add_subdirectory(MyFunction)
add_subdirectory(Rastrigin2)
add_subdirectory(GKLS)
add_subdirectory(Grishagin)
add_subdirectory(StronginC3)


#target_link_libraries(src ${EXAMIN_MPFR_LIB} ${EXAMIN_MPIR_LIB} ${DISLIN_LIB})


#-----------------------------------------------------------------------------------
#file(GLOB hdrs "*.h*") #��������� � ���������� hdrs ��� ����� � ����������� .h
#file(GLOB srcs "*.cpp") #��������� � ���������� srcs ��� ����� � ����������� .cpp

#include_directories(lib1)
#add_executable(${PROJECT_NAME} ${srcs} ${hdrs}) #������� ������ � ������ �� ���������� ${PROJECT_NAME}, � ���� ����������� ����� �� ���������� ${srcs} � ${hdrs}


 #�������� ������ � ������������� lib1
#target_link_libraries(${PROJECT_NAME} ${MyLib1}) #���������� � ������� � ������ ${PROJECT_NAME} ���������� � ������ ${MyLib1}