 cmake_minimum_required(VERSION 2.8)
PROJECT (sceneVtk)

INCLUDE (${CMAKE_ROOT}/Modules/FindVTK.cmake)
IF (USE_VTK_FILE)
  INCLUDE(${USE_VTK_FILE})
ENDIF (USE_VTK_FILE)

ADD_EXECUTABLE(main main.cxx)
TARGET_LINK_LIBRARIES(main vtkRendering)
