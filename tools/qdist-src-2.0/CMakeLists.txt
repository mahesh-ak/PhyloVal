CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


SET(includeBlasFile include_blas.hpp)

IF(APPLE)

  PROJECT(QDist2 CXX C)

  SET(BLAS_LIBRARIES "")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=core2")
  FILE(WRITE ${includeBlasFile} "#include<vecLib/vBLAS.h>")

  SET(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -framework vecLib")

ELSE(APPLE)

  PROJECT(QDist2 CXX C Fortran)

  FIND_PACKAGE(BLAS REQUIRED)
  FILE(WRITE ${includeBlasFile} "#include<cblas.h>")

  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

ENDIF(APPLE)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -Wall")



SET(SOURCE_FILES
  DirectedEdge.hpp
  InternalNode.hpp
  LeafNode.hpp
  Matrix.hpp
  NewickParser.hpp
  NewickParser.cpp
  Node.hpp
  QDist.hpp
  QDist.cpp
  Tree.hpp
  TreeUtil.hpp
  TreeUtil.cpp
  Util.hpp
  Util.cpp)



ADD_EXECUTABLE(               qdist main.cpp         ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(        qdist                  ${BLAS_LIBRARIES})
INSTALL(TARGETS qdist RUNTIME DESTINATION bin)

ENABLE_TESTING()

ADD_EXECUTABLE(testMatrix testMatrix.cpp Matrix.hpp)
TARGET_LINK_LIBRARIES(testMatrix ${BLAS_LIBRARIES})
ADD_TEST(testMatrix testMatrix)

ADD_EXECUTABLE(testQDist testQDist.cpp ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(testQDist ${BLAS_LIBRARIES})
ADD_TEST(testQDist testQDist)



SET(zipFileName ${CMAKE_CURRENT_BINARY_DIR}/qdist-src.zip)
SET(zipFileContents
  ${SOURCE_FILES}
  main.cpp
  cubic-main.cpp
  quartic-main.cpp
  testMatrix.cpp
  testQDist.cpp
  CMakeLists.txt
  README
  testdata/*)

ADD_CUSTOM_COMMAND(
  OUTPUT  ${zipFileName}
  COMMAND zip -9 ${zipFileName} ${zipFileContents}
  COMMENT "Zipping QDist to ${zipFileName}"
)

ADD_CUSTOM_TARGET(
  zip
  DEPENDS ${zipFileName}
)
