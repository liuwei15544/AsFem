cmake_minimum_required(VERSION 3.8)
project(AsFem)

set(CMAKE_CXX_STANDARD 14)
#######################################
### For debug mode                  ###
#######################################
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fopenmp -O2 ")

#######################################
### For release mode                ###
#######################################
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -Ofast -march=native")

#######################################
### For release mode                ###
#######################################
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
#######################################
### Set the PETSc path with your own###
#######################################
set(PETSc "/home/by/Programs/petsc3.9.3")
include_directories("${PETSc}/include")
link_libraries("${PETSc}/lib/libpetsc.so")
#######################################
### Set the mpi path with your own  ###
#######################################
set(MPI "/opt/moose/mpich-3.2/gcc-7.3.0")
include_directories("${MPI}/include")
link_libraries("${MPI}/lib/libmpi.so")



################################################################
### In most cases, you don't need to edit the following line!###
### In most cases, you don't need to edit the following line!###
### In most cases, you don't need to edit the following line!###
################################################################
#######################################
### For src and header file         ###
#######################################
set(SRC src/main.cpp)

# For Utils
#  for string utils
set(INC ${INC} include/Utils/StringUtils.h)
set(SRC ${SRC} src/Utils/StringUtils/RemoveSpace.cpp src/Utils/StringUtils/StrToLower.cpp)
set(SRC ${SRC} src/Utils/StringUtils/StrVecToLower.cpp src/Utils/StringUtils/StrToUpper.cpp)
set(SRC ${SRC} src/Utils/StringUtils/StrVecToUpper.cpp src/Utils/StringUtils/SplitStr.cpp)
set(SRC ${SRC} src/Utils/StringUtils/SplitNum.cpp src/Utils/StringUtils/SplitNumAfter.cpp)
set(SRC ${SRC} src/Utils/StringUtils/IsBracketMatch.cpp src/Utils/StringUtils/IsSubBracketMatch.cpp)
set(SRC ${SRC} src/Utils/StringUtils/SplitStrFromBracket.cpp src/Utils/StringUtils/GotoLine.cpp)


# For mesh
set(INC include/Mesh/Mesh.h)
set(SRC ${SRC} src/Mesh/Mesh.cpp src/Mesh/MeshSettings.cpp)
set(SRC ${SRC} src/Mesh/CreateMesh.cpp src/Mesh/Create1DMesh.cpp src/Mesh/Create2DMesh.cpp)
set(SRC ${SRC} src/Mesh/GetMeshInfo.cpp src/Mesh/Create3DMesh.cpp)
set(SRC ${SRC} src/Mesh/PrintAndSaveMesh.cpp)

# For dof handler
set(INC ${INC} include/DofHandler/DofHandler.h)
set(SRC ${SRC} src/DofHandler/DofHandler.cpp src/DofHandler/DofHandlerSetting.cpp)
set(SRC ${SRC} src/DofHandler/GetDofHandlerInfo.cpp src/DofHandler/PrintDofHandlerInfo.cpp)

# For equation system
set(INC ${INC} include/EquationSystem/EquationSystem.h)
set(SRC ${SRC} src/EquationSystem/EquationSystem.cpp src/EquationSystem/EquationSystemSetting.cpp)
set(SRC ${SRC} src/EquationSystem/EquationSystemInitAndReinit.cpp src/EquationSystem/EquationSystemMathUtils.cpp)
set(SRC ${SRC} src/EquationSystem/GetEquationSystemInfo.cpp)

# For boundary system
set(INC ${INC} include/BCSystem/BCBlockInfo.h include/BCSystem/BCInfo.h)
set(SRC ${SRC} src/BCSystem/BCInfo.cpp)
set(INC ${INC} include/BCSystem/BCSystem.h)
set(SRC ${SRC} src/BCSystem/BCSystem.cpp src/BCSystem/ApplyBoundaryCondition.cpp)
set(SRC ${SRC} src/BCSystem/ApplyNeumannBC.cpp src/BCSystem/ApplyDirichletBC.cpp src/BCSystem/ApplyConstraint.cpp)


# For input system
set(INC ${INC} include/InputSystem/InputSystem.h)
set(SRC ${SRC} src/InputSystem/InputSystem.cpp src/InputSystem/ReadMeshBlock.cpp)
set(SRC ${SRC} src/InputSystem/ReadDofsName.cpp src/InputSystem/ReadKernelBlock.cpp src/InputSystem/ReadBoundaryBlock.cpp)

# For element system (AsFem's uel)
set(INC ${INC} include/ElementSystem/KernelBlockInfo.h)
set(SRC ${SRC} src/ElementSystem/KernelBlockInfo.cpp)
set(INC ${INC} include/ElementSystem/ElementSystem.h)
set(SRC ${SRC} src/ElementSystem/ElementSystem.cpp src/ElementSystem/ElementSystemSetting.cpp)
set(SRC ${SRC} src/ElementSystem/RunElmtLib.cpp src/ElementSystem/elements/Poisson.cpp)

# For shapefuns
set(INC ${INC} include/ShapeFuns/ShapeFuns.h)
set(SRC ${SRC} src/ShapeFuns/Shp1D.cpp src/ShapeFuns/Shp2D.cpp)


add_executable(ASFEM ${INC} ${SRC})