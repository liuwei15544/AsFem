cmake_minimum_required(VERSION 3.8)
project(AsFem)

set(CMAKE_CXX_STANDARD 14)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include/Eigen")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
##############################################
### For debug mode                         ###
##############################################
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O2 -Wall -Wextra")

##############################################
### For release mode                       ###
##############################################
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3 -Wall -march=native -DNDEBUG")

#####################################################################
### In most cases, you don't need to modify the following lines!!!###
### In most cases, you don't need to modify the following lines!!!###
### In most cases, you don't need to modify the following lines!!!###
#####################################################################
set(SRC src/main.cpp)
set(INC include/Welcome.h)

#################################################
### Message print for input step              ###
#################################################
set(INC ${INC} include/MsgPrint/MsgPrintForInput.h)
set(SRC ${SRC} src/MsgPrint/Input/Msg_InputFileNameWrong.cpp)

### For program
set(INC ${INC} include/MsgPrint/MsgPrintForProgram.h)
set(SRC ${SRC} src/MsgPrint/Program/Msg_ExitProgram.cpp)

#################################################
### String utils                              ###
#################################################
#  for string utils
set(INC ${INC} include/Utils/StringUtils.h)
set(SRC ${SRC} src/Utils/StringUtils/RemoveSpace.cpp src/Utils/StringUtils/StrToLower.cpp)
set(SRC ${SRC} src/Utils/StringUtils/StrVecToLower.cpp src/Utils/StringUtils/StrToUpper.cpp)
set(SRC ${SRC} src/Utils/StringUtils/StrVecToUpper.cpp src/Utils/StringUtils/SplitStr.cpp)
set(SRC ${SRC} src/Utils/StringUtils/SplitNum.cpp src/Utils/StringUtils/SplitNumAfter.cpp)
set(SRC ${SRC} src/Utils/StringUtils/IsBracketMatch.cpp src/Utils/StringUtils/IsSubBracketMatch.cpp)
set(SRC ${SRC} src/Utils/StringUtils/SplitStrFromBracket.cpp src/Utils/StringUtils/GotoLine.cpp)



#################################################
### Input system                              ###
#################################################
set(INC ${INC} include/InputSystem/InputSystem.h)
set(SRC ${SRC} src/InputSystem/InputSystem.cpp src/MsgPrint/Input/Msg_InputSystemInitWarning.cpp)

#################################################


add_executable(asfem ${INC} ${SRC})

