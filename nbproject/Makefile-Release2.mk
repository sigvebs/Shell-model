#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_CONF=Release2
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/includes/lib.o \
	${OBJECTDIR}/Lanczos.o \
	${OBJECTDIR}/HamiltonMatrix.o \
	${OBJECTDIR}/ImportSlater.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-m64 -O3 -DMPICH_IGNORE_CXX_SEEK -llapack -lblas -larmadillo
CXXFLAGS=-m64 -O3 -DMPICH_IGNORE_CXX_SEEK -llapack -lblas -larmadillo

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-larmadillo

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/shellmodel

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/shellmodel: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/shellmodel ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O3 -w -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/includes/lib.o: includes/lib.cpp 
	${MKDIR} -p ${OBJECTDIR}/includes
	${RM} $@.d
	$(COMPILE.cc) -O3 -w -MMD -MP -MF $@.d -o ${OBJECTDIR}/includes/lib.o includes/lib.cpp

${OBJECTDIR}/Lanczos.o: Lanczos.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O3 -w -MMD -MP -MF $@.d -o ${OBJECTDIR}/Lanczos.o Lanczos.cpp

${OBJECTDIR}/HamiltonMatrix.o: HamiltonMatrix.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O3 -w -MMD -MP -MF $@.d -o ${OBJECTDIR}/HamiltonMatrix.o HamiltonMatrix.cpp

${OBJECTDIR}/ImportSlater.o: ImportSlater.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O3 -w -MMD -MP -MF $@.d -o ${OBJECTDIR}/ImportSlater.o ImportSlater.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/shellmodel

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
