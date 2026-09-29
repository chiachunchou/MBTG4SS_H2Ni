TARGET_NAME=qtr_serial

#-----------------------------------------------------------

#C compiler settings

CC=gcc
CXX=g++
LD=ld
AR=ar cru
RANLIB=ranlib

ifndef DEBUG
        NDEBUG=1
        OPT=-O3 -fopenmp -flto -ffast-math
        ifeq ($(shell uname -m),aarch64)
                OPT += -mcpu=neoverse-v2
        else
                OPT += -march=native
        endif
endif

CXXFLAGS += -Wall -Wno-sign-compare -Wno-unused-function -Wfatal-errors $(OPT) -std=c++11

LDFLAGS += $(OPT)

#-----------------------------------------------------------

#MPI settings

ifdef QTRMPI
        CXX=mpicxx
        CXXFLAGS += -DQTRMPI
        TARGET_NAME=qtr_mpi
endif

#-----------------------------------------------------------

#Ensures all code is statically linked on a Linux machine

UNAME := $(shell uname)
MACHINE := $(shell uname -m)

ifeq ($(UNAME), Linux)

        LINUX=1

        ifeq ($(MACHINE), x86_64)
                LINUX64=1
        else
                LINUX32=1
        endif
endif


ifdef LINUX
        CXXFLAGS += -DLINUX
endif

include Rules.mk
