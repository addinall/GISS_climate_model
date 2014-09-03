#makefile for the component Ent
.SUFFIXES:

SUPPORTED_OPTIONS = ESMF

FSRCS = 
F90SRCS = $(wildcard *.F90)
F90OBJS = $(F90SRCS:%.F90=%.o)

ifeq ($(ESMF),YES)
$(F90OBJS) : F90FLAGS+=-DUSE_MPI
endif

include $(CONFIG_DIR)/base.mk
include $(CONFIG_DIR)/rules.mk




