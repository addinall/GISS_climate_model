# makefile for the component dd2d
.SUFFIXES:

SUPPORTED_OPTIONS = ESMF NC_IO

FSRCS = cdl_mod.f

ifeq ($(NC_IO),PNETCDF)
FSRCS += pario_pnc.f
else
FSRCS += pario_nc.f
endif

ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

include $(CONFIG_DIR)/base.mk
include $(CONFIG_DIR)/rules.mk

