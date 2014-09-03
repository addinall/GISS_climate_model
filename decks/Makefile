.SUFFIXES:
.PHONY: rundeck depend gcm setup clean vclean newstart exe cmp .R htmldoc setup_script setup_nocomp tests check_sanity

# modelE directory structure
MODEL_E_ROOT = ..
MODEL_DIR = $(MODEL_E_ROOT)/model
AUX_DIR = $(MODEL_E_ROOT)/aux
SCRIPTS_DIR = $(MODEL_E_ROOT)/exec
DECKS_DIR = $(shell pwd)
BIN_DIR = $(DECKS_DIR)/$(RUN)_bin
CONFIG_DIR = $(MODEL_E_ROOT)/config
TEMPLATES_DIR = $(MODEL_E_ROOT)/templates
HTMLDOC_DIR = $(DECKS_DIR)/$(RUN)_html
export HTMLDOC_DIR

RUNSRC  = E1M20
RUN     =
RUNDECKSRC = $(TEMPLATES_DIR)/$(RUNSRC).R
RUNDECK = $(DECKS_DIR)/$(RUN).R

USER = $(shell whoami)
HOMEDIR = $(wildcard ~)
MODELERC ?= $(HOMEDIR)/.modelErc
#default number of processors for setup
NPROC = 1
NPES = 1
# getting options from configuration files
# default modelE configuration file:
include $(SCRIPTS_DIR)/modelErc
# configuration file with user-defined options:
sinclude $(MODELERC)
export COMPILER
export MODELERC

V =

ifeq ($(VERBOSE_OUTPUT),NO)
  MAKEFLAGS=-s
endif

# SETUP = $(SCRIPTS_DIR)/setup_e

RUNSRC = E1M20

ALIAS_LIST = CMPE001 qc

ifeq ($(MPP),YES)
  MPI=YES
endif

ifeq ($(ESMF),YES)
  MPI=YES
endif


help:
	$(MAKE) -C $(MODEL_DIR) help

sinclude $(CONFIG_DIR)/rules.mk
sinclude $(DECKS_DIR)/$(RUN).mk

ifort_version:
	@echo $(IFORT_RELEASE)

rundeck:
	@if [ -f $(DECKS_REPOSITORY)/$(RUN).R -a \
	-z "`find $(DECKS_REPOSITORY)/$(RUN).R -user $(USER)`" ]; then \
	echo 'name $(RUN) is already used by somebody else'; \
	echo 'please choose a unique one'; exit 1  ; fi
	@if [ -s $(DECKS_REPOSITORY)/$(RUN).R -a $(OVERWRITE) != YES ];\
	then echo '$(RUN).R already exists in repository'; \
	echo 'use OVERWRITE=YES to overwrite'; \
	exit 1; fi
	#$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) RUNSRC=$(RUNSRC)
	@if [ ! -s $(RUNDECKSRC) ]; then \
	echo 'Source rundeck $(RUNDECKSRC) does not exist'; exit 1; fi
	@if [ "$(RUN)" = "" ]; then \
	echo 'You should specify run name on a command line (i.e. RUN=...)'; \
	exit 1; fi
	@echo '---- Regenerating Rundeck $(RUNDECK) from default ($(RUNDECKSRC)) ----'
	@if [ $(RUNDECK) != $(RUNDECKSRC) ]; then  \
	sed 's/^'$(RUNSRC)'\([:\. ]\)/'$(RUN)'\1/g' $(RUNDECKSRC) \
	| $(SCRIPTS_DIR)/protect_ccp_options.pl | $(CPP) -I$(TEMPLATES_DIR) \
	| $(SCRIPTS_DIR)/protect_ccp_options.pl -u \
	> $(RUNDECK); \
	fi;
	@echo '---------             Rundeck  created               ---------' 
	@cp -f $(RUN).R $(DECKS_REPOSITORY)/$(RUN).R
	@chmod 644 $(DECKS_REPOSITORY)/$(RUN).R      #just in case ...
	@echo '### Information below describes your run. Do not delete! ###'\
	    >> $(DECKS_REPOSITORY)/$(RUN).R

depend: $(RUN).mk
	$(MAKE) check_sanity RUN=$(RUN)
	$(MAKE) -C $(MODEL_DIR) depend_all RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)


gcm $(BIN_DIR)/$(RUN).exe: $(RUN).mk $(BIN_DIR) depend
	$(MAKE) -C $(MODEL_DIR) gcm RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)
	-mv $(MODEL_DIR)/$(RUN).bin $(BIN_DIR)/$(RUN).exe

setup: gcm
	$(MAKE) setup_script RUN=$(RUN)

tests: 
	$(MAKE) -C $(MODEL_DIR) unit_tests

setup_nocomp:
	@[ -d $(BIN_DIR) ] || mkdir $(BIN_DIR)
	@if [ ! -e $(BIN_DIR)/$(RUN).exe ]; then \
	  cp $(CMRUNDIR)/$(RUN)/$(RUN).exe $(BIN_DIR)/$(RUN).exe || \
	  { echo "No executable in run directory"; \
	  echo "You have to run 'gmake setup ...' first"; exit 1; } \
	fi
	$(MAKE) setup_script RUN=$(RUN)

setup_script:
	@echo '---------       Looks like it was compiled OK        ---------'
	@echo '----- Saving Rundeck and other info to global repository -----'
	@if [ ! -d $(DECKS_REPOSITORY) ]; then \
	 echo '--- $(DECKS_REPOSITORY) does not exist. Rundeck not saved ---';\
	else \
	  if [ -f $(DECKS_REPOSITORY)/$(RUN).R -a \
	    -z "`find $(DECKS_REPOSITORY)/$(RUN).R -user $(USER)`" ]; then \
	    echo 'name $(RUN) is already used by somebody else'; \
	    echo 'please choose a unique one'; exit 1; \
	  else \
	    cp -f $(RUN).R $(DECKS_REPOSITORY)/$(RUN).R ;\
	    chmod 644 $(DECKS_REPOSITORY)/$(RUN).R  ;    \
	    echo '### Information below describes your run. Do not delete! ###'\
	    >> $(DECKS_REPOSITORY)/$(RUN).R ; \
	    date >> $(DECKS_REPOSITORY)/$(RUN).R ; \
	    echo $(F90_VERSION) >> $(DECKS_REPOSITORY)/$(RUN).R ;\
	    if [ -f $(MODEL_DIR)/CVS/Tag ]; then \
	      echo "CVS Repository Tag: `cat $(MODEL_DIR)/CVS/Tag`" \
	      >> $(DECKS_REPOSITORY)/$(RUN).R; else \
	      echo 'CVS Repository: MAIN Branch' \
	      >> $(DECKS_REPOSITORY)/$(RUN).R ;\
	    fi ;\
	  fi ;\
	fi
	@echo '---------        Starting setup for $(RUN)        ----------'
	@echo '--------------------------------------------------------------'
ifeq ($(MPI),YES)
	@$(SETUP) $(SETUP_FLAGS) -mpidistr "$(MPIDISTR)"  -mpi $(NPES) $(RUN)
else
	@$(SETUP) $(SETUP_FLAGS) -omp $(NPROC) $(RUN)
endif

$(SETUP_GFDL):
	$(SCRIPTS_DIR)/create_setup_gfdl

setup_gfdl: gcm $(SETUP_GFDL)
	@echo '---------       Looks like it was compiled OK        ---------'
	@echo '----- Saving Rundeck and other info to global repository -----'
	@if [ ! -d $(DECKS_REPOSITORY) ]; then \
	  echo '--- $(DECKS_REPOSITORY) does not exist. Rundeck not saved ---';\
	else \
	  if [ -f $(DECKS_REPOSITORY)/$(RUN).R -a \
	    -z "`find $(DECKS_REPOSITORY)/$(RUN).R -user $(USER)`" ]; then \
	    echo 'name $(RUN) is already used by somebody else'; \
	    echo 'please choose a unique one'; exit 1; \
	  else \
	    cp -f $(RUN).R $(DECKS_REPOSITORY)/$(RUN).R ;\
	    chmod 644 $(DECKS_REPOSITORY)/$(RUN).R  ;    \
	    echo '### Information below describes your run. Do not delete! ###'\
	    >> $(DECKS_REPOSITORY)/$(RUN).R ; \
	    date >> $(DECKS_REPOSITORY)/$(RUN).R ; \
	    echo $(F90_VERSION) >> $(DECKS_REPOSITORY)/$(RUN).R ;\
	    if [ -f $(MODEL_DIR)/CVS/Tag ]; then \
	      echo "CVS Repository Tag: `cat $(MODEL_DIR)/CVS/Tag`" \
	      >> $(DECKS_REPOSITORY)/$(RUN).R; else \
	      echo 'CVS Repository: MAIN Branch' \
	      >> $(DECKS_REPOSITORY)/$(RUN).R ;\
	    fi ;\
	  fi ;\
	fi
	@echo '---------        Starting setup for $(RUN)          ---------'
	@echo '--------------------------------------------------------------'
	@$(SETUP_GFDL) $(RUN)


aux rvr csrvr CMPE002: depend $(BIN_DIR)
	$(MAKE) -C $(MODEL_DIR) gcmlib RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)
	$(MAKE) -C $(AUX_DIR) $@ RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)

auxinput: aux
	$(MAKE) input -C $(AUX_DIR) RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)

auxqflux: aux
	$(MAKE) qflux -C $(AUX_DIR) RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)

auxdeep: aux
	$(MAKE) odeep -C $(AUX_DIR) RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)

clean vclean clean_all:
	$(MAKE) -C $(MODEL_DIR) $@ $(OPTS_MAIN) RUN=$(RUN) DECKS_DIR=$(DECKS_DIR)
	$(MAKE) -C $(AUX_DIR) $@ $(OPTS_MAIN)

exe: gcm
	@if [ -n "`find $(DECKS_REPOSITORY)/$(RUN).R -user $(USER)`" ] ; then \
	echo 'Rundeck in repository is ok - will add info'; \
	else \
	if [ -f $(DECKS_REPOSITORY)/$(RUN).R ] ; then \
	  echo 'repository file for $(RUN) belongs to somebody else'; \
	  echo 'Info not written. Please check!'; exit 1 ; fi ; \
	echo 'repository file for $(RUN) is missing'; \
	echo 'Either you didn''t run setup yet or repository is corrupted'; \
	echo 'PLEASE CHECK! You may be loosing important data!'; echo ; \
	echo 'trying to create a new repository file'; \
	echo '--- no information on original rundeck ---' \
	  > $(DECKS_REPOSITORY)/$(RUN).R ; \
	cat $(RUN).R >> $(DECKS_REPOSITORY)/$(RUN).R ; \
	echo '### Information below describes your run.' \
	  >> $(DECKS_REPOSITORY)/$(RUN).R ; \
	if [ ! -f $(DECKS_REPOSITORY)/$(RUN).R ] ; then \
	  echo 'no luck...' ; exit 1 ; fi ; \
	echo 'new $(DECKS_REPOSITORY)/$(RUN).R created'; fi
	@if [ ! -s $(CMRUNDIR)/$(RUN)/I ] ; then \
	echo \
	  'directory $(CMRUNDIR)/$(RUN) or input files inside it are missing';\
	echo 'you should make "setup" first'; \
	exit 1; fi
	mv -f $(BIN_DIR)/$(RUN).exe $(CMRUNDIR)/$(RUN)/
	@echo 'adding corresponding info to rundeck repository'
	@perl -e "while(<>){ if(/^### Information/){exit;} print; }" \
	<$(DECKS_REPOSITORY)/$(RUN).R  >$(DECKS_REPOSITORY)/$(RUN).R.tmp
	@echo '### Executable was recompiled ' >> $(DECKS_REPOSITORY)/$(RUN).R
	@date >> $(DECKS_REPOSITORY)/$(RUN).R
	@echo $(F90_VERSION) >> $(DECKS_REPOSITORY)/$(RUN).R
	@if [ -f $(MODEL_DIR)/CVS/Tag ]; then \
	echo "CVS Repository Tag: `cat $(MODEL_DIR)/CVS/Tag`" \
	>> $(DECKS_REPOSITORY)/$(RUN).R; else \
	echo 'CVS Repository: MAIN Branch' \
	>> $(DECKS_REPOSITORY)/$(RUN).R; fi
	@echo 'diff old.R new.R :' >> $(DECKS_REPOSITORY)/$(RUN).R
	@diff $(DECKS_REPOSITORY)/$(RUN).R.tmp $(RUN).R \
	>> $(DECKS_REPOSITORY)/$(RUN).R || echo 'rundeck changed'
	@rm -f $(DECKS_REPOSITORY)/$(RUN).R.tmp
	@echo done

newstart:
	@if [ "$(RUN)" = "" ]; then \
	echo 'You should specify run name on a command line (i.e. RUN=...)'; \
	exit 1; fi
	rm -i $(CMRUNDIR)/$(RUN)/*

alias: aux
	for i in $(ALIAS_LIST); do ln -sf $(BIN_DIR)/$$i .; done

unalias:
	for i in $(ALIAS_LIST); do rm -f $$i; done

htmldoc: depend
	-rm -rf $(HTMLDOC_DIR)
	$(MAKE) -C $(MODEL_DIR) htmldoc_all RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)

config configure:
	@if [ -f $(MODELERC) -a $(OVERWRITE) != YES ]; then \
	echo '~/.modelErc already exists, use OVERWRITE=YES to overwrite'; \
	exit 1; fi
	@cp -f $(SCRIPTS_DIR)/modelErc $(MODELERC)
	@if [ `hostname` = *'giss.nasa.gov' ]; then \
	sed "s/MAILTO=/MAILTO=`whoami`@giss.nasa.gov/" < $(MODELERC) > \
	$(MODELERC).tmp; mv -f $(MODELERC).tmp $(MODELERC); fi
	@echo '$(MODELERC) created.'
	@echo 'You may want to edit it manually to customize the options.'

update:
	$(MAKE) -C $(MODEL_E_ROOT) $@ RELEASE=$(RELEASE)

src_cpp_dir:
	$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)

etags:
	$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)

$(RUN).mk: $(RUN).R $(SCRIPTS_DIR)/r_to_mk.pl
	$(SCRIPTS_DIR)/r_to_mk.pl $< > $@

#sinclude $(DECKS_DIR)/$(RUN).mk

FSRCS = $(addsuffix .f,$(strip $(OBJ_LIST)))

echo:
	echo XX$(strip $(OBJ_LIST))XX
	echo XX$(FSRCS)XX

$(BIN_DIR):
	mkdir $(BIN_DIR)

check_sanity:
	@if [ "$(COMPONENTS)"xx = xx ] ; then \
	  echo "Your rundeck contains no Components" ; \
	  echo "Old style rundeck ?  ... aborting ..." ; \
	  exit 1 ; fi
	@if [ $(MPI)x = YESx -a $(ESMF)x != YESx ] ; then \
	  echo "You need ESMF=YES for MPI runs. Aborting" ; \
	  exit 1 ; fi

