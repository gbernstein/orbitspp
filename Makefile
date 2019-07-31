# These site-dependent items should be defined in environment:

# CXX = g++-6 -fopenmp
# CXXFLAGS

# GBUTILS_DIR
# ASTROMETRY_DIR
# SPICE_DIR

# MKL_DIR (optional, used with Eigen)

# orbits code is NOT being written to support TMV linear algebra.
# So make sure not to use it.
unexport TMV_DIR

INCLUDES := 

LIBS := -lm

EXTDIRS := 

# Collect the includes and libraries we need
ifdef SPICE_DIR
INCLUDES += -I $(SPICE_DIR)/include
LIBS += -L $(SPICE_DIR)/lib -lcspice
else
$(error Require SPICE_DIR in environment)
endif

ifdef GBUTIL_DIR
INCLUDES += -I $(GBUTIL_DIR)/include
EXTDIRS += $(GBUTIL_DIR)
GBUTIL_OBJ = $(GBUTIL_DIR)/obj
else
$(error Require GBUTIL_DIR in environment)
endif

ifdef ASTROMETRY_DIR
INCLUDES += -I $(ASTROMETRY_DIR)/include
EXTDIRS += $(ASTROMETRY_DIR)
ASTROMETRY_OBJ := $(ASTROMETRY_DIR)/obj
else
$(error Require ASTROMETRY_DIR in environment)
endif

ifdef EIGEN_DIR
INCLUDES += -I $(EIGEN_DIR) -D USE_EIGEN
else
$(error Require EIGEN_DIR in environment)
endif

ifdef MKL_DIR
INCLUDES += -I $(MKL_DIR)/include -D USE_MKL
endif

# Not really using this but Astrometry might, so
# keep it for linking if we have it.
ifdef YAML_DIR
LIBS += -L $(YAML_DIR)/lib -lyaml-cpp
endif

ifdef CFITSIO_DIR
INCLUDES += -I $(CFITSIO_DIR)/include
LIBS += -L $(CFITSIO_DIR)/lib -lcfitsio -lcurl
else
$(error Require CFITSIO_DIR in environment)
endif

ifneq ($(origin CFITSIO_NEEDS_CURL),undefined)
LIBS += -lcurl
endif

ifdef GBFITS_DIR
INCLUDES += -I $(GBFITS_DIR)
EXTDIRS += $(GBFITS_DIR)
else
$(error Require GBFITS_DIR in environment)
endif

# Object files found in external packages:
EXTOBJS =$(GBUTIL_OBJ)/StringStuff.o $(GBUTIL_OBJ)/Table.o $(GBUTIL_OBJ)/Pset.o \
	$(GBUTIL_OBJ)/Expressions.o \
	$(ASTROMETRY_OBJ)/Astrometry.o \
	$(GBFITS_DIR)/FITS.o $(GBFITS_DIR)/Header.o $(GBFITS_DIR)/Hdu.o $(GBFITS_DIR)/FitsTable.o \
	$(GBFITS_DIR)/FTable.o $(GBFITS_DIR)/FTableExpression.o \
	$(GBFITS_DIR)/Image.o $(GBFITS_DIR)/FitsImage.o

##### 
BINDIR = bin
OBJDIR = obj
SRCDIR = src
SUBDIR = src/subs
INCLUDEDIR = include
TESTDIR = tests
TESTBINDIR = testbin

# INCLUDES can be relative paths, and will not be exported to subdirectory makes.
INCLUDES += -I $(INCLUDEDIR)

# Executable C++ programs
EXECS :=  $(wildcard $(SRCDIR)/*.cpp)
TARGETS := $(EXECS:$(SRCDIR)/%.cpp=$(BINDIR)/%)
OBJS := $(EXECS:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

# Python executables
PYEXECS :=  $(wildcard $(SRCDIR)/*.py)
PYTARGETS :=  $(PYEXECS:$(SRCDIR)/%.py=$(BINDIR)/%.py)
# C++ subroutines
SUBS :=  $(wildcard $(SUBDIR)/*.cpp)
SUBOBJS := $(SUBS:$(SUBDIR)/%.cpp=$(OBJDIR)/%.o)

CP = /bin/cp -p
RM = /bin/rm -f

#######################
# Rules - ?? dependencies on INCLUDES ??
#######################

all: cpp python

cpp: exts $(TARGETS)

python: $(PYTARGETS)

# No setup.py to do here ... python ./setup.py install

# Compilation
$(OBJS):  $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(SUBOBJS): $(OBJDIR)/%.o : $(SUBDIR)/%.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Linking
$(TARGETS): $(BINDIR)/% : $(OBJDIR)/%.o $(SUBOBJS) $(EXTOBJS)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@

# Python executables - copy into bin directory
$(PYTARGETS): $(BINDIR)/% : $(SRCDIR)/%
	$(CP) $^ $@

# External objects - call external makefiles.  ?? Conditional??
# Semicolon prevents calling default rule for .o/.cpp
$(EXTOBJS): exts ;

######### Test programs

TESTSRC := $(wildcard $(TESTDIR)/*.cpp)
TESTINCLUDE := -I $(TESTDIR)
TESTOBJS := $(TESTSRC:$(TESTDIR)/%.cpp=$(OBJDIR)/%.o)
TESTTARGETS := $(TESTSRC:$(TESTDIR)/%.cpp=$(TESTBINDIR)/%)
TESTSPY := $(wildcard $(TESTDIR)/*.py)

tests: $(TESTTARGETS)

$(TESTOBJS):  $(OBJDIR)/%.o : $(TESTDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(TESTINCLUDE) -c $^ -o $@

$(TESTTARGETS): $(TESTBINDIR)/% : $(OBJDIR)/%.o $(SUBOBJS) $(EXTOBJS)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@

###############################################################
## Standard stuff:
###############################################################

exts:
	for dir in $(EXTDIRS); do (cd $$dir && $(MAKE)); done

depend: local-depend
	for dir in $(EXTDIRS); do (cd $$dir && $(MAKE) depend); done

local-depend:
	$(RM) .depend
	for src in $(SUBS:$(SUBDIR)/%.cpp=%); \
	 do $(CXX) $(CXXFLAGS) $(INCLUDES) -MM $(SUBDIR)/$$src.cpp -MT obj/$$src.o >> .depend; \
	done
	for src in $(EXECS:%.cpp=%); \
	 do $(CXX) $(CXXFLAGS) $(INCLUDES) -MM $$src.cpp -MT obj/$$src.o >> .depend; \
        done

clean: local-clean
	for dir in $(EXTDIRS); do (cd $$dir && $(MAKE) clean); done

local-clean:
	rm -rf $(OBJDIR)/*.o $(BINDIR)/* $(TESTBINDIR)/* *~ *.dvi *.aux core .depend

ifeq (.depend, $(wildcard .depend))
include .depend
endif

.PHONY: all install dist depend clean local-clean local-depend exts tests cpp python

