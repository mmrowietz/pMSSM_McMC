# --- variables defined by configure ---

SRC = ./src
QUADSUFFIX = 
PREFIX = /nfs/dust/cms/user/mrowietm/python_scan/packages/FeynHiggs-2.16.1
LIBDIRSUFFIX = 64

EXE = 
DEF = -D
NOUNDERSCORE = 0
QUAD = 0
QUADSIZE = 16
GM2FULL = 0

FC = gfortran
FFLAGS = -O3 -ffixed-line-length-none \
  $(DEF)QUAD=$(QUAD) $(DEF)QUADSIZE=$(QUADSIZE) \
  $(DEF)U77EXT=0 \
  $(DEF)VT100=1 \
  $(DEF)SLHAPARA=0 \
  $(DEF)GM2FULL=$(GM2FULL)

CC = clang
CFLAGS = -O3 -g -fomit-frame-pointer -ffast-math -Wall -m64 \
  -DQUAD=$(QUAD) -DQUADSIZE=$(QUADSIZE) \
  -DNOUNDERSCORE=$(NOUNDERSCORE) \
  -DBIGENDIAN=0 \
  -DVT100=1

CXX = clang++
CXXFLAGS = -O3 -g -fomit-frame-pointer -ffast-math -Wall -m64 -stdlib=libstdc++

ML = 0
MCC = mcc
MCFLAGS = -st -b64
MCLIBS = -lpthread -lrt

LDFLAGS = -L/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/gcc/8.2.0/bin/../libexec/gcc/x86_64-unknown-linux-gnu/8.2.1/liblto_plugin.so -L/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/gcc/8.2.0/bin/../lib/gcc/x86_64-unknown-linux-gnu/8.2.1 -L/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/gcc/8.2.0/bin/../lib/gcc -L/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/gcc/8.2.0/bin/../lib/gcc/x86_64-unknown-linux-gnu/8.2.1/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/gcc/8.2.0/bin/../lib/gcc/x86_64-unknown-linux-gnu/8.2.1/../../../../x86_64-unknown-linux-gnu/lib -L/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/gcc/8.2.0/bin/../lib/gcc/x86_64-unknown-linux-gnu/8.2.1/../../.. -lgfortran -lm -lgcc_s -lgcc -lquadmath -lm -lgcc_s -lgcc -lgcc_s -lgcc

PARALLEL = -j 12

# --- end defs by configure ---

FHVERSION = 2.16.1

BLD = build$(QUADSUFFIX)
VERSION = $(BLD)/version.h

LIBDIR = $(PREFIX)/lib$(LIBDIRSUFFIX)
INCLUDEDIR = $(PREFIX)/include
BINDIR = $(PREFIX)/bin

ARGS = $(PARALLEL) \
  GM2FULL="$(GM2FULL)" \
  NOUNDERSCORE="$(NOUNDERSCORE)" \
  EXE="$(EXE)" \
  FC="$(FC)" \
  FFLAGS="$(FFLAGS)" \
  CC="$(CC)" \
  CFLAGS="$(CFLAGS)" \
  CXX="$(CXX)" \
  CXXFLAGS="$(CXXFLAGS)" \
  QUADSUFFIX="$(QUADSUFFIX)" \
  QUADSIZE="$(QUADSIZE)" \
  ML="$(ML)" \
  MCC="$(MCC)" \
  MCFLAGS="$(MCFLAGS)" \
  LDFLAGS="$(LDFLAGS)" \
  AR="$(AR)" \
  RANLIB="$(RANLIB)" \
  NM="$(NM)" \
  DLLTOOL="$(DLLTOOL)"


default all lib frontend mma: force
	cd $(BLD) && $(MAKE) $(ARGS) $@

FHBIN = $(BLD)/FeynHiggs$(EXE) $(BLD)/table$(EXE)
MFHBIN = $(BLD)/MFeynHiggs$(EXE)

install: default
	-mkdir "$(PREFIX)"
	-mkdir "$(LIBDIR)" "$(INCLUDEDIR)" "$(BINDIR)"
	cp -p $(BLD)/libFH.a "$(LIBDIR)"
	cp -p $(BLD)/CFeynHiggs.h $(BLD)/CSLHA.h $(BLD)/FHCouplings.h \
	  $(BLD)/FHRecord.h $(BLD)/SLHA.h $(BLD)/SLHADefs.h \
	  $(BLD)/PDG.h "$(INCLUDEDIR)"
	strip $(FHBIN)
	cp -p $(FHBIN) $(BLD)/fcc "$(BINDIR)"
	test ! -f $(MFHBIN) || { \
	  test -z "$(EXE)" || ldd $(MFHBIN) | awk '/\/usr\/bin\/cyg/ { system("cp -p " $$3 " $(BINDIR)/"); }' ; \
	  strip $(MFHBIN) ; \
	  cp -p $(MFHBIN) $(BINDIR); \
	}

force: $(BLD)/version.h

$(BLD)/version.h:
	-mkdir $(BLD)
	find "`cd $(SRC) ; pwd`" -follow -exec ln -sf {} $(BLD) \;
	( echo "#define FHVERSION \"$(FHVERSION)\"" ; \
	  echo "#define FHBUILDDATE \"`date '+%b %d, %Y'`\"" ) > $(BLD)/version.h


FH = FeynHiggs-$(FHVERSION)
FHTAR = $(FH).tar.gz

dist:
	rm -f $(FHTAR)
	ln -s . $(FH)
	tar cvfz $(FHTAR) --owner=root --group=root \
	  $(FH)/README $(FH)/configure $(FH)/makefile.in \
	  $(FH)/src $(FH)/extse $(FH)/gen \
	  $(FH)/man $(FH)/example $(FH)/web \
	  $(FH)/COPYING
	rm $(FH)

clean:
	rm -fr $(BLD) $(FHTAR)

.PHONY: default all install force clean dist
