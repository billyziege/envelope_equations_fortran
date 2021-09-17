#Compiler and Linker
FC           := gfortran

#The Target Binary Program
TARGET1      := spheroidal_model
TARGET2      := long_chirp_scan
TARGET3      := spheroidal_model_adaptive
TARGETS      := $(TARGET2) $(TARGET3) 

#The Directories, Source, Includes, Objects, Binary and Resources
SRCDIR      := src
INCDIR      := inc
BUILDDIR    := obj
TARGETDIR   := bin
SRCEXT      := f90
OBJEXT      := o
MODEXT      := mod

#Flags, Libraries and Includes
FLAGS       := -O3
MODFLAGS    := -I$(BUILDDIR) -J$(BUILDDIR)

SOURCES        := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS        := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))
FOREIGNOBJS    := $(shell find $(INCDIR) -type f -name *.$(OBJEXT))
FOREIGNMODS    := $(shell find $(INCDIR) -type f -name *.$(MODEXT))
SPHENVDEP      := $(BUILDDIR)/myarray.o $(BUILDDIR)/constants.o
SPHMODELDEP    := $(BUILDDIR)/spheroidal_envelope.o $(BUILDDIR)/myarray.o $(BUILDDIR)/constants.o

#Defauilt Make
all: directories $(TARGETS)

#Remake
remake: cleaner all

#Make the Directories
directories:
	@mkdir -p $(TARGETDIR)
	@mkdir -p $(BUILDDIR)

#Clean only Objecst
clean:
	@$(RM) -rf $(BUILDDIR)

#Full Clean, Objects and Binaries
cleaner: clean
	@$(RM) -rf $(TARGETDIR)

#Link
$(TARGET1): $(BUILDDIR)/$(TARGET1)_main.o $(SPHMODELDEP)
	$(FC) -o $(TARGETDIR)/$(TARGET1) $^ $(FOREIGNOBJS) $(MODFLAGS)
$(TARGET2): $(BUILDDIR)/$(TARGET2)_main.o $(SPHMODELDEP)
	$(FC) -o $(TARGETDIR)/$(TARGET2) $^ $(FOREIGNOBJS) $(MODFLAGS)
$(TARGET3): $(BUILDDIR)/$(TARGET3)_main.o $(SPHMODELDEP)
	$(FC) -o $(TARGETDIR)/$(TARGET3) $^ $(FOREIGNOBJS) $(MODFLAGS)

#Compile
$(BUILDDIR)/$(TARGET1)_main.o: $(SRCDIR)/$(TARGET1)_main.f90 $(SPHMODELDEP)
	@mkdir -p $(dir $@)
	@cp -f $(FOREIGNMODS) $(BUILDDIR)
	$(FC) $(FLAGS) -c -o $@ $< $(MODFLAGS)

$(BUILDDIR)/$(TARGET2)_main.o: $(SRCDIR)/$(TARGET2)_main.f90 $(SPHMODELDEP)
	@mkdir -p $(dir $@)
	@cp -f $(FOREIGNMODS) $(BUILDDIR)
	$(FC) $(FLAGS) -c -o $@ $< $(MODFLAGS)

$(BUILDDIR)/$(TARGET3)_main.o: $(SRCDIR)/$(TARGET3)_main.f90 $(SPHMODELDEP)
	@mkdir -p $(dir $@)
	@cp -f $(FOREIGNMODS) $(BUILDDIR)
	$(FC) $(FLAGS) -c -o $@ $< $(MODFLAGS)

$(BUILDDIR)/spheroidal_envelope.o: $(SRCDIR)/spheroidal_envelope.f90 $(SPHENVDEP)
	@mkdir -p $(dir $@)
	$(FC) $(FLAGS) -c -o $@ $< $(MODFLAGS)

$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(FC) $(FLAGS) -c -o $@ $< $(MODFLAGS)

print-%  : ; @echo $* = $($*)
#Non-File Targets
.PHONY: all remake clean cleaner print-
