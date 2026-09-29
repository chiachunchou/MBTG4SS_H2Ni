
#------------------------------------

# QTR source code

OBJECTS += Qtr.o Error.o Main.o InputFile.o Job.o Job_Scatter2d.o Log.o Parameters.o RandNum.o Scatter2d.o 

TEMPOBJ := $(OBJECTS)
DEPOBJECTS := $(addprefix ../,$(TEMPOBJ))
DEPLIBS := $(addprefix ../,$(LIBS))

#------------------------------------

# Build rules

all: $(POTDIRS) $(FPOTDIRS) qtr
        #@echo
        #@echo "QTR Compilation"

qtr: $(OBJECTS) $(LIBS)
	$(CXX) -o $(TARGET_NAME) $^ $(LDFLAGS)

$(LIBS):
	$(MAKE) -C $@

$LIBS: $(POTDIRS) $(FPOTDIRS)

$(POTDIRS):
	$(MAKE) -C $@ CC="$(CC)" CXX="$(CXX)" LD="$(LD)" AR="$(AR)" RANLIB="$(RANLIB)" CXXFLAGS="$(CXXFLAGS)"

$(FPOTDIRS):
	$(MAKE) -C $@ CC="$(CC)" CXX="$(CXX)" LD="$(LD)" AR="$(FAR)" FC="$(FC)" FFLAGS="$(FFLAGS)" RANLIB="$(RANLIB)" CXXFLAGS="$(CXXFLAGS)"

clean:
	rm -f $(OBJECTS) $(DEPENDS) $(TARGET_NAME)

clean-all: clean

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $<

DEPENDS= $(wildcard *.d)
-include $(DEPENDS)

.PHONY : all $(POTDIRS) $(FPOTDIRS) clean clean-all version.h
# DO NOT DELETE
