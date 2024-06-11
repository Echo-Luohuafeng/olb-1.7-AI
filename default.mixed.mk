include $(OLB_ROOT)/rules.mk

export LC_MESSAGES=C

all: dependencies core $(EXAMPLE)

INCLUDE_DIRS := $(OLB_ROOT)/src

ifeq ($(USE_EMBEDDED_DEPENDENCIES), ON)
INCLUDE_DIRS += \
	$(OLB_ROOT)/external/zlib \
	$(OLB_ROOT)/external/tinyxml

LDFLAGS := -L$(OLB_ROOT)/external/lib $(LDFLAGS)

dependencies:
	$(MAKE) -C $(OLB_ROOT)/external
else
.PHONY: dependencies
endif

LDFLAGS += -L$(OLB_ROOT)/build/lib

core:
	$(MAKE) -C $(OLB_ROOT) core

CPP_FILES := $(wildcard *.cpp)
OBJ_FILES := $(CPP_FILES:.cpp=.o)
DEP_FILES := $(CPP_FILES:.cpp=.d)

INCLUDE_FLAGS := -I$(subst $(EMPTY) $(EMPTY), -I,$(INCLUDE_DIRS))

%.d: $(CPP_FILES)
	@$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $< -MM -MT $(@:.d=.o) >$@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c -o $@ $<

build/missing.txt: $(OBJ_FILES)
	mkdir -p build
	$(CXX) $^ $(LDFLAGS) -lolbcore 2>&1 | grep -oP ".*undefined reference to \`\K[^']+\)" | sort | uniq > $@

EXPLICIT_METHOD_INSTANTIATION := \
	FieldTypeRegistry \
	StatisticsPostProcessor \
	checkPlatform

EXPLICIT_CLASS_INSTANTIATION := \
	Column \
	CyclicColumn \
	ConcreteBlockO \
	ConcreteBlockCollisionO \
	ConcreteBlockCouplingO \
	ConcreteCommunicatable \
	ConcreteBlockCommunicator \
	HeterogeneousCopyTask \
	ConcreteParametersD \
	FieldTypeRegistry \
	unique_ptr \
	StatisticsPostProcessor

OLB_MIXED_MODE_INCLUDE := \#include <olb.h>\n
ifdef OLB_MIXED_MODE_INCLUDE_CPP
	OLB_MIXED_MODE_INCLUDE := $(CPP_FILES:%=\n#include "../%")
endif

build/olbcuda.cu: build/missing.txt
	printf '$(OLB_MIXED_MODE_INCLUDE)' > $@
# Transform missing symbols into explicit template instantiations by:
# - filtering for a set of known and automatically instantiable methods
# - excluding destructors
# - dropping resulting empty lines
# - adding the explicit instantiation prefix (all supported methods are void, luckily)
	cat build/missing.txt \
	| grep '$(subst $() $(),\|,$(EXPLICIT_METHOD_INSTANTIATION))' \
	| grep -wv '.*\~.*\|FieldTypeRegistry()' \
	| xargs -0 -n1 | grep . \
	| sed -e 's/.*/template void &;/' -e 's/void void/void/' >> $@
# - filtering for a set of known and automatically instantiable classes
# - dropping method cruft and wrapping into explicit class instantiation
# - removing duplicates
	cat build/missing.txt \
	| grep '.*\($(subst $() $(),\|,$(EXPLICIT_CLASS_INSTANTIATION))\)<' \
	| sed -e 's/\.*>::.*/>/' -e 's/.*/template class &;/' -e 's/class void/class/' \
	| sort | uniq >> $@

CUDA_CPP_FILES := build/olbcuda.cu
CUDA_OBJ_FILES := $(CUDA_CPP_FILES:.cu=.o)

%.o: %.cu
	$(CUDA_CXX) $(CUDA_CXXFLAGS) $(INCLUDE_FLAGS) $(PARALLEL_FLAGS) -DPLATFORM_CPU_SISD -DPLATFORM_GPU_CUDA -Xcompiler -fPIC -c -o $@ $<

build/olbcuda.version: $(CUDA_OBJ_FILES)
	echo 'LIBOLBCUDA { global: ' > $@
# Declare exposed explicitly instantiated symbols to prevent duplicate definitions by:
# - filtering for the set of automatically instantiated classes
# - excluding CPU_SISD symbols (we only instantiate GPU_CUDA-related symbols)
# - dropping the shared library location information
# - postfixing by semicolons
	nm $(CUDA_OBJ_FILES) \
	| grep '$(subst $() $(),\|,$(EXPLICIT_CLASS_INSTANTIATION))\|cuda.*device\|checkPlatform' \
	| grep -wv '.*sisd.*' \
	| cut -c 20- \
	| sed 's/$$/;/' >> $@
	echo 'local: *; };' >> $@

libolbcuda.so: $(CUDA_OBJ_FILES) build/olbcuda.version
	$(CUDA_CXX) $(CUDA_CXXFLAGS) -Xlinker --version-script=build/olbcuda.version -shared $(CUDA_OBJ_FILES) -o $@

$(EXAMPLE): $(OBJ_FILES) libolbcuda.so
	$(CXX) $(OBJ_FILES) -o $@ $(LDFLAGS) -L . -lolbcuda -lolbcore $(CUDA_LDFLAGS)

$(EXAMPLE)-no-cuda-recompile: $(OBJ_FILES)
	$(CXX) $^ -o $(EXAMPLE) $(LDFLAGS) -L . -lolbcuda -lolbcore $(CUDA_LDFLAGS)

.PHONY: onlysample
onlysample: $(EXAMPLE)

.PHONY: no-cuda-recompile
no-cuda-recompile: $(EXAMPLE)-no-cuda-recompile

.PHONY: clean-tmp
clean-tmp:
	rm -f tmp/*.* tmp/vtkData/*.* tmp/vtkData/data/*.* tmp/imageData/*.* tmp/imageData/data/*.* tmp/gnuplotData/*.* tmp/gnuplotData/data/*.*

.PHONE: clean-core
clean-core:
	make -C $(OLB_ROOT) clean-core

.PHONY: clean
clean: clean-tmp
	rm -f $(OBJ_FILES) $(DEP_FILES) $(EXAMPLE) $(CUDA_OBJ_FILES) $(CUDA_CPP_FILES) libolbcuda.so

-include $(DEP_FILES)
