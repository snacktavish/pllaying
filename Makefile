TAG = GCC-SEQ
DIM=3D

#CONFIGURE BUILD SYSTEM
TARGET	   = MAIN$(DIM)-$(TAG)
BUILD_DIR  = ./$(TAG)$(DIM)
SRC_DIR    = ./include
MAKE_DIR   = ./system
Q         ?= @

DEFINES  +=  -DEPOT -DDIFFUSION -DTHERMO

#DO NOT EDIT BELOW
include $(MAKE_DIR)/include_$(TAG).mk

INCLUDES += -I./include/

ifeq ($(DIM),2D)
	DEFINES  += -DTWODIM
endif


VPATH     = $(SRC_DIR)
OBJ       = $(patsubst $(SRC_DIR)/%.C, $(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.C))
CPPFLAGS := $(CPPFLAGS) $(DEFINES) $(INCLUDES) 

${TARGET}: $(BUILD_DIR) $(OBJ)
	@echo "===>  LINKING  $(TARGET)"
	$(Q)${CXX} ${LFLAGS} -o $(TARGET) $(OBJ) $(LIBS)

$(BUILD_DIR)/%.o:  %.C
	@echo "===>  COMPILE  $@"
	$(Q)$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@
	$(Q)$(CXX) $(CPPFLAGS) -MT $(@:.d=.o) -MM  $< > $(BUILD_DIR)/$*.d

$(BUILD_DIR):
	@mkdir $(BUILD_DIR)

ifeq ($(findstring $(MAKECMDGOALS),clean),)
-include $(OBJ:.o=.d)
endif

.PHONY: clean

clean:
	@echo "===>  CLEAN"
	@rm -rf $(BUILD_DIR)
	@rm -f $(TARGET)

