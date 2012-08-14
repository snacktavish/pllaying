TAG = REORDER
 

#CONFIGURE BUILD SYSTEM
TARGET	   = RAxML-$(TAG)
BUILD_DIR  = ./$(TAG)
SRC_DIR    = ./
MAKE_DIR   = ./system
#Q         ?= @


#DO NOT EDIT BELOW
include $(MAKE_DIR)/include_$(TAG).mk

SRC += axml.c bipartitionList.c evaluateGenericSpecial.c evaluatePartialGenericSpecial.c makenewzGenericSpecial.c models.c optimizeModel.c restartHashTable.c searchAlgo.c topologies.c trash.c treeIO.c newviewGenericSpecial.c fastDNAparsimony.c randomTree.c recom.c

CXX_SRC += sim_reorder.cpp

LIBS += -lm

#OBJ       = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.c))
OBJ     += $(patsubst %.c, $(BUILD_DIR)/%.o,$(SRC))
SOBJ    += $(patsubst %.c, $(BUILD_DIR)/%_special.o,$(SPECIAL_SRC))
OBJ +=      $(patsubst %.cpp, $(BUILD_DIR)/%.o,$(CXX_SRC))

CCFLAGS := $(CCFLAGS) $(DEFINES) $(INCLUDES) 

${TARGET}: $(BUILD_DIR) $(OBJ) $(SOBJ)
	@echo "===>  LINKING  $(TARGET)"
	$(Q)${CXX} ${LFLAGS} -o $(TARGET) $(OBJ) $(SOBJ) $(LIBS)

$(BUILD_DIR)/%.o:  %.c
	@echo "===>  COMPILE  $@"
	$(Q)$(CC) -c $(CCFLAGS) $< -o $@
	$(Q)$(CC) $(CCFLAGS) -MT $(@:.d=.o) -MM  $< > $(BUILD_DIR)/$*.d

$(BUILD_DIR)/%_special.o:  %.c
	@echo "===>  COMPILE  $@"
	$(Q)$(CC) -c $(SPECIAL_FLAGS) $(CCFLAGS) $< -o $@
	$(Q)$(CC) $(CCFLAGS) -MT $(@:.d=.o) -MM  $< > $(BUILD_DIR)/$*.d


$(BUILD_DIR)/%.o:  %.cpp
	@echo "===>  COMPILE  $@"
	$(Q)$(CXX) -c $(CCFLAGS) $< -o $@
	$(Q)$(CXX) $(CCFLAGS) -MT $(@:.d=.o) -MM  $< > $(BUILD_DIR)/$*.d

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

