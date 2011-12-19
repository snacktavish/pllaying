TAG = GCC-SS3

#CONFIGURE BUILD SYSTEM
TARGET	   = RAxML-$(TAG)
BUILD_DIR  = ./$(TAG)
SRC_DIR    = ./
MAKE_DIR   = ./system
Q         ?= @


#DO NOT EDIT BELOW
include $(MAKE_DIR)/include_$(TAG).mk

SRC += axml.c bipartitionList.c evaluateGenericSpecial.c evaluatePartialGenericSpecial.c makenewzGenericSpecial.c models.c optimizeModel.c parsePartitions.c restartHashTable.c searchAlgo.c topologies.c trash.c treeIO.c newviewGenericSpecial.c
LIBS += -lm

#OBJ       = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.c))
OBJ     += $(patsubst %.c, $(BUILD_DIR)/%.o,$(SRC))
SOBJ    += $(patsubst %.c, $(BUILD_DIR)/%_special.o,$(SPECIAL_SRC))
CCFLAGS := $(CCFLAGS) $(DEFINES) $(INCLUDES) 

${TARGET}: $(BUILD_DIR) $(OBJ) $(SOBJ)
	@echo "===>  LINKING  $(TARGET)"
	$(Q)${CC} ${LFLAGS} -o $(TARGET) $(OBJ) $(SOBJ) $(LIBS)

$(BUILD_DIR)/%.o:  %.c
	@echo "===>  COMPILE  $@"
	$(Q)$(CC) -c $(CCFLAGS) $< -o $@
	$(Q)$(CC) $(CCFLAGS) -MT $(@:.d=.o) -MM  $< > $(BUILD_DIR)/$*.d

$(BUILD_DIR)/%_special.o:  %.c
	@echo "===>  COMPILE  $@"
	$(Q)$(CC) -c $(SPECIAL_FLAGS) $(CCFLAGS) $< -o $@
	$(Q)$(CC) $(CCFLAGS) -MT $(@:.d=.o) -MM  $< > $(BUILD_DIR)/$*.d

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

