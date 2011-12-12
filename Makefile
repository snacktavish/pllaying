TAG = GCC-SS3

#CONFIGURE BUILD SYSTEM
TARGET	   = RAxML-$(TAG)
BUILD_DIR  = ./$(TAG)
SRC_DIR    = ./
MAKE_DIR   = ./system
Q         ?= @


#DO NOT EDIT BELOW
include $(MAKE_DIR)/include_$(TAG).mk


VPATH     = $(SRC_DIR)
OBJ       = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.c))
CCFLAGS := $(CCFLAGS) $(DEFINES) $(INCLUDES) 

${TARGET}: $(BUILD_DIR) $(OBJ)
	@echo "===>  LINKING  $(TARGET)"
	$(Q)${CC} ${LFLAGS} -o $(TARGET) $(OBJ) $(LIBS)

$(BUILD_DIR)/%.o:  %.c
	@echo "===>  COMPILE  $@"
	$(Q)$(CC) -c $(CCFLAGS) $< -o $@
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

