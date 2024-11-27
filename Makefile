BUILD_DIR = build

all: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)
	./$(BUILD_DIR)/FlowCommutators

$(BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ ..

clean:
	rm -rf $(BUILD_DIR)