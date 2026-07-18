BUILD_DIR = build

TeXOptions = -lualatex \
			 -interaction=nonstopmode \
			 -halt-on-error \
			 -output-directory=$(BUILD_DIR)

all: $(BUILD_DIR)/Makefile $(BUILD_DIR)/main.pdf

output.tex: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)
	./$(BUILD_DIR)/FlowCommutators >output.tex

$(BUILD_DIR)/Makefile: CMakeLists.txt sources/*
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake ..

$(BUILD_DIR)/main.pdf: output.tex main.tex
	latexmk $(TeXOptions) main.tex 1> $(BUILD_DIR)/log || cat $(BUILD_DIR)/log

clean:
	rm -rf $(BUILD_DIR)