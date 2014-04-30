check: PSI_Plotter.R
	@echo "OK"

build: preprocess_sample_colors.pdf

preprocess_sample_colors.pdf: preprocess_sample_colors.Rd
	R CMD Rd2pdf --force $^

test: PSI_Plotter.R test_data/INCLUSION_LEVELS-ALL3m-Mmu89-SELECTED.test.tab test_data/Tissues.Mmu.txt
	./$^

PSI_Plotter.R: preprocess_sample_colors.R
