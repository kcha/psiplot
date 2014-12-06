install:
	R --no-save -e 'library(devtools);library(roxygen2);document();install()'

test:
	./tests/run_tests.R
