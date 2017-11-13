FILES := results/kamvar2017population.md

all : results/bootstrap.txt $(FILES)

results/bootstrap.txt : code/boot.R
	R --slave -e "options(repos = 'https://cran.rstudio.com'); source('$<', echo = TRUE)" &> $@

results/%.md : code/%.Rmd
	R --slave -e "rmarkdown::render(input = '$<',                 \
	                                output_file = '$@',           \
	                                knit_root_dir = here::here(), \
	                                )"
results/%.pdf : code/%.Rmd
	R --slave -e "rmarkdown::render(input = '$<',                 \
	                                output_file = '$@',           \
	                                knit_root_dir = here::here(), \
	                                )"
results/%.html : code/%.Rmd
	R --slave -e "rmarkdown::render(input = '$<',                 \
	                                output_file = '$@',           \
	                                knit_root_dir = here::here(), \
	                                )"