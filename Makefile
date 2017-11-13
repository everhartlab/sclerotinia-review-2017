FILES := results/kamvar2017population.md

all : results/bootstrap.txt $(FILES)

results/bootstrap.txt : code/boot.R
	R --slave -e "options(repos = 'https://cran.rstudio.com'); source('$<', echo = TRUE)" &> $@

results/%.md : code/%.Rmd
	R --slave -e "rmarkdown::render(input = '$<',                              \
	                                output_dir = '$(@D)',                      \
	                                output_format = rmarkdown::md_document(),  \
	                                knit_root_dir = here::here(),              \
	                                )"
results/%.pdf : code/%.Rmd
	R --slave -e "rmarkdown::render(input = '$<',                              \
	                                output_dir = '$(@D)',                      \
	                                output_format = rmarkdown::pdf_document(), \
	                                knit_root_dir = here::here(),              \
	                                )"
results/%.html : code/%.Rmd
	R --slave -e "rmarkdown::render(input = '$<',                               \
	                                output_dir = '$(@D)',                       \
	                                output_format = rmarkdown::html_document(), \
	                                knit_root_dir = here::here(),               \
	                                )"
results/%.docx : code/%.Rmd
	R --slave -e "rmarkdown::render(input = '$<',                               \
	                                output_dir = '$(@D)',                       \
	                                output_format = rmarkdown::word_document(), \
	                                knit_root_dir = here::here(),               \
	                                )"

.PHONY : clean

clean :
	$(RM) $(FILES)
	$(RM) results/bootstrap.txt
