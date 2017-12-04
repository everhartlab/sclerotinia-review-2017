FILES := results/kamvar2017population.md

all : results/bootstrap.txt $(FILES)

results/bootstrap.txt : code/boot.R DESCRIPTION code/kop
	R --slave -e "options(repos = 'https://cran.rstudio.com'); source('$<', echo = TRUE)" &> $@

results/%.md : code/%.Rmd
	R --slave -e "ezknitr::ezknit(file = '$<',                                 \
	                              out_dir = '$(@D)',                           \
	                              fig_dir = './figures/$*',                    \
	                              keep_html = FALSE                            \
	                              )"
results/%.html : code/%.Rmd
	R --slave -e "ezknitr::ezknit(file = '$<',                                 \
	                              out_dir = '$(@D)',                           \
	                              fig_dir = './figures/$*'                     \
	                              )"
results/%.pdf : code/%.Rmd
	R --slave -e "rmarkdown::render(input = '$<',                              \
	                                output_dir = '$(@D)',                      \
	                                output_format = rmarkdown::pdf_document(), \
	                                knit_root_dir = '.',                       \
	                                )"
results/%.docx : code/%.Rmd code/template.docx
	R --slave -e "rmarkdown::render(input = '$<',                               \
	                                output_dir = '$(@D)',                       \
	                                output_format = rmarkdown::word_document(reference_docx = here::here('code/template.docx')), \
	                                knit_root_dir = '.',                        \
	                                )"

.PHONY : clean

clean :
	$(RM) $(FILES)
	$(RM) results/bootstrap.txt
