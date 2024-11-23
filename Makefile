.PHONY : extraction

EXTRACTION_MARKDOWN = $(wildcard data/*/*-extraction.md)
# EXTRACTION_IPYNB = ${EXTRACTION_MARKDOWN:.md=.tmp.ipynb}
EXTRACTION_HTML = ${EXTRACTION_MARKDOWN:.md=.html}

extraction : ${EXTRACTION_HTML}

# Run the conversion in two steps: markdown to ipynb, ipynb to html with execution. We
# manually change the working directory using `cd` to avoid having to write a temporary
# ipynb file.
${EXTRACTION_HTML} : %.html : %.md
	cd $(dir $<) \
	&& jupytext --to ipynb --output - $(notdir $<) \
	| jupyter nbconvert --stdin --execute --to html --output $(notdir $@)
