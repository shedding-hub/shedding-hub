.PHONY : backup_data assert_data_unchanged extraction

DATA_FILES = $(wildcard data/*/*.yaml)
DATA_BACKUPS = $(addprefix ${TMPDIR},$(notdir ${DATA_FILES}))
DATA_CHECKS = ${DATA_BACKUPS:.yaml=.null}
EXTRACTION_MARKDOWN = $(wildcard data/*/*-extraction.md)
EXTRACTION_HTML = ${EXTRACTION_MARKDOWN:.md=.html}
EXTRACTION_PY = $(wildcard data/*/*-extraction.py)
EXTRACTION_LOGS = ${EXTRACTION_PY:.py=.log}

extraction : ${EXTRACTION_HTML} ${EXTRACTION_LOGS}

# Run the conversion in two steps: markdown to ipynb, ipynb to html with execution. We
# manually change the working directory using `cd` to avoid having to write a temporary
# ipynb file.
${EXTRACTION_HTML} : %.html : %.md
	cd $(dir $<) \
	&& jupytext --to ipynb --output - $(notdir $<) \
	| jupyter nbconvert --stdin --execute --to html --output $(notdir $@)

# For Python extraction scripts, simply run them in one step and log to output file.
${EXTRACTION_LOGS} : %.log : %.py
	cd $(dir $@) && python $(notdir $<) > $(notdir $@)

backup_data : ${DATA_BACKUPS}

${DATA_BACKUPS} : ${TMPDIR}%.yaml :
	mv data/$*/$*.yaml $@

assert_data_unchanged : ${DATA_CHECKS}

${DATA_CHECKS} : ${TMPDIR}%.null : ${TMPDIR}%.yaml
	python .github/workflows/compare.py data/$*/$*.yaml $<
