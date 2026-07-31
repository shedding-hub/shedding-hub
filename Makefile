.PHONY : backup_data assert_data_unchanged extraction catalog parameters review review_range

EXTRACTION_MARKDOWN = $(wildcard data/*/*-extraction.md)
EXTRACTION_HTML = ${EXTRACTION_MARKDOWN:.md=.html}
EXTRACTION_PY = $(wildcard data/*/*-extraction.py)
EXTRACTION_LOGS = ${EXTRACTION_PY:.py=.log}
DATA_FILES_PY = ${EXTRACTION_PY:%-extraction.py=%.yaml}
DATA_FILES_MARKDOWN = ${EXTRACTION_MARKDOWN:%-extraction.md=%.yaml}
DATA_FILES = ${DATA_FILES_PY} ${DATA_FILES_MARKDOWN}
DATA_BACKUPS = $(addprefix ${TMPDIR},$(notdir ${DATA_FILES}))
DATA_CHECKS = ${DATA_BACKUPS:.yaml=.null}

extraction : ${DATA_FILES}

# Run the conversion in two steps: markdown to ipynb, ipynb to html with execution. We
# manually change the working directory using `cd` to avoid having to write a temporary
# ipynb file.
${DATA_FILES_MARKDOWN} : %.yaml : %-extraction.md
	cd $(dir $<) \
	&& jupytext --to ipynb --output - $(notdir $<) \
	| jupyter nbconvert --stdin --execute --to html --output $(notdir $*)-extraction.html

# For Python extraction scripts, simply run them in one step and log to output file.
${DATA_FILES_PY} : %.yaml : %-extraction.py
	cd $(dir $@) && python $(notdir $<) > $(notdir $*)-extraction.log

backup_data : ${DATA_BACKUPS}

${DATA_BACKUPS} : ${TMPDIR}%.yaml :
	mv data/$*/$*.yaml $@

assert_data_unchanged : ${DATA_CHECKS}

${DATA_CHECKS} : ${TMPDIR}%.null : ${TMPDIR}%.yaml
	python .github/workflows/compare.py data/$*/$*.yaml $<

# Refit every analyte in data/ and rewrite the shipped catalog. Slow by design;
# run it whenever datasets are added or changed.
catalog :
	python scripts/build_shedding_catalog.py

# Write the fitted parameters for every dataset, JSON and CSV, from the shipped
# catalog. Fits nothing; run it after `make catalog`.
parameters :
	python scripts/export_parameter_table.py

# Render every catalog fit against the data behind it, one page each, for
# review. The PDF is regenerable and deliberately untracked.
review :
	python scripts/build_catalog_review.py

# The same pages, but shading the full range of the simulated cohort rather than
# its central 90%, with the y axis widened to fit. Shows what each fit considers
# possible rather than typical.
review_range :
	python scripts/build_catalog_review.py --range
