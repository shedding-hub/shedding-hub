.PHONY : backup_data assert_data_unchanged extraction catalog parameters review review_range catalog_ct review_ct review_ct_range catalog_ct_gate2 review_ct_gate2 review_ct_gate2_range

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

# Fit the cycle-threshold analytes, which the shipped catalog excludes, to a
# catalog of their own. Their peak heights are cycles below CT_REFERENCE rather
# than log10 concentrations, so they are kept in a separate file.
catalog_ct :
	python scripts/build_shedding_catalog.py --value-types ct --output shedding_catalog_ct.yaml

# Review pages for those fits. Run after `make catalog_ct`. The y axis carries
# real Ct numbers, so peaks read high on the page and low in cycles.
review_ct :
	python scripts/build_catalog_review.py --catalog shedding_catalog_ct.yaml --output shedding_catalog_review_ct.pdf

review_ct_range :
	python scripts/build_catalog_review.py --catalog shedding_catalog_ct.yaml --range --output shedding_catalog_review_ct_range.pdf

# The Ct catalog under a 2-unit over-extrapolation gate. Note that 2 here means
# 2 *cycles*, not the 2 log10 of the concentration gate2 build: at a slope near
# 3.5 those are about 0.57 log10 against 2, so this is by far the stricter of
# the two. That asymmetry is the point of building it -- it is the same gate
# scale-dependence documented in modeling-methods.md, made visible.
catalog_ct_gate2 :
	python scripts/build_shedding_catalog.py --value-types ct --max-peak-above-observed 2 --output shedding_catalog_ct_gate2.yaml

review_ct_gate2 :
	python scripts/build_catalog_review.py --catalog shedding_catalog_ct_gate2.yaml --output shedding_catalog_review_ct_gate2.pdf

review_ct_gate2_range :
	python scripts/build_catalog_review.py --catalog shedding_catalog_ct_gate2.yaml --range --output shedding_catalog_review_ct_gate2_range.pdf
