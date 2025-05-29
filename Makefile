.PHONY : backup_data assert_data_unchanged extraction metadata create_metadata

EXTRACTION_MARKDOWN = $(wildcard data/*/*-extraction.md)
EXTRACTION_HTML = ${EXTRACTION_MARKDOWN:.md=.html}
EXTRACTION_PY = $(wildcard data/*/*-extraction.py)
EXTRACTION_LOGS = ${EXTRACTION_PY:.py=.log}
DATA_FILES_PY = ${EXTRACTION_PY:%-extraction.py=%.yaml}
DATA_FILES_MARKDOWN = ${EXTRACTION_MARKDOWN:%-extraction.md=%.yaml}
DATA_FILES = ${DATA_FILES_PY} ${DATA_FILES_MARKDOWN}
DATA_BACKUPS = $(addprefix ${TMPDIR},$(notdir ${DATA_FILES}))
DATA_CHECKS = ${DATA_BACKUPS:.yaml=.null}
YAML_FILES := $(wildcard data/*/*.yaml)
DATASET_NAMES := $(basename $(notdir $(YAML_FILES)))
METADATA_FILES := $(patsubst %,data/%/metadata.jsonld,$(DATASET_NAMES))

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

# Master target to trigger all metadata generation
create_metadata: data/summary.yaml $(METADATA_FILES)

# Single script that generates both metadata.jsonld and summary.yaml
data/summary.yaml $(METADATA_FILES): .github/workflows/metadata.py $(YAML_FILES)
	python .github/workflows/metadata.py