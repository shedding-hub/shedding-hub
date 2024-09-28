# Shedding Hub [![Shedding Hub](https://github.com/shedding-hub/shedding-hub/actions/workflows/build.yaml/badge.svg)](https://github.com/shedding-hub/shedding-hub/actions/workflows/build.yaml)

The Shedding Hub collates data and statistical models for biomarker shedding (such as viral RNA or drug metabolites) in different human specimen (such as stool or sputum samples). Developing wastewater-based epidemiology into a quantitative, reliable epidemiological monitoring tool motivates the project.

Datasets are extracted from appendices, figures, and supplementary materials of peer-reviewed studies. Each dataset is stored as a [`.yaml`](https://en.wikipedia.org/wiki/YAML) file and validated against our [data schema](data/.schema.yaml) to verify its integrity.

## 📊 Getting the Data

You can obtain the data by [downloading it from GitHub](https://github.com/shedding-hub/shedding-hub/tree/main/data). We also provide a [convenient Python package](http://pypi.org/project/shedding-hub/) so you can download the most recent data directly in your code or obtain a specific version of the data for reproducible analysis. Install the package by running `pip install shedding-hub` from the command line. The example below downloads the [data from Wölfel et al. (2020)](https://shedding-hub.github.io/datasets/woelfel2020virological.html) as of the commit [`259ca0d`](https://github.com/shedding-hub/shedding-hub/commit/259ca0d).

```python
>>> import shedding_hub as sh

>>> sh.load_dataset('woelfel2020virological', ref='259ca0d')
{'title': 'Virological assessment of hospitalized patients with COVID-2019',
 'doi': '10.1038/s41586-020-2196-x',
 ...}

```

## 🤝 Contributing

Thank you for contributing your data to the Shedding Hub and supporting wastewater-based epidemiology! If you hit a bump along the road, [create a new issue](https://github.com/shedding-hub/shedding-hub/issues/new) and we'll sort it out together.

We use [pull requests](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests) to add and update data, allowing for review and quality assurance. Learn more about the general workflow [here](https://docs.github.com/en/get-started/using-github/github-flow). To contribute your data, follow these easy steps (if you're already familiar with pull requests, steps 2 and 3 are for you):

1. Create a [fork](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo) of the Shedding Hub repository by clicking [here](https://github.com/shedding-hub/shedding-hub/fork) and [clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) the fork to your computer. You only have to do this once.
2. Create a new `my_cool_study/my_cool_study.yaml` file in the [`data`](data) directory and populate it with your data. See [here](data/woelfel2020virological/woelfel2020virological.yaml) for a comprehensive example from [Wölfel et al. (2020)](https://www.nature.com/articles/s41586-020-2196-x). A minimal example for studies with a single analyte (e.g., SARS-CoV-2 RNA concentration in stool samples) is available [here](tests/examples/valid_single_analyte.yaml), and a minimal example for studies with multiple analytes (e.g., crAssphage RNA concentration in stool samples and caffeine metabolites in urine) is available [here](tests/examples/valid_multiple_analytes.yaml).
3. Optionally, if you have a recent version of [Python](https://www.python.org) installed, you can validate your data to ensure it has the right structure before contributing it to the Shedding Hub.
    - Run `pip install -r requirements.txt` from the command line to install all the Python packages you need.
    - Run `pytest` from the command line to validate all datasets, including the one you just created.
4. Create a new [branch](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-branches) by running `git checkout -b my_cool_study`. Branches let you isolate changes you are making to the data, e.g., if you're simultaneously working on adding multiple studies–much appreciated! You should create a new branch from the `main` branch for each dataset you contribute; see [here](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow) for more information.
5. Add your changes by running `git add data/my_cool_study/my_cool_study.yaml` and commit them by running `git commit -m "Add data from Someone et al. (20xx)."`. Feel free to pick another commit message if you prefer.
6. Push the dataset to your fork by running `git push origin my_cool_study`. This will send the data to GitHub, and the output of the command will include a line `Create a pull reuqest for 'my_cool_study' on GitHub by visiting: https://github.com/[your-username]/shedding-hub/pull/new/my_cool_study`. Click on the link and follow the next steps to create a new pull request.

Congratulations, you've just created your first pull request to contribute a new dataset! We'll now [review the changes](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/about-pull-request-reviews) you've made to make sure everything looks good. Once any questions have been resolved, we'll [merge your changes](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/incorporating-changes-from-a-pull-request/merging-a-pull-request) into the repository. You've just contributed your first dataset to help make wastewater-based epidemiology a more quantitative public health monitoring tool–thank you!
