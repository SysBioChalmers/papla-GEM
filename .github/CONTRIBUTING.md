## Contributor guidelines

Thank you for contributing to papla-GEM! Anybody is welcome to contribute, but please abide by the following guidelines.

You can contribute in two main ways: by creating issues, and by sending pull requests (PRs) with additions, deletions, or corrections to the model.

### Reporting issues in the model

Report an issue at https://github.com/SysBioChalmers/papla-GEM/issues if you note any of the following:

* Incorrect annotation for any model component.
* Missing feature or field you would like the model to have.
* Bug/weird simulation results.
* Lacking documentation.
* Any type of feedback.

When creating the issue, please make sure:

* You tested your code (if any) with all requirements for running the model.
* You did your analysis in the `main` branch of the repository.
* You provide any necessary files/links needed for understanding the issue.
* You checked that a similar issue does not exist already.

Please comply with our [code of conduct](https://github.com/SysBioChalmers/papla-GEM/blob/main/.github/CODE_OF_CONDUCT.md) when commenting on issues.

### Contributing to the model

Here is how to set up papla-GEM for local development:

1. Make sure you have [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN) installed and working in MATLAB.

2. Fork the papla-GEM repository on GitHub.

3. Clone your fork locally:
    ```
    git clone https://github.com/<your GitHub name>/papla-GEM.git
    ```

4. Check out the `develop` branch, and create a branch for local development from it:
    ```
    git checkout develop
    git checkout -b name-of-your-branch
    ```

5. Make your changes to the model in MATLAB, loading it with `importModel('model/papla-GEM.xml')`.

6. From the `code` directory, run `newCommit(model)` to export the updated model files (`.txt`, `.xml`, `.yml`) before committing. Binary formats (`.mat`, `.xlsx`) are only exported on `main`, as part of a release (see `newRelease.m`).

7. Commit your changes and push your branch to GitHub:
    ```
    git add .
    git commit -m "Title of your commit"
    git push origin name-of-your-branch
    ```

8. Submit a pull request to the `develop` branch of the original SysBioChalmers repository (not your fork).

Thank you very much for contributing to papla-GEM!

#### Branching model

* `develop`: The branch all pull requests should be based on, and target.
* `main`: Only touched by the administrator; contains the tested & reviewed model that is released or ready for the next release.

## Administrator guidelines

The main duties of the administrator are:
* To make sure conventions and standards in the model are kept.
* To keep the repository clean and organized.
* To review and merge pull requests into `develop`.
* To generate new releases of the model on `main` using `newRelease.m`, updating `version.txt` and `history.md` accordingly.
