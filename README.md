# Project AAAA_2019_Project-Template

## Project overview

[Write a concise overview of the project]

## Contributors

[List contributors to the project, affiliations, and if appropriate contact details]

## Project setup

### Define project name and identifier

To use this template project for a new project, first copy the directory or clone the repository to a new location.

The convention for project names is that they are identified by a randomly-generated four-letter code (in the format `[A-Z][AEIOU][A-Z][A-Z][A-Z]`), the year in which the project commenced and a project name (title case with words separated by hyphens). These elements are separated by underscores.

Thus, this template project has a valid identifier ("project ID"): `AAAA_2019_Project-Template`.

R code to generate a unique four-letter project code (setting a random seed based on the date it is run):
```
set.seed(as.numeric(Sys.Date()))
c(sample(LETTERS, 1), sample(c("A","E","I","O","U"), 1), sample(LETTERS, 2))
```

### Update template for new project



Specify the way to cite the project in the `CITATION` file. This will likely need updating over the course of the project.

Change the seed for random number generation for the project in `_workflowr.yml` to allow reproducibility of anlayses that use "random" numbers. A sensible choice is the date that the project is setup in YYYYMMDD format (i.e. a numeric value).

Check that the `LICENSE` file is appropriate. Unless there is good reason to prefer something different, we strongly prefer open, permissive licenses to make our work as widely accessible and useable as possible.

The `Snakefile` file that defines the Snakemake workflow will need to be rewritten, but the included file in this template provides a useful starting point and some handy patterns to exploit.

The `analysis/about.Rmd`, `analysis/index.Rmd` and `analysis/license.Rmd` files need to be completed, but the existing templates from a previous project provide a useful guide and starting point.

## Project organisation and management

### File naming conventions

Files should have names that follow this pattern: `<project code>_<date file created>_<descriptive name>.<file extension>`. The four-letter project code should be used (e.g. AAAA). The YYYY-MM-DD format should be used (e.g. 2019-01-01 for January 1 2019). Hyphens should be used to separate words in the descriptive name (e.g. exploratory-data-analysis). 

Thus a valid file name is `AAAA_2019-01-01_exploratory-data-analysis.Rmd`.

Applying these conventions ensures that user-created files are easily identifiable and uniquely named.

### Correspondence about the project

All emails sent regarding the project should have the four-letter project code at the start of the email subject (e.g. `AAAA - <topic of email>`). This ensures that email correspondence for the project is easily searchable in overstuffed email archives and inboxes.

Particularly important emails should be exported to PDF and saved in a subfolder (e.g. `correspondence`) of the `org` folder.

### To-dos and notes

The `org` folder is intended to house files used for organising and managing the project. These might include markdown (`.md`) files for notes, and `.org` mode files (another plaintext document format compatible with Emacs org-mode and org-mode emulators in other modern text editors) for to-do lists and outlines.

There is an included org-mode file `org/project_management.org` that can act as a template for organising the project with a roadmap and to-do lists.

## Reproducibility and version control

We aim to make all of our projects and analyses completely open and reproducible. The [workflowr][] package makes this aim easier to achieve by providing a set of tools and conventions for reproducible analysis workflows in R and the simultaneous building of a website that presents the analyses, with source code, in a readable way. 

The system integrates seamlessly with the git version control system.

As many project files as possible should be under version control. All user-created text files (e.g. `.md`, `.Rmd`, etc) and code (`.R`, `.py`, etc. files) should be under version control. Files that are large in size should not be under version control (too many large files in the repository makes git become unwieldy), and, in general, there is no need for files produced by code or analysis files to be versioned as running the workflow will produce them.

Docker containers encapsulating software environments enable reproducibility by allowing the same software and environments to be easily used by different people across different platforms. As such, we use them extensively.

Workflow management software makes a very big difference when trying to run complicated computational workflows and making them portable across local, cluster, and other HPC computing environments.

## Acknowledgements

This project is a [workflmiowr][] project. Making use of the workflowr package for reproducible analyses dictates certain structures for the project file.

[workflowr]: https://github.com/jdblischak/workflowr
