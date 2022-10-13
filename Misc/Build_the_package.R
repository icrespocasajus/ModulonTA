# creates description and namespace files
usethis::use_description()
usethis::use_namespace()

# Create R directory
# Go to the R package main directory
base::dir.create("R")

# creates Package-level documentation so you can run ?nameofpackage
usethis::use_package_doc()

# created README.Rmd for Github landing page
# an .Rbuildignore file gets created
usethis::use_readme_rmd()

#Update README
devtools::build_readme()

# creates news file
usethis::use_news_md()

# setup continuous integration via travis-ci
usethis::use_github_actions()

# sets up testing infrastructure
usethis::use_testthat()


# See available functions
pacman::p_functions(ModulonTA)


usethis::use_data_raw()
# this will setup the folders needed for the data and raw-data
usethis::use_data_raw()
#usethis::use_data(data, overwrite = TRUE)# This line should be included in an DATASET.R, which has to be executed and kept in data-raw: source("./DATASET.R")
# A second script (data.R) has to be created in R/data.R for te data documentation





# Vignettes


devtools::check()# WARNING! Be aware that using check from RStudio build panel is going to wipe out the doc folder; vignettes may stop working properly
usethis::use_vignette(name = "INTRODUCTION")

tools::buildVignettes(dir = ".", tangle=TRUE)
# Go to the project main directory
dir.create("inst")
dir.create("inst/doc")
file.copy(dir("vignettes", full.names=TRUE), "inst/doc", overwrite=TRUE)

# Commit and re-build

browseVignettes('ModulonTA')
vignette('INTRODUCTION','ModulonTA')

# Dependencies NOT WORKING!
usethis::use_package("corrplot", type = "Imports")

# Check the package
devtools::check()


#Update Readme.Rmd
devtools::build_readme()
