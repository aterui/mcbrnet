
# build package -----------------------------------------------------------

#usethis::use_mit_license(copyright_holder = "Akira Terui")
usethis::use_roxygen_md()
usethis::use_package_doc()
devtools::document()
devtools::load_all()
devtools::check(vignettes=FALSE)


# check syntax ------------------------------------------------------------

# lintr::lint_package()

#usethis::use_coverage()
# devtools::build(path='.')
# covr::package_coverage()
# file.remove("mcbrnet_1.2.3.tar.gz")


# build website -----------------------------------------------------------

# Run once to configure package to use pkgdown
#usethis::use_pkgdown()

# Run to build the website
pkgdown::build_site()
