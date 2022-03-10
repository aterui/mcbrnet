
# build package -----------------------------------------------------------

  #usethis::use_mit_license(name = "Akira Terui")
  usethis::use_roxygen_md()
  usethis::use_package_doc()
  devtools::document()
  devtools::load_all()
  devtools::check()


# check syntax ------------------------------------------------------------

  lintr::lint_package()

  #usethis::use_coverage()
  devtools::build(path='.')
  covr::package_coverage()
  file.remove("mcbrnet_1.2.1.tar.gz")
