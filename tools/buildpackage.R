
# build package -----------------------------------------------------------

  #usethis::use_mit_license(name = "Akira Terui")
  usethis::use_roxygen_md()
  usethis::use_package_doc()
  devtools::document()
  devtools::load_all()
  devtools::check()


# check syntax ------------------------------------------------------------

  #usethis::use_coverage()
  covr::package_coverage()
