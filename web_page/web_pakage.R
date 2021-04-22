library(pkgdown)
sink("_pkgdown.yml")
template_navbar()
template_reference()
template_articles()
#init_site()
sink()
pkgdown::build_site()


build_news()
build_home()
build_articles() 


# usethis::use_pkgdown()
# pkgdown::template_navbar()
# pkgdown::build_site()

