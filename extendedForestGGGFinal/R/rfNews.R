rfNews <- function() {
    newsfile <- file.path(system.file(package="extendedForestGGGFinal
tGGG"), "NEWS") #GGGmod
    file.show(newsfile)
}
