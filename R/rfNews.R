rfNews <- function() {
    newsfile <- file.path(system.file(package="extendedForestGGG"), "NEWS") #GGGmod
    file.show(newsfile)
}
