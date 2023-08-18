cat("Installing R libraries\n")

install_package <- function(package_name) {
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Install", package_name, "R package failed! You may have to install it manually.\n"))
    }
  }
}

install_package("genio")
install_package("ranger")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("logicFS")


# # # check_and_install_package("logicFS")
required_packages <- c("genio", "ranger", "logicFS")

recheck_packages <- function(package_names) {
  missing_packages <- package_names[!sapply(package_names, requireNamespace, quietly = TRUE)]
  return(missing_packages)
}


missing_packages <- recheck_packages(required_packages)

if (length(missing_packages) == 0) {
  cat("Required R packages have been successfully installed!\n")
} else {
  cat("The installation of the following packages has failed. Please install them manually.\n")
  cat(paste(missing_packages, collapse = ", "), "\n")
}

