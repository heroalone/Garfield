cat("Loading R libraries\n")

if(!suppressMessages(library(genio, warn.conflicts=F))){
	install.packages("genio", dep=T)
	if(!suppressMessages(library(genio, warn.conflicts=F))) {
		stop("Install genio R package failed! You may have to install it manually.\n")
	}
}


if(!suppressMessages(library(ranger, warn.conflicts=F))){
	install.packages("ranger", dep=T)
	if(!suppressMessages(library(ranger, warn.conflicts=F))) {
		stop("Install ranger R package failed! You may have to install it manually.\n")
	}
}


if(!suppressMessages(library(logicFS, warn.conflicts=F))){
	if (!require("BiocManager", quietly = TRUE))
		install.packages("BiocManager")

	BiocManager::install("logicFS")

	if(!suppressMessages(library(logicFS, warn.conflicts=F))) {
		stop("Loading package cowplot failed! Please install manually.\n")
	}
}

cat("R packages installed succesfully\n")

