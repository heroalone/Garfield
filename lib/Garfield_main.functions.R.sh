#!/usr/bin/Rscript

args <- commandArgs()
PLINKFILE <- args[6]
TRAITFILE <- args[7]
PREFIX <- args[8]
NONEGATIVETRAIT = as.numeric(args[9])


# # # module load r/3.5.1-foss-2018b
# FILENAME <- paste(INDEX, ".", TRAIT, ".bed", BED, sep = "")
# PLINKFILE2 <- gsub(paste(TRAIT, ".", sep = ""), "", PLINKFILE)
FILENAME2 <- paste(PREFIX, sep = "")

# Function to read genotype in PLINK format
read_plink_file <- function(PLINKFILE) {
	suppressMessages(library("genio"))
	suppressMessages(data_genio <- read_plink(paste(PLINKFILE, sep = "")))
	id <- colnames(data_genio$X)
	geno <- as.data.frame(t(data_genio$X) * 0.5, row.names = FALSE)
	geno_data <- as.data.frame(cbind(id, geno))
	CHR <- unique(data_genio$bim$chr)
	return(list(geno_data = geno_data, CHR = CHR))
}

# Function to read phenotype and process it into binary types
read_phenotype_file <- function(TRAITFILE, NONEGATIVETRAIT = 1) {
	pheno_data <- read.table(paste(TRAITFILE, sep = ""), header = FALSE, sep = " ", skipNul = TRUE)
	if (NONEGATIVETRAIT){
		if ((!any(is.na(pheno_data[, 3]))) && any(pheno_data[, 3] < 0)) {
			suppressWarnings(pheno_data[which(pheno_data[, 3] < 0), ][, 3] <- as.numeric("NA"))
		}
	}
	names(pheno_data) <- c("fid", "id", "phe")
	phenotype <- na.omit(pheno_data[, -1])
	cluster_assignment <- kmeans(phenotype$phe, 2)
	phenotype <- data.frame(id = phenotype$id, phe = cluster_assignment$cluster - 1)
	return(phenotype)
}

# Function to perform random forest analysis and get feature importance
random_forest_getImp <- function(data, minnode = 5, seed = 222) {
	suppressMessages(library("ranger"))
	rf_model <- ranger(factor(data$phe) ~ ., data = data, importance = "permutation", write.forest = TRUE, min.node.size = minnode, seed = seed)
	feature_importance <- data.frame(variables = names(importance(rf_model)), IMP = importance(rf_model))
	feature_importance <- feature_importance[order(feature_importance$IMP, decreasing = TRUE), ]
	positive_imp_features <- feature_importance[which(feature_importance$IMP > 0), ]
	return(list(rf_ImpPositive = positive_imp_features, rf_rawImp = feature_importance))
}

# Function to select get logic expressions
logicFS_getImp <- function(data, pheno_ID, B=100, nleaves = 10, ntrees = 1, seed = 222, start = 2, end = -2, iter = 1000){
	suppressMessages(library("logicFS"))
	binary_snps <- as.matrix(data[, -1])
	suppressWarnings(logic_bag <- logic.bagging(binary_snps, factor(data$phe), B = B, nleaves = nleaves, ntrees = ntrees, rand = seed, anneal.control = logreg.anneal.control(start = start, end = end, iter = iter), addMatImp = TRUE))
	predicted_genotype <- cbind(data.frame(id = pheno_ID, predict = predict(logic_bag, binary_snps))) #### predicGeno by LogicFS
	variable_importance <- data.frame(importance = logic_bag$vim$vim, expression = logic_bag$vim$primes)
	variable_importance <- variable_importance[order(-variable_importance$importance), ][variable_importance$importance > 0, ]
	# variable_importance <- variable_importance[variable_importance$importance > 0, ]
	variable_importance <- variable_importance[!grepl("TEMP_ABCDEFG", variable_importance$expression), ]
	return(list(variable_importance = variable_importance, predicted_genotype = predicted_genotype))
}

# Function to select top variable based on importance
select_top_variable <- function(data) {
	clustering <- kmeans(data, 2)
	min_change_point3 <- which(abs(diff(clustering$cluster)) == max(abs(diff(clustering$cluster))))
	min_change_point1 <- which(abs(diff(data, difference = 1)) == max(abs(diff(data, difference = 1))))
	
	if (max(clustering$cluster) > 1) {
		min_change_point2 <- which(abs(diff(data, difference = 2)) == max(abs(diff(data, difference = 2))))
	}
	min_change_point <- ifelse(min_change_point3 > 1, min_change_point3, max(min_change_point1, min_change_point2))
	return(min_change_point)
}

select_top_variable2 <- function(data, maxLogic = 3) {
	max_change_point1 <- which(abs(diff(data)) == max(abs(diff(data))))
	max_change_point2 <- which(abs(diff(data, difference = 2)) == max(abs(diff(data, difference = 2))))
	max_change_point <- max(max_change_point1, max_change_point2)
	min_change_point <- ifelse(max_change_point < maxLogic, max_change_point, maxLogic)
	return(min_change_point)
}

# Function to write TPED file
write_TPED <- function(CHR, FILENAME2, Value) {
	tped_OUT <- paste(CHR, FILENAME2, "0", "1", paste(Value, Value, sep = " ", collapse = " "), sep = " ")
	return(tped_OUT)
}

write_bestDNF <- function(CHR, FILENAME2, Value) {
	dnf_OUT <- paste(CHR, FILENAME2, Value, sep = "\t")
	return(dnf_OUT)
}


OUT_TPED="NULL"
# run Garfield to produce new tped and bestDNF files
Garfield_main <- function(PLINKFILE, TRAITFILE, FILENAME2) {
	# Read files
	geno <- read_plink_file(PLINKFILE)
	genotype <- geno$geno_data
	CHR <- geno$CHR
	phenotype <- read_phenotype_file(TRAITFILE)
	
	# Combine phenotype and genotype data
	data <- cbind(phenotype, genotype)[, c(-1, -3)]
	
	# Perform random forest analysis and get feature importance
	Imp_RF <- random_forest_getImp(data)
	importance_RF <- Imp_RF$rf_ImpPositive
	raw_importance <- Imp_RF$rf_rawImp
	rm(list=c('Imp_RF', 'geno'))


	if (length(importance_RF$IMP) < 3) {
		if (length(importance_RF$IMP) == 0) {
			importance_RF <- raw_importance[1, ]
		}
		# # OUT_TPED <- write_TPED(CHR, FILENAME2, genotype[, as.vector(importance_RF$variables[1])] + 1)
		OUT_DNF <- write_bestDNF(CHR, FILENAME2, importance_RF$variables[1])
	} else {
		if (length(importance_RF$IMP) > 5) {
			importance_RF <- importance_RF[order(importance_RF$IMP, decreasing = TRUE), ][c(1:5), ]
		}
		min_change_point <- select_top_variable(importance_RF$IMP)
		if (min_change_point == 1) {
			# # OUT_TPED <- write_TPED(CHR, FILENAME2, genotype[, as.vector(importance_RF$variables[1])] + 1)
			OUT_DNF <- write_bestDNF(CHR, FILENAME2, importance_RF$variables[1])
		} else {
			top_best <- data[, c("phe", as.vector(importance_RF$variables[1:min_change_point]))]
			if (min_change_point == 2) {
				top_best$TEMP_ABCDEFG <- 1
			}
			# Perform logic analysis to get feature importance
			LFS <- logicFS_getImp(top_best, phenotype$id)
			vim.result <- LFS$variable_importance
			predicted_genotype <- LFS$predicted_genotype
			
			if (nrow(vim.result) == 0 || unique(factor(predicted_genotype$predict)) == 1) {
				# # OUT_TPED <- write_TPED(CHR, FILENAME2, genotype[, as.vector(importance_RF$variables[1])] + 1)
				OUT_DNF <- write_bestDNF(CHR, FILENAME2, importance_RF$variables[1])
			} else {
				if (nrow(vim.result) == 1) {
					sig_logic_exp = as.vector(vim.result$expression[1])
					OUT_DNF <- write_bestDNF(CHR, FILENAME2, sig_logic_exp)
				} else if (nrow(vim.result) == 2) {
					min_change_point = 2
					sig_logic_exp <- paste(vim.result$expression[1:min_change_point], sep = "", collapse = " || ")
					OUT_DNF <- write_bestDNF(CHR, FILENAME2, sig_logic_exp)
				} else {
					min_change_point <- select_top_variable2(vim.result$importance)
					sig_logic_exp <- paste(vim.result$expression[1:min_change_point], sep = "", collapse = " || ")
					OUT_DNF <- write_bestDNF(CHR, FILENAME2, sig_logic_exp)
				}
				
				if (nrow(vim.result)> 1) {
					merged_data <- merge(phenotype, predicted_genotype, by = "id", sort = FALSE, row.names = phenotype$id, all = TRUE)
					merged_data <- merged_data[match(phenotype$id, merged_data$id), ]
					merged_data$predict <- as.numeric(merged_data$predict)
					if (any(is.na(merged_data$predict))) {
						merged_data[is.na(merged_data$predict), ]$predict <- -1
					}
					
					OUT_TPED <- write_TPED(CHR, FILENAME2, merged_data$predict + 1)
				}
			}
		}
	}
	return(list=c(OUT_TPED, OUT_DNF))
}

Garfield_results <- Garfield_main(PLINKFILE, TRAITFILE, FILENAME2)
# Garfield_results <- Garfield_main(PLINKFILE, TRAITFILE, FILENAME2)
cat(Garfield_results, sep = "\n")
