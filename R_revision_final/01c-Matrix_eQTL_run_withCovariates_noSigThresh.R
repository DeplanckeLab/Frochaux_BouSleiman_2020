library(MatrixEQTL)

### This is an eQTL analysis that does not enforce a significance threshold on the cis associations (pvOutputThreshold.cis = 1)
### This is used for the power calculation script

#########################################
###### For NAIVE/Control condition ######
#########################################

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = "./Data/genoForMatrixeQTL.noS14.noS54.Naive.txt";

# Gene expression file name
expression_file_name = "./Data/cpm.table.Naive.txt";

# Gene position file name
geneposition_file_name = "./Data/genlocs.txt";

# snp position file name
snpposition_file_name = "./Data/snpsloc.txt"

# Covariates file name
# Set to character() for no covariates
covariates_file_name = "./Data/covariates.pc.wolb.inv.Naive.txt";

# Output file name
output_file_name = "./Data/MatrixEQTL_noSigThresh/MatrixeQTL.Naive.output.withCovariates";

# Only associations significant at this level will be saved
#pvOutputThreshold = 1e-5; #when set to zero, only cis is performed
pvOutputThreshold = 0;
pvOutputThreshold.cis = 1;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = " ";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 5000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);
## transforming to standard normal distribution
for( sl in 1:length(gene) ) {
  mat = gene[[sl]];
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(gene)+1));
  gene[[sl]] = mat;
}
rm(sl, mat);


## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = " ";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 5000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene position file
geneposition = read.table(geneposition_file_name, header=TRUE)
## Load snp position file
snpposition = read.table(snpposition_file_name, header=TRUE)

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = ",";      # the separator character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}
## Run the analysis
#memory.limit(size= 4000)


me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = paste(output_file_name, "cis", sep="."),
  pvOutputThreshold.cis = pvOutputThreshold.cis,
  snpspos = snpposition,
  genepos = geneposition,
  cisDist = 10e3, #10 kb
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);
#unlink(output_file_name);
## Results:
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)
## Plot the histogram of all p-values
# plot(me)


##########################################
##### For TREATED/INFECTED condition #####
##########################################

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = "./Data/genoForMatrixeQTL.noS14.noS54.Treated.txt";

# Gene expression file name
expression_file_name = "./Data/cpm.table.Treated.txt";

# Gene position file name
geneposition_file_name = "./Data/genlocs.txt";

# snp position file name
snpposition_file_name = "./Data/snpsloc.txt"

# Covariates file name
# Set to character() for no covariates
covariates_file_name = "./Data/covariates.pc.wolb.inv.Treated.txt";

# Output file name
output_file_name = "./Data/MatrixEQTL_noSigThresh/MatrixeQTL.Treated.output.withCovariates";

# Only associations significant at this level will be saved
#pvOutputThreshold = 1e-5; #when set to zero, only cis is performed
pvOutputThreshold = 0;
pvOutputThreshold.cis = 1;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = " ";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 5000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);
## transforming to standard normal distribution
for( sl in 1:length(gene) ) {
  mat = gene[[sl]];
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(gene)+1));
  gene[[sl]] = mat;
}
rm(sl, mat);


## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = " ";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 5000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);


## Load gene position file
geneposition = read.table(geneposition_file_name, header=TRUE)
## Load snp position file
snpposition = read.table(snpposition_file_name, header=TRUE)

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = ",";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
#memory.limit(size= 4000)

###association
me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = paste(output_file_name, "cis", sep="."),
  pvOutputThreshold.cis = pvOutputThreshold.cis,
  snpspos = snpposition,
  genepos = geneposition,
  cisDist = 10e3, #10 kb
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

#unlink(output_file_name);


## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)

## Plot the histogram of all p-values

# plot(me)