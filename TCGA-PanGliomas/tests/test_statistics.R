# Run from any directory: Rscript --vanilla tests/test_statistics.R
suppressPackageStartupMessages({library(dplyr);library(readr);library(jsonlite)})
arg <- grep('^--file=',commandArgs(),value=TRUE)
repo <- normalizePath(file.path(dirname(sub('^--file=','',arg)), '..'))
source(file.path(repo,'functions/R/statistics_helpers.R'))
# Load only the validator definition, without executing the production analyses.
expr <- parse(file.path(repo,'functions/R/recompute_extended.R'))
for(e in expr) if(is.call(e) && identical(e[[1]],as.name('<-')) && identical(e[[2]],as.name('regenerate'))) eval(e)
work <- file.path(repo,'outputs/test-statistics'); dir.create(work,recursive=TRUE,showWarnings=FALSE)
out <- file.path(work,'derived'); dir.create(out,showWarnings=FALSE)
rd <- function(p)read_tsv(p,show_col_types=FALSE,progress=FALSE)
checks <- coverage <- list(); results <- list()
test <- function(name,expr) {force(expr);results[[name]] <<- 'PASS'}
reject <- function(expr)stopifnot(inherits(tryCatch({force(expr);NULL},error=identity),'error'))
path <- file.path(work,'reference.tsv')
reference <- tibble(metric=c('A','B'),pvalue=c(1e-140,.03),q_value=c(2e-140,.03),value=c(NA_real_,4))
write_tsv(reference,path)
validate <- function(x)regenerate(path,x,'metric','test fixture','2 tests','synthetic fixture')
test('relative tolerance accepts TSV rounding',stopifnot(numeric_agreement(1e-140*(1+1e-10),1e-140,TRUE)))
test('zero cannot replace tiny probability',stopifnot(!numeric_agreement(0,1e-140,TRUE)))
test('zero reference remains exact',stopifnot(!numeric_agreement(1e-300,0,TRUE)))
test('missingness and finite checks',stopifnot(!numeric_agreement(c(NA,1),c(0,1)),!numeric_agreement(Inf,Inf),numeric_agreement(c(NA,1),c(NA,1))))
test('all probability column variants detected',stopifnot(all(probability_column(c('P','p.value','p_value','pvalue','padj_bh','padj_bh_all_tests','q_value','q_wilcox_primary','welch_pvalue','permutation_pvalue','global_range_permutation_pvalue','fdr','fdr_within_class'))),!probability_column('pct_group1')))
test('table rejects tiny P corruption before writing',{
  before<-tools::md5sum(path);a<-reference;a$pvalue[1]<-0;reject(validate(a));stopifnot(identical(tools::md5sum(path),before))
})
test('table rejects missing or duplicate comparison keys',{
  reject(validate(reference[1,]));a<-reference;a$metric[2]<-'A';reject(validate(a))
})
test('table rejects changed missingness',{
  a<-reference;a$value[1]<-0;reject(validate(a))
})
test('table aligns intact comparison keys',{
  a<-validate(reference[2:1,]);stopifnot(identical(a$metric,reference$metric))
})
# Very different distributions across metrics make a pooled BH correction wrong.
v <- data.frame(metric=rep(c('A','B'),each=15),group=rep(rep(c('G1','G2','G3'),each=5),2),
                value=c(1:5,6:10,11:15,1:5,2:6,3:7))
a <- pairwise_metrics(v,'metric','group',combn(c('G1','G2','G3'),2,simplify=FALSE),c('A','B'))
test('BH family boundaries include all three pairs',{
  stopifnot(nrow(a)==6,all(a$n_group_1==5),all(a$n_group_2==5),
            all(abs(a$padj_bh[a$metric=='A']-p.adjust(a$pvalue[a$metric=='A'],'BH'))<1e-14),
            !isTRUE(all.equal(a$padj_bh,p.adjust(a$pvalue,'BH'))))
})
d <- data.frame(Tumor_Barcode=letters[1:6],group=rep(c('G1','G2','G3'),each=2),
                weighted_time_sum=c(1,2,3,6,8,10),segment_weight_sum=rep(10,6),n_segments=1)
p <- weighted_permutation(d,'group',11,n_perm=99)
q <- weighted_permutation(d,'group',11,n_perm=99,plus_one=TRUE)
test('permutation weighting and plus-one convention',{
  stopifnot(abs(p$observed-.75)<1e-14,p$exceedances==q$exceedances,
            q$pvalue==(p$exceedances+1)/100,p$pvalue==p$exceedances/99)
})
test('permutation rejects duplicate specimens',{
  bad<-rbind(d,d[1,]);reject(weighted_permutation(bad,'group',11,n_perm=9))
})
write_json(list(status='PASS',tests=results),file.path(work,'results.json'),pretty=TRUE,auto_unbox=TRUE)
cat(length(results),'statistical regression tests PASS\n')
