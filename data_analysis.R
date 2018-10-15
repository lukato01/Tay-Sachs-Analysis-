suppressMessages(library(janitor))
suppressMessages(library(vcfR))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(rmarkdown))
suppressMessages(library(rJava))
suppressMessages(library(openxlsx))
suppressMessages(library(lubridata))
suppressMessages(library(rio))
suppressMessages(library(data.table))
suppressMessages(library(genetics))
suppressMessages(library(lawstat))
suppressMessages(library(car))
suppressMessages(library(rpart))
suppressMessages(library(gdata))
suppressMessages(library(rio)) 
suppressMessages(library(tidyr))
suppressMessages(library(lubridate))
suppressMessages(library(rlist))
suppressMessages(library(janitor))
suppressMessages(library(SNPassoc))
suppressMessages(library(data.table))

#plink 1.9 to convert and plink 1.07 for analysis 

custom_fun <- function(myvector) {
  return(myvector[1])
}

custom_fun2 <- function(myvector) {
  return(myvector[2])
}
custom_fun3<- function(myvector) {
  return(myvector[3])
}

#function to perform GLM test for every SNP
GLM_function15 <- function(variant){
  local1<- lm(pheno_15~(variant)+sex_15+race_15+year_15, na.action=na.exclude)
  a<-summary(local1)
  return(a$coefficients[2,4])
}

beta_15 <- function(variant){
  local1<- lm(pheno_15~(variant)+sex_15+race_15+year_15, na.action=na.exclude)
  a<-summary(local1)
  return(a$coefficients[2,1])
}

SE_15 <- function(variant){
  local1<- lm(pheno_15~(variant)+sex_15+race_15+year_15, na.action=na.exclude)
  a<-summary(local1)
  return(a$coefficients[2,2])
}


#function to perform anova test for every SNP
GLM_function5 <- function(par2){ 
  local1<- lm(pheno_5~(par2)+sex_5+race_5+year_5, na.action=na.exclude)
  a<-summary(local1)
  return(a$coefficient[2,4])
}


beta_5 <- function(variant){
  local1<- lm(pheno_5~(variant)+sex_5+race_5+year_5, na.action=na.exclude)
  a<-summary(local1)
  return(a$coefficients[2,1])
}

SE_5 <- function(variant){
  local1<- lm(pheno_5~(variant)+sex_5+race_5+year_5, na.action=na.exclude)
  a<-summary(local1)
  return(a$coefficients[2,2])
}

check_integer_mean <- function(vec,dataframe) {
  for(i in vec) {
    a <- as.data.table(table(dataframe[,i]), keep.rownames=TRUE)
    
    ## mean column, calculate and add
    average <- aggregate( as.formula(paste("Result ~ `", i, "`", sep="")), dataframe, mean )
    
    # Rename the column V1 
    setnames(a, old = 'V1', new = i)
    
    # bind the 'Result' (avg) from blah into a
    new <- merge(average, a,  by =i)
    
    
    # if ( !(1 %in% new$N) & !(2 %in% new$N) & !(3 %in% new$N) ) {
      print (new)
    }
  }
# }

check_integer_mean2 <- function(vec,dataframe) {
  df_combined <- data.frame(Variant=character(),
                            Genotype=character(), 
                            EnzymeActivity=character(), 
                            N=character()) 
  
  for(i in vec) {
    a <- as.data.table(table(dataframe[,i]), keep.rownames=TRUE) # A is the dataframe with single column
    ## mean column, calculate and add
    average <- aggregate( as.formula(paste("Result ~ `", i, "`", sep="")), dataframe, mean )
    # Rename the column V1 
    setnames(a, old = 'V1', new = i)
    # bind the 'Result' (avg) from blah into a
    new <- merge(average, a,  by =i)
    ## Some more table manipulation here...
    # add a column Variant, as a value put the name of the Genotype column
    new <- cbind(Variant = colnames(new)[1] , new)
    # rename the Genotype columns 
    colnames(new)[2] <- "Genotype"
    colnames(new)[3] <- "EnzymeActivity"
    colnames(new)[4] <- "N"
    
    df_combined <- rbind(df_combined, new)
  }
  return(df_combined)
}


check_integer_mean2_sd <- function(vec,dataframe) {
  df_combined <- data.frame(Variant=character(),
                            Genotype=character(), 
                            SD=character(),
                            AvgEnzymeActivity=character(), 
                            N=character())
                            
  
  for(i in vec) {
    a <- as.data.table(table(dataframe[,i]), keep.rownames=TRUE)
    
    ## mean column, calculate and add
    average <- aggregate( as.formula(paste("Result ~ `", i, "`", sep="")), dataframe, mean )
    ## sd column, calculate and add
    mysd <- aggregate( as.formula(paste("Result ~ `", i, "`", sep="")), dataframe, sd )
    
    # Rename the column V1 
    setnames(a, old = 'V1', new = i)
    # bind the 'Result' (avg) from blah into a
    other <- merge(mysd, average,  by =i)
    new <- merge( other, a , by =i)
    
   
    
    ## Some more table manipulation here...
    # add a column Variant, as a value put the name of the Genotype column
    new <- cbind(Variant = colnames(new)[1] , new)
    # rename the Genotype columns 
    colnames(new)[2] <- "Genotype"
    colnames(new)[3] <- "SD"
    colnames(new)[4] <- "AvgEnzymeActivity"
    colnames(new)[5] <- "N"
    
    df_combined <- rbind(df_combined, new)
  }
  return(df_combined)
}



#### Subset the molecular data into table with heterozygous mutation & homozygous mutations ----

heterozygous_mutation<-subset(molecular, (
  gt_GT == "0/1"|
    gt_GT == "0/2"|
    gt_GT == "0/3"|
    gt_GT == "0/4"|
    gt_GT == "0/5"|
    gt_GT == "0/6"|
    gt_GT == "0/7"|
    gt_GT == "0/8"|
    gt_GT == "0/9"|
    gt_GT == "0/10" ))

homozygous_mutations<-subset(molecular, (
  gt_GT!= "0/0" &
    gt_GT!= "0/1" &
    gt_GT!= "0/2" &
    gt_GT != "0/3" &
    gt_GT != "0/4" &
    gt_GT != "0/5" &
    gt_GT != "0/6" &
    gt_GT != "0/7" &
    gt_GT!= "0/8" &
    gt_GT != "0/9" &
    gt_GT != "0/10"))

# Heterozygous mutations will be merged to the classification table by looking for the common position and genotype

mergedData01 <- merge(heterozygous_mutation, complete_list_of_classified_variants_within_target, by.x=c("ChromPOS","gt_GT_alleles"),
                      by.y=c('ChromPOS', "Genotype"))
# homozygous mutations will be merged to the classification table by the common position and commor REF or ALT 
###seperating the genotype column into ref and alt so that the # homozygous mutation table could be compared and combined by these newly formed columns
homozygous_mutations_REFALT<-homozygous_mutations %>% separate(gt_GT_alleles, c('Ref', 'Alt'), sep="/")

mergedData_1_a <- merge(homozygous_mutations_REFALT, complete_list_of_classified_variants_within_target, by.x=c("ChromPOS","Alt"),
                        by.y=c('ChromPOS', "ALT" ) )

mergedData_1_b <- merge(homozygous_mutations_REFALT, complete_list_of_classified_variants_within_target, by.x=c("ChromPOS","Ref"),
                        by.y=c('ChromPOS', "ALT"))
#select unquie rows and get rid of all the duplicates 
mergedData_1_1<-rbind(mergedData_1_a, mergedData_1_b) %>% distinct

#removing and renaming columns in order to combine two tables 
mergedData_1_1<- subset(mergedData_1_1, select=-Ref)
mergedData_1_1<- rename.vars( mergedData_1_1, c("Genotype" ,"Alt"), c("gt_GT_alleles", "ALT"  ))
combined_tables<-rbind(mergedData01, mergedData_1_1)

#### Number of patients with mutations in the classified variant table ----

# get frequency for counts for cDNA by ChromPOS
table<-data.table(table(combined_tables$cDNA, combined_tables$ChromPOS))
#get rid of all the rows that have 0
table[table == 0] <- NA
s<-na.omit(table)
#merge two tables in order to get the frequencies
pleasework <- merge(complete_list_of_classified_variants_within_target, s , by.x=c("ChromPOS","cDNA"),
                    by.y=c('V2', "V1"), all.x=TRUE)
# substitute all NA for 0
pleasework[is.na(pleasework)] <- 0
# order the rows
d<-pleasework %>% separate(ChromPOS, c('Chrom', 'POS'), sep=":")
d<-arrange(d,c(Chrom))
d<-arrange(d,c(POS))
d<- d %>% unite(ChromPOS,Chrom,POS, sep=":", remove=TRUE)
d<-subset(d, N!=0)

#remove rows

export(list(d=d), "summary_of_variants.xlsx")

t1 <- subset(molecular, (ChromPOS == "15:72636373" & gt_GT != "0/0"))


#### Summary Table ----

# HEXA gene is on chromosome 15
# HEXB gene is on chromosome 5


merged_data= merge(biochem, combined_tables, by='Indiv',all = FALSE, sort = FALSE)
merged_data<-merged_data%>% separate(ChromPOS, c('Chrom', 'POS'), sep=":")
merged_data$Classification <- as.character(merged_data$Classification)

reclassified_vars <- within(merged_data, Classification[Classification == 'Likely pathogenic'] <- 'Pathogenic')
reclassified_vars <- within(reclassified_vars, Classification[Classification == 'pathogenic'] <- 'Pathogenic')
reclassified_vars <- within(reclassified_vars, Classification[Classification == 'Likely benign'] <- 'Benign')
reclassified_vars <- within(reclassified_vars, Classification[Classification == 'Likely Benign'] <- 'Benign')

# No HEXA/HEXB variants or only benign/likely benign variants found
not_benign<-subset(reclassified_vars, Classification=="Pathogenic"|Classification=="Pseudodeficienct"|Classification=="VOUS" )
# samples that are in non benign are excluded and we are left with samples that are only have benigh mutations
benign <- reclassified_vars[ ! reclassified_vars$Indiv %in% not_benign$Indiv, ]
kk<-subset(benign,   select=c(Indiv, biochem_groups, Classification))
# unique the samples that are benign
only_benign<-kk[!duplicated(kk),]
table(only_benign$biochem_groups)


# Sandhoff_Carrier_only_benign<- subset(only_benign, (biochem_groups =="Sandhoff-Carrier"))
# H5 <- molecular[molecular$Indiv %in% Sandhoff_Carrier_only_benign$Indiv, ]
# H5<- subset(H5, (gt_GT !="0/0"))
# write.xlsx(H5, file= "H5.xlsx",  row.names=TRUE)

# figure out which samples are duplicated in the non benign dataframe and these are the ones
# that have more then one mutation and exclude those from the list
bb<-subset(not_benign,   select=c(Indiv))
more_then_one_mutation<-bb[duplicated(bb),]


#exclude sample that have more then one mutation. these samples only 1 specific classification 
not_benign_with_one_mutation <- not_benign[ ! not_benign$Indiv %in% more_then_one_mutation$Indiv, ]

#15
Pathogenic_15<-subset(not_benign_with_one_mutation, Classification=="Pathogenic"  & Chrom=="15")
table(Pathogenic_15$biochem_groups)

#getting inconclusive 
# Inconclusive_Pathogenic_15<- subset(Pathogenic_15, (biochem_groups =="Inconclusive"))
# B3 <- molecular[molecular$Indiv %in% Inconclusive_Pathogenic_15$Indiv, ]
# B3<- subset(B3, (gt_GT !="0/0"))
# write.xlsx(B3, file= "B3.xlsx",  row.names=TRUE)

# Non_Carrier_Pathogenic_15<- subset(Pathogenic_15, (biochem_groups =="Non-Carrier"))
# B4 <- molecular[molecular$Indiv %in% Non_Carrier_Pathogenic_15$Indiv, ]
# B4<- subset(B4, (gt_GT !="0/0"))
# write.xlsx(B4, file= "B4.xlsx",  row.names=TRUE)


Pseudodeficienct_15<-subset(not_benign_with_one_mutation, Classification=="Pseudodeficienct"  & Chrom=="15")
table(Pseudodeficienct_15$biochem_groups)

# Non_Carrier_Pseudodeficienct_15<- subset(Pseudodeficienct_15, (biochem_groups =="Non-Carrier"))
# C4 <- molecular[molecular$Indiv %in% Non_Carrier_Pseudodeficienct_15$Indiv, ]
# C4<- subset(C4, (gt_GT !="0/0"))
# write.xlsx(C4, file= "C4.xlsx",  row.names=TRUE)

VOUS_15<-subset( not_benign_with_one_mutation, Classification=="VOUS"  & Chrom=="15")
table(VOUS_15$biochem_groups)

# Tay_Sachs_Carrier_VOUS_15<- subset(VOUS_15, (biochem_groups =="Tay-Sachs Carrier"))
# D2 <- molecular[molecular$Indiv %in% Tay_Sachs_Carrier_VOUS_15$Indiv, ]
# D2 <- subset(D2, (gt_GT !="0/0"))
# write.xlsx(D2, file= "D2.xlsx",  row.names=TRUE)

#5
Pathogenic_5<-subset(not_benign_with_one_mutation, Classification=="Pathogenic"  & Chrom=="5")
table(Pathogenic_5$biochem_groups)

# non_Carrier_Pathogenic_5<- subset(Pathogenic_5, (biochem_groups =="Non-Carrier"))
# E4 <- molecular[molecular$Indiv %in% non_Carrier_Pathogenic_5$Indiv, ]
# E4<- subset(E4, (gt_GT !="0/0"))
# write.xlsx(E4, file= "E4.xlsx",  row.names=TRUE)


Pseudodeficienct_5<-subset(not_benign_with_one_mutation, Classification=="Pseudodeficienct"  & Chrom=="5")
table(Pseudodeficienct_5$biochem_groups)

VOUS_5<-subset(not_benign_with_one_mutation, Classification=="VOUS"  & Chrom=="5")
table(VOUS_5$biochem_groups)


# Sandhoff_Carrier_VOUS_5<- subset(VOUS_5, (biochem_groups =="Sandhoff-Carrier"))
# G5 <- molecular[molecular$Indiv %in% Sandhoff_Carrier_VOUS_5$Indiv, ]
# G5<- subset(G5, (gt_GT !="0/0"))
# write.xlsx(G5, file= "G5.xlsx",  row.names=TRUE)

# samples that have at least 2 mutations 
multiple_mutations <- not_benign[not_benign$Indiv %in% more_then_one_mutation$Indiv, ]

multiple_mutations_HEXA<-subset(multiple_mutations, grepl("^15", Chrom))
multiple_mutations_HEXB<-subset(multiple_mutations, grepl("^5", Chrom))

# samples that have one hex a and hexb mutation 
common_Indiv<-data.frame(intersect(multiple_mutations_HEXA$Indiv,multiple_mutations_HEXB$Indiv))
mutations_in_HEX_A_and_B <- multiple_mutations [multiple_mutations$Indiv %in% common_Indiv$intersect.multiple_mutations_HEXA.Indiv..multiple_mutations_HEXB.Indiv., ]

ww<-subset(mutations_in_HEX_A_and_B,   select=c(Indiv, biochem_groups))
ff<-ww[!duplicated(ww),]
table(ff$biochem_groups)
#
# multiple_mutations <- multiple_mutations[!multiple_mutations$Indiv %in% mutations_in_HEX_A_and_B $Indiv, ]
# mutations in either hex a or hex b
multiple_mutations <- multiple_mutations [!multiple_mutations$Indiv %in% common_Indiv$intersect.multiple_mutations_HEXA.Indiv..multiple_mutations_HEXB.Indiv., ]


multiple_mutations_HEXA<-subset(multiple_mutations, grepl("^15", Chrom))
multiple_mutations_HEXA<-unique(multiple_mutations_HEXA, by=c(1:2))
table(multiple_mutations_HEXA$biochem_groups)

multiple_mutations_HEXB<-subset(multiple_mutations, grepl("^5", Chrom))
multiple_mutations_HEXB<-unique(multiple_mutations_HEXB, by=c(1:2))
table(multiple_mutations_HEXB$biochem_groups)

  

#### Preparing data for analysis(removing Vous, pathogenic and HEXA polymorpisms) ----
  
# Removing patients containing two HEXA polymorphisms("T/T" or "T/C" genotype) chr15:72638892 and chr15:72637795 that decrease biochemical enzyme activity 
  
# sub1_special is all those rows that match the conditions
sub1_special <- subset(molecular, 
                         (gt_GT_alleles == "T/T" & ChromPOS == "15:72637795") |
                           (gt_GT_alleles == "T/C" & ChromPOS == "15:72637795") |
                           (gt_GT_alleles == "T/T" & ChromPOS == "15:72638892") |
                           (gt_GT_alleles == "T/C" & ChromPOS == "15:72638892") 
)
# Find the unique elements in 'Indiv' column of sub1_special data frame
unique_indivs_to_delete <- sub1_special[!duplicated(sub1_special[,c('Indiv')]),]
# Now go back and filter the original table based on these unique Indivs
molecular_without_bad_indivs <- molecular[ ! molecular$Indiv %in% unique_indivs_to_delete$Indiv, ]
  
# Selecting patience that have genotype 0/1 for all the variance 
try_0_1<-subset(molecular_without_bad_indivs, (gt_GT == "0/1"|
                                                   gt_GT == "0/2"| 
                                                   gt_GT == "0/3"|
                                                   gt_GT == "0/4"| 
                                                   gt_GT == "0/5"|
                                                   gt_GT == "0/6"|
                                                   gt_GT == "0/7"|
                                                   gt_GT == "0/8"|
                                                   gt_GT == "0/9"|
                                                   gt_GT == "0/10" ))
  
# Merging patience that have genotype 0/1 to the variance  with classification combining these two tables by the locus and  specific getoype in that position 
mergedData_0_1 <- merge(try_0_1, complete_list_of_classified_variants_within_target, by.x=c("ChromPOS","gt_GT_alleles"),
                          by.y=c('ChromPOS', "Genotype"))
  
# Filtering the table with genotype 0/1 by getting rid of all the patience that have "Pathogenic", "Likely pathogenic", or "VOUS" mutations
sub_special_0_1 <- subset(mergedData_0_1, 
                            ( Classification == "Pathogenic") |
                              (Classification == "Likely pathogenic") |
                              (Classification == "VOUS")) 
# Now go back and filter the original table based on these unique Indivs that were removed 
molecular_without_path_or_vous_n1 <- molecular_without_bad_indivs[ ! molecular_without_bad_indivs$Indiv %in% sub_special_0_1$Indiv, ]
  
# Selecting patience that have genotype 1/1 for all the variance
try_1_1<-subset(molecular_without_bad_indivs, (gt_GT!= "0/0" & 
                                                   gt_GT!= "0/1" & 
                                                   gt_GT!= "0/2" &
                                                   gt_GT != "0/3" & 
                                                   gt_GT != "0/4" & 
                                                   gt_GT != "0/5" &
                                                   gt_GT != "0/6" & 
                                                   gt_GT != "0/7" &
                                                   gt_GT!= "0/8" & 
                                                   gt_GT != "0/9" & 
                                                   gt_GT != "0/10"))
# Seperating the genotype column into ref and alt
new_1_1<-try_1_1 %>% separate(gt_GT_alleles, c('Ref', 'Alt'), sep="/")
  
# merging patience that have genotype 1/1 for all the variance  with classification of variance and combining these two tables by the chromosome position and by alt genotype
mergedData_1_a <- merge(new_1_1, complete_list_of_classified_variants_within_target, by.x=c("ChromPOS","Alt"),
                          by.y=c('ChromPOS', "ALT" ) )
  
mergedData_1_b <- merge(new_1_1, complete_list_of_classified_variants_within_target, by.x=c("ChromPOS","Ref"),
                          by.y=c('ChromPOS', "ALT"))
ii<-rbind(mergedData_1_a, mergedData_1_b)
mergedData_1_1<-ii %>% distinct
  
sub_special_1_1 <- subset(mergedData_1_1, 
                            ( Classification == "Pathogenic") |
                              (Classification == "Likely pathogenic") |
                              (Classification == "VOUS")) 
# Now go back and filter the original table based on these unique Indivs that were removed 
molecular_without_path_or_vous <- molecular_without_path_or_vous_n1[ ! molecular_without_path_or_vous_n1$Indiv %in% sub_special_1_1$Indiv, ]
  
# molecular_without_path_or_vous_5 <- subset(molecular_without_path_or_vous, grepl("^5", ChromPOS))
# molecular_without_path_or_vous_15 <- subset(molecular_without_path_or_vous, grepl("^15", ChromPOS))
 
# merging molecular data without patients containing two HEXA polymorphisms and classifeid variants
# datamerge2_5<-merge(biochem,molecular_without_path_or_vous_5, by = "Indiv", all = FALSE, sort = FALSE)
# datamerge2_15<-merge(biochem,molecular_without_path_or_vous_15, by = "Indiv", all = FALSE, sort = FALSE)

molecular_5 <- subset(molecular, grepl("^5", ChromPOS))
molecular_15 <- subset(molecular, grepl("^15", ChromPOS))

datamerge2_5<-merge(biochem,molecular_5, by = "Indiv", all = FALSE, sort = FALSE)
datamerge2_15<-merge(biochem,molecular_15, by = "Indiv", all = FALSE, sort = FALSE)


#### Preparing data for analysis(removing sample that have pathogenic, likely pathogenic, and Pseudodeficienct polymorphism) ----

# Selecting patience that have genotype 0/1 for all the variance 
try_0_1<-subset(molecular, (gt_GT == "0/1"|
                                                 gt_GT == "0/2"| 
                                                 gt_GT == "0/3"|
                                                 gt_GT == "0/4"| 
                                                 gt_GT == "0/5"|
                                                 gt_GT == "0/6"|
                                                 gt_GT == "0/7"|
                                                 gt_GT == "0/8"|
                                                 gt_GT == "0/9"|
                                                 gt_GT == "0/10"))

# Merging patience that have genotype 0/1 to the variance  with classification combining these two tables by the locus and  specific getoype in that position 
mergedData_0_1 <- merge(try_0_1, complete_list_of_classified_variants_within_target, by.x=c("ChromPOS","gt_GT_alleles"),
                        by.y=c('ChromPOS', "Genotype"))

# Filtering the table with genotype 0/1 by getting rid of all the patience that have "Pathogenic", "Likely pathogenic" mutations
# sub_special_0_1 <- subset(mergedData_0_1, 
#                           ( Classification == "Pathogenic") |
#                             (Classification == "Likely pathogenic"))

sub_special_0_1 <- subset (mergedData_0_1, 
                          (Classification == "Pathogenic") |
                            (Classification == "Likely pathogenic")|
                               (Classification == "Pseudodeficienct"))

# Now go back and filter the original table based on these unique Indivs that were removed 
molecular_without_path_n1 <- molecular[ ! molecular$Indiv %in% sub_special_0_1$Indiv, ]

# Selecting patience that have genotype 1/1 for all the variance
try_1_1<-subset(molecular, (gt_GT!= "0/0" & 
                                                 gt_GT!= "0/1" & 
                                                 gt_GT!= "0/2" &
                                                 gt_GT != "0/3" & 
                                                 gt_GT != "0/4" & 
                                                 gt_GT != "0/5" &
                                                 gt_GT != "0/6" & 
                                                 gt_GT != "0/7" &
                                                 gt_GT!= "0/8" & 
                                                 gt_GT != "0/9" & 
                                                 gt_GT != "0/10"))
# Seperating the genotype column into ref and alt
new_1_1<-try_1_1 %>% separate(gt_GT_alleles, c('Ref', 'Alt'), sep="/")

# merging patience that have genotype 1/1 for all the variance  with classification of variance and combining these two tables by the chromosome position and by alt genotype
mergedData_1_a <- merge(new_1_1, complete_list_of_classified_variants_within_target, by.x=c("ChromPOS","Alt"),
                        by.y=c('ChromPOS', "ALT" ) )

mergedData_1_b <- merge(new_1_1, complete_list_of_classified_variants_within_target, by.x=c("ChromPOS","Ref"),
                        by.y=c('ChromPOS', "ALT"))
ii<-rbind(mergedData_1_a, mergedData_1_b)
mergedData_1_1<-ii %>% distinct

sub_special_1_1 <- subset(mergedData_1_1, 
                          ( Classification == "Pathogenic") |
                            (Classification == "Likely pathogenic")|
                            (Classification == "Pseudodeficienct"))
                           
# Now go back and filter the original table based on these unique Indivs that were removed 
molecular_without_path <- molecular_without_path_n1[ ! molecular_without_path_n1$Indiv %in% sub_special_1_1$Indiv, ]

molecular_5 <- subset(molecular_without_path, grepl("^5", ChromPOS))
molecular_15 <- subset(molecular_without_path, grepl("^15", ChromPOS))

datamerge2_5<-merge(biochem,molecular_5, by = "Indiv", all = FALSE, sort = FALSE)
datamerge2_15<-merge(biochem,molecular_15, by = "Indiv", all = FALSE, sort = FALSE)

#### Preparing data for analysis(include all samples in analysis)----

molecular_5 <- subset(molecular, grepl("^5", ChromPOS))
molecular_15 <- subset(molecular, grepl("^15", ChromPOS))

datamerge2_5<-merge(biochem,molecular_5, by = "Indiv", all = FALSE, sort = FALSE)
datamerge2_15<-merge(biochem,molecular_15, by = "Indiv", all = FALSE, sort = FALSE)


#### Convert data from long to wide format ----
#5
all1<-dcast(datamerge2_5, Indiv+Result+Sex+Race+year ~ ChromPOS, value.var='gt_GT_alleles', drop=TRUE, fun.aggregate = custom_fun)
all1.2<-dcast(datamerge2_5, Indiv+Result+Sex+Race+year ~ ChromPOS, value.var='gt_GT_alleles', drop=TRUE, fun.aggregate = custom_fun2)

all1.2 <- remove_empty_cols(all1.2)

mergedTable5 <- merge(all1, all1.2,  by ="Indiv",all=TRUE, suffixes = c("","_2"))
mergedTable5<- subset(mergedTable5, select = -c(Result_2, Sex_2, Race_2, year_2))


# mergedTable5<-subset(mergedTable5, Indiv!=1674745)

#____________________

#15
all1_15<-dcast(datamerge2_15, Indiv+Result+Sex+Race+year ~ ChromPOS, value.var='gt_GT_alleles', drop=TRUE, fun.aggregate = custom_fun)
all1_15.2<-dcast(datamerge2_15, Indiv+Result+Sex+Race+year ~ ChromPOS, value.var='gt_GT_alleles', drop=TRUE, fun.aggregate = custom_fun2)


all1_15.2 <- remove_empty_cols(all1_15.2)


mergedTable15 <- merge(all1_15, all1_15.2 ,  by ="Indiv",all=TRUE, suffixes = c("","_2"))
mergedTable15<- subset(mergedTable15, select = -c(Result_2, Sex_2, Race_2, year_2))

# mergedTable15<-subset(mergedTable15, Indiv!=1611839)
# mergedTable15<-subset(mergedTable15, Indiv!=17009280)


####Analysis_15 ----
# Replace NA's with the most frequent value in that column
df1 <- apply(mergedTable15, 2, function(x){ 
    x[is.na(x)] <- names(which.max(table(x)))
    return(x) })
df2_15 <- as.data.frame(df1)
  
pheno_15<-as.numeric(as.character(df2_15$Result))
sex_15<-df2_15$Sex
race_15<-df2_15$Race
year_15<-df2_15$year
  
# convert columns to factor 
try1 = data.frame(sapply(df2_15, as.factor))
# get rid of the columns that have less then one factor
try1<-try1[, sapply(try1, nlevels)  >1]
  
#recode genotype to 0,1,2
# try1 = data.frame(sapply(try1[,6:264], additive))
try1 = data.frame(sapply(try1[,6:224], additive))
#apply function 
# results15<-apply(try1[,1:259], 2, GLM_function15 )
# beta15<-data.frame(apply(try1[,1:259], 2, beta_15 ))
# SE15<-data.frame(apply(try1[,1:259], 2, SE_15 ))

results15<-apply(try1[,1:219], 2, GLM_function15 )
beta15<-data.frame(apply(try1[,1:219], 2, beta_15 ))
SE15<-data.frame(apply(try1[,1:219], 2, SE_15 ))


# FDR adjustment   
pval.corrected.fdr_15 <- p.adjust(results15, method="fdr")
res_15<-pval.corrected.fdr_15
# res_15 <- which(pval.corrected.fdr_15 <= 0.05)
# length(res_15)
res_15<-data.frame(res_15)


#subset variants that are statistically significant
# res_15 <- as.data.frame(subset(pval.corrected.fdr_15, pval.corrected.fdr_15 <=0.05))
# res_15 <- as.data.frame(subset(pval.corrected.fdr_15, pval.corrected.fdr_15 >0.05))

make_row_names_into_a_column_15<-setDT(res_15, keep.rownames = TRUE)[]
beta15<-setDT(beta15, keep.rownames = TRUE)[]
SE15<-setDT(SE15, keep.rownames = TRUE)[]

#merge
make_row_names_into_a_column_15<-merge(SE15, make_row_names_into_a_column_15, by="rn")
make_row_names_into_a_column_15<-merge(beta15, make_row_names_into_a_column_15, by="rn")
###
subset_variant_column_15<-subset(make_row_names_into_a_column_15,select=rn)
ab<-subset_variant_column_15 %>% separate(rn, into=c("3", "2"),sep="X")
significant_variants_15<-(data.table(gsub("[.]", ":", ab$`2`)))$V1

# abd<-res_15%>% separate(rn, into=c("3", "2"),sep="X")
abd<-make_row_names_into_a_column_15%>% separate(rn, into=c("3", "2"),sep="X")
abd$`2`<-data.table(gsub("[.]", ":", abd$`2`))
abr<- subset(abd, select=-c(`3`))
# abr<- subset(abd, select=c(`2`, "subset(pval.corrected.fdr_15, pval.corrected.fdr_15 <= 0.05)"))
# abr<- subset(abd, select=c(`2`, "subset(pval.corrected.fdr_15, pval.corrected.fdr_15 > 0.05)"))
other <- merge(complete_list_of_classified_variants_within_target, abr,  by.x ="ChromPOS", by.y= "2")
other<-data.frame(subset(other, select = -c(REF, ALT)))
other<- rename.vars( other, c("ChromPOS") , c("Variant" ))
df2_15$Result<-as.numeric(as.character(df2_15$Result))

table<- check_integer_mean2_sd(significant_variants_15,df2_15)
# complete_table_15<-merge(table, other, by.x= "Variant", by.y= "ChromPOS")

complete_table_15<-merge(table, other, by.x=c("Variant","Genotype"),
                        by.y=c('Variant', "Genotype"),all.x = TRUE)

export(list(complete_table_15=complete_table_15), "complete_table_15.xlsx")



# res_15_dataframe<-data.frame(res_15)
# res_15_dataframe<-setDT(res_15_dataframe, keep.rownames = TRUE)[]
# abd<-res_15_dataframe%>% separate(rn, into=c("3", "2"),sep="X")
# abd$`2`<-data.table(gsub("[.]", ":", abd$`2`))
# abr<- subset(abd, select=c(`2`, "subset.pval.corrected.fdr_15..pval.corrected.fdr_15....0.05."))
# other <- merge(complete_list_of_classified_variants_within_target, abr,  by.x ="ChromPOS", by.y= "2")
# other<-data.frame(subset(other, select=- c(REF, ALT, Genotype)))

#### Analysis_5 ----
# Replace NA's with the most frequent value in that column
df2 <- apply(mergedTable5, 2, function(x){ 
  x[is.na(x)] <- names(which.max(table(x)))
  return(x) })
df2_5 <- as.data.frame(df2)
  
  
pheno_5<-as.numeric(as.character(df2_5$Result))
sex_5<-df2_5$Sex
race_5<-df2_5$Race
year_5<-df2_5$year
  
# convert columns to factor 
try1_5 = data.frame(sapply(df2_5, as.factor))
# get rid of the columns that have less then one factor
try1_5<-try1_5[, sapply(try1_5, nlevels)  >1]
  
#recode genotype to 0,1,2
# try1_5 = data.frame(sapply(try1_5[,6:308], additive))
try1_5 = data.frame(sapply(try1_5[,6:276], additive))
#apply function 
# b_5<-apply(try1_5[,1:303], 2, GLM_function5)
# beta5<-data.frame(apply(try1_5[,1:303], 2, beta_5 ))
# SE5<-data.frame(apply(try1_5[,1:303], 2, SE_5 ))

b_5<-apply(try1_5[,1:271], 2, GLM_function5)
beta5<-data.frame(apply(try1_5[,1:271], 2, beta_5 ))
SE5<-data.frame(apply(try1_5[,1:271], 2, SE_5 ))

# FDR adjustment   
pval.corrected.fdr_5 <- p.adjust(b_5, method="fdr")
res_5<-as.data.frame(pval.corrected.fdr_5)

#subset variants that are statistically significant
 # res_5 <- as.data.frame(subset(pval.corrected.fdr_5, pval.corrected.fdr_5 <=0.05))
# res_5 <- as.data.frame(subset(pval.corrected.fdr_5, pval.corrected.fdr_5 >0.05))
make_row_names_into_a_column_5<-setDT(res_5, keep.rownames = TRUE)[]
beta5<-setDT(beta5, keep.rownames = TRUE)[]
SE5<-setDT(SE5, keep.rownames = TRUE)[]

# merge
make_row_names_into_a_column_5<-merge(SE5, make_row_names_into_a_column_5, by="rn")
make_row_names_into_a_column_5<-merge(beta5, make_row_names_into_a_column_5, by="rn")

subset_variant_column_5<-subset(make_row_names_into_a_column_5,select=rn)
ab<-subset_variant_column_5 %>% separate(rn, into=c("3", "2"),sep="X")
significant_variants_5<-(data.table(gsub("[.]", ":", ab$`2`)))$V1

# res_5_dataframe<-data.frame(res_5)
# res_5_dataframe<-setDT(res_5_dataframe, keep.rownames = TRUE)[]
# ab<-res_5%>% separate(rn, into=c("3", "2"),sep="X")
ab<-make_row_names_into_a_column_5%>% separate(rn, into=c("3", "2"),sep="X")
ab$`2`<-data.table(gsub("[.]", ":", ab$`2`))
ar<- subset(ab, select=-c(`3`))
# ar<- subset(ab, select=c(`2`, "subset(pval.corrected.fdr_5, pval.corrected.fdr_5 > 0.05)"))
other_2<- merge(complete_list_of_classified_variants_within_target, ar,  by.x ="ChromPOS", by.y= "2")
# other_2<-data.frame(subset(other_2, select=- c(REF, ALT, Genotype)))
other_2<-data.frame(subset(other_2, select=- c(REF, ALT)))
other_2<- rename.vars( other_2, c("ChromPOS") , c("Variant" ))
df2_5$Result<-as.numeric(as.character(df2_5$Result))

table3<- check_integer_mean2_sd(significant_variants_5,df2_5)
# complete_table_5<-merge(table3, other_2, by.x="Variant", by.y= "Genotype")
complete_table_5<-merge(table3, other_2, by.x=c("Variant","Genotype"),
                        by.y=c('Variant', "Genotype"),all.x = TRUE)
                        
export(list(complete_table_5=complete_table_5), "complete_table_5.xlsx")

#### LD ----
asd1 <- try1_5[,grep("(X15|X5)", colnames(try1_5))]
mdf <- NULL
  for(i in names(asd1)) { 
  something <- genotype(asd1[[i]]) # the current computed genotype
  mdf[[i]] <- something # add to it to list of thingies
}
ldall1 <- genetics::LD(data.frame(mdf))
a<-ldall1$"R^2"
  
b<-a[rowSums(is.na(a))!=ncol(a), ]
  
export(list(b=b), "LD.xlsx")
  
  
asd1 <- try1[,grep("(X15|X5)", colnames(try1))]
mdf <- NULL
  for(i in names(asd1)) { 
    something <- genotype(asd1[[i]]) # the current computed genotype
    mdf[[i]] <- something # add to it to list of thingies
}
ldall1 <- genetics::LD(data.frame(mdf))
a<-ldall1$"R^2"
  
c<-a[rowSums(is.na(a))!=ncol(a), ]
  
export(list(c=c), "LD5.xlsx")

####check----
t1 <- subset(df2_5, select=c("5:73981144"))
table(t1$`5:73981144`)

5:74009338_2

t1 <- subset(df2_5, (ChromPOS == "5:74009338" & gt_GT !="0/0"))
table(t1$gt_GT_alleles)
View((t1))

t2 <- subset(molecular1, (POS == "72668305" & gt_GT !="0/0" ))
table(t2$gt_GT_alleles)
View((t2))

t3 <- subset(molecular2, (POS == "72668305" & gt_GT !="0/0"  ))
table(t3$gt_GT_alleles)
View((t3))

t1 <- subset(df2_5, select=c("5:74009338"))
table(t1$gt_GT_alleles)


#### boxplot ----



ww<-dcast(datamerge2_15, Indiv+Result ~ ChromPOS, value.var='gt_GT_alleles', drop=TRUE, fun.aggregate = custom_fun)                                 
ww$Result<-as.numeric(ww$Result)
ww<- rename.vars( ww, c("Result") , c("Enzyme_Activity" ))

rects <- data.frame(ystart = c(30.5,50,55,75), yend= c(50,55,75,100), col = c( "Tay-Sachs Carrier","Inconclusive","Non-Carrier","Sandhoff-Carrier"))

pdf("boxplot1.pdf")
ggplot(ww) +
geom_boxplot(aes(x = `15:72637795`, y = Enzyme_Activity, fill = `15:72637795`)) +
annotate("text", x = 1, y = 30.5, label = "P-value = 0") +
scale_fill_brewer(palette="Dark2")+
scale_y_continuous(breaks=c(50,55,75))+
geom_rect(data = rects, aes(xmin = -Inf, xmax = Inf ,ymin = ystart, ymax = yend, fill = col), alpha = 0.4)+
# geom_jitter(width=0.07,aes(x = `15:72637795`, y = Enzyme_Activity))+
# ggtitle('Boxplot of genome size by citrate mutant type') +
# xlab('citrate mutant') +
# ylab('genome size') +
theme(panel.grid.major = element_line(size = .7, color = "grey"),
  axis.text.x = element_text(angle=45, hjust=1),
  axis.title = element_text(size = rel(1.5)),
  axis.text = element_text(size = rel(1.25)), legend.position="F")
  # axis.text = element_text(size = rel(1.25)))
dev.off()
                 
      

pdf("15_72638892_smaller_dots_boxplot2.pdf")
ggplot(ww) +
  geom_boxplot(aes(x = `15:72638892`, y = Enzyme_Activity, fill = `15:72638892`), show.legend = F) +
  annotate("text", x = 1, y = 30.5, label = "P-value = 0") +
  scale_fill_brewer(palette="Dark2")+
  scale_y_continuous(breaks=c(50,55,75))+
  # geom_rect(data = rects, aes(xmin = -Inf, xmax = Inf ,ymin = ystart, ymax = yend, fill = col), alpha = 0.4)+
  # geom_jitter(width=0.07,aes(x = `15:72638892`, y = Enzyme_Activity))+
  # ggtitle('Boxplot of genome size by citrate mutant type') +
  # xlab('citrate mutant') +
  # ylab('genome size') +
  theme(panel.grid.major = element_line(size = .7, color = "grey"),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_text(size = rel(1.5)),
        # axis.text = element_text(size = rel(1.25)), legend.position="F")
axis.text = element_text(size = rel(1.25)))
dev.off()    
# make dot smalles? done              
# add cutoff range for carrier, noncarrier                    
# make boxplot without pathogenic and pseudofivient                 
# add p value 
               
# figure legend()

                              
                 
                 
                 
                 
                 
