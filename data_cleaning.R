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


#### Molecular data;  ----

# HEXA gene is on chromosome 15
# HEXB gene is on chromosome 5

# Importing vcf  file
data1 <- read.vcfR( "/Users/olgalukatskaya/Desktop/work_projects/SR-2738/DataFiles/work/molecular-results_backup_modded.vcf.gz", verbose = FALSE )
data2 <- read.vcfR( "/Users/olgalukatskaya/Desktop/work_projects/SR-2738/DataFiles/merged_large_selected_second_modded.vcf.gz", verbose = FALSE )


# Cleaning the molecular data and converting it to long format
# 1st mole data
longdata<- vcfR2tidy(data1, toss_INFO_column = FALSE,format_fields = c("GT"))
molecular1<-longdata$gt
# 2nd mole data
longdata2<- vcfR2tidy(data2, toss_INFO_column = FALSE,format_fields = c("GT"))
molecular2<-longdata2$gt


# Combining two files together
molecular<-rbind(molecular1, molecular2)


# Filtering: getting rid of characters before and after underscore
molecular$Indiv <- gsub("_.{2}$", "", molecular$Indiv)
molecular$Indiv <- gsub("(^.{3,5}_)", "", molecular$Indiv)


# Renaming chromosomes
molecular$ChromKey<-replace(molecular$ChromKey, molecular$ChromKey==2, 5)
molecular$ChromKey<-replace(molecular$ChromKey, molecular$ChromKey==1, 15)


# Merging chromkey and pos column to create a new column ChromPOS and getting rid of Chromkey and Pos
molecular<-within(molecular, ChromPOS <- paste(ChromKey,POS,sep=':'))
molecular<- subset(molecular, select = -c(ChromKey,POS) ) 

# removing patience 
molecular<-subset(molecular, Indiv!= "1551180" & Indiv!="1604474"& Indiv!="1614437"& Indiv!="1616165"&Indiv!=" 1646388"&
                    Indiv!="1653689")

t1 <- subset(molecular, (Indiv == "1604474" & gt_GT !="0/0")) 

# Export the results
save(molecular,file="data.molecular")

#### Removing non-target varaints from molecular data ---- 

#editing the table with the hexa/b targe regions
target_regions<-read.xlsx("/Users/olgalukatskaya/Desktop/work_projects/SR-2738/Jun_analysis_details /HEXA_HEXB_Target_Regions.xlsx",colNames=FALSE,1)
target_regions<-target_regions %>% separate(X1, c('h', 'Chrom'), sep="r")
target_regions<- subset(target_regions, select=-c(h))

#separating the molecular table and target region  table based on 2 chromosomes 
regions_5 <- subset(target_regions, grepl("^5", Chrom))
regions_15 <- subset(target_regions, grepl("^15", Chrom))
molecular_chrompos_seperated<-molecular %>% separate(ChromPOS, c('Chrom', 'POS'), sep=":")
molecular_15 <- subset(molecular_chrompos_seperated, grepl("^15", Chrom))
molecular_5 <- subset(molecular_chrompos_seperated, grepl("^5", Chrom))

# getting all the values that are not in the target regions
#listing all the target sequences from start(X2) to end(X3) and stuffing them into a vector and  excluding regions from the 
# molecular table so that now we are left with non target regions that we want to remove
#15
non_target_reigons_15<-as.data.frame(setdiff(molecular_15$POS,unlist(with(regions_15, Map( `:`, X2, X3)))))
# rename vars 
non_target_reigons_15<- rename.vars( non_target_reigons_15, from="setdiff(molecular_15$POS, unlist(with(regions_15, Map(`:`, X2, X3))))", to = "var1")
# removing non target reigons from the molecular table   
molecular_15_target_regions <- molecular_15[ ! molecular_15$POS %in% non_target_reigons_15$var1, ]

#5
non_target_reigons_5<-as.data.frame(setdiff(molecular_5$POS,unlist(with(regions_5, Map( `:`, X2, X3)))))
# rename vars 
non_target_reigons_5<- rename.vars( non_target_reigons_5, from="setdiff(molecular_5$POS, unlist(with(regions_5, Map(`:`, X2, X3))))", to = "var1")
# removing non target reigons from the molecular table   
molecular_5_target_regions <- molecular_5[ ! molecular_5$POS %in% non_target_reigons_5$var1, ]

molecular_within_target_region <- rbind(molecular_15_target_regions, molecular_5_target_regions)
# molecular_within_target_region<- molecular_within_target_region %>% unite(ChromPOS, Chrom,POS, sep=":", remove=TRUE)
molecular<- molecular_within_target_region %>% unite(ChromPOS, Chrom,POS, sep=":", remove=TRUE)

#### Biochemical data ----

# Importing enzyme file
biochem1 <- readr::read_tsv("/Users/olgalukatskaya/Desktop/work_projects/SR-2738/DataFiles/biochemical-results-oldtx.txt", na=".")
biochem2 <- readr::read_tsv("/Users/olgalukatskaya/Desktop/work_projects/SR-2738/DataFiles/biochemical-results-newtx.txt", na=".")

#delete a row in biochem1
biochem1<-subset(biochem1, `Lab #_1`!="1634147NG"  )
#delete  rows in biochem2
biochem2<-subset(biochem2, `Lab #_1`!="EPIC" |  `System #`!="P0388233"   )
biochem2<-subset(biochem2, `System #`!="P0239572"   )

# Combining the two tables into one
biochem <- rbind(biochem1, biochem2)

# Filtering: in the "Lab #_1" column, only keep numeric characters; 
biochem$`Lab #_1` <- gsub("[^0-9]", "", biochem$`Lab #_1`) 

# getting rid of duplicates - retain only the most recent date for every ID
biochem <- 
  biochem %>%
  group_by(`System #`) %>%
  slice(which.max(as.Date(Collected, '%m/%d/%Y')))

# get rid of duplicate samples that has different system # 
biochem<-subset(biochem, `System #`!="P0349522")

# Rename Lab #_1 to Indiv and get rid of  Patient and Lab #_1 columns 
biochem$Indiv<-biochem$`Lab #_1`

# change the results value from 0 or from 178 to thier true enzyme activity 
biochem <- within(biochem, Result[Result == '0' & Indiv == '1650769'] <- '60.3')
biochem <- within(biochem, Result[Result == '0' & Indiv == '17040667'] <- '82')
biochem <- within(biochem, Result[Result == '0' & Indiv == '17040875'] <- '84.3')
biochem <- within(biochem, Result[Result == '178' & Indiv == '1659152'] <- '69.6')

# use data.table reference to re-categorize 'Result' into new categories and save that as 'biochem_groups' column
# No need to use the biochem_table after this, can go back to using the biochem data.frame
# biochem_table <- setDT(biochem)
# biochem_table[Result < 50.0, biochem_groups := "Tay-Sachs Carrier"]
# biochem_table[Result >= 50.0 & Result <=55.0, biochem_groups := "Inconclusive"]
# biochem_table[Result > 55.0 & Result <=72.0, biochem_groups := "Non-Carrier"]
# biochem_table[Result >72.0, biochem_groups:= "Sandhoff-Carrier"]

###
biochem_table <- setDT(biochem)
biochem_table[Result < 50.0, biochem_groups := "Tay-Sachs Carrier"]
biochem_table[Result >= 50.0 & Result <=55.0, biochem_groups := "Inconclusive"]
biochem_table[Result > 55.0 & Result <=75.0, biochem_groups := "Non-Carrier"]
biochem_table[Result >75.0, biochem_groups:= "Sandhoff-Carrier"]

# biochem<-subset(biochem, select = c(Indiv,biochem_groups,Result) )
biochem<-subset(biochem, select = c(Indiv,biochem_groups,Result,Collected,DOB,Sex, Race, Ethnicity) )


#### Make additional changes to the biochem data file ----

library(qdap)
race<-as.data.frame(table(biochem$Race))
race$Var1 <-  genX(race$Var1, " (", ")")
a<-data.frame(do.call(rbind, strsplit(race$"Var1", split="(?<=[a-z])(?=[A-Z])", perl=T)))
X1<-data.frame(table(a$X1))
X2<-data.frame(table(a$X2))
X3<-data.frame(table(a$X3))
X4<-data.frame(table(a$X4))
X5<-data.frame(table(a$X5))
X6<-data.frame(table(a$X6))
X7<-data.frame(table(a$X7))
X8<-data.frame(table(a$X8))
X9<-data.frame(table(a$X9))
x10<-rbind(X1,X2,X3,X4,X5,X6,X7,X8,X9)
b<-setDT(x10)[, list(Freq=paste(Freq, collapse=";")) ,c("Var1")]
__________

caucasian<-list("European","Caucasian","White","Albanian","Italian","German","Russian",
                "Armenian","Austrian","Bosnian","Belarusian","Belgian","Croatian",
                "Cypriot", "Czech","Danish","Dutch","Finnish","Georgian","Lithuanian",
                "Southern European","Maltese","Portuguese","Scandanvian","Serbian",
                "Tatar","Welsh", "Eastern European","English","French","Polish",
                "Western European", "Northern European","Ukranian","Yugoslavian", "Irish","Hungarian","Greek",
                "Swiss","Swedish","Spaniard", "Slovakian","Scottish","Norwegian","Bulgarian", "Roman","Romanian",
                "Latvian", "French Canadian", "French Canadian or Cajun","Macedonian","Sicilian", "Australian")

ashkenazi_jewish<-list("Ashkenazi Jewish")
asian<-list("Asian","Chinese","Korean","Japanese","Vietnamese", "South Asian","Southeast Asian","Phillipino","Thai","Nepali",
            "Pakistani","Indian","East Asian","Uzbekistan","Bengali")

african_american<-list("African","African American","North African","South African Afrikaner","Tunisian", "Black", "Jamaican", "Haitian")

middle_eastern<-list("Arabic","Egyptian","Middle Eastern", "Palestinian", 
                     "Persian","Israeli", "Iraqi","Iranian", "Sephardic Jewish",
                     "Sephardic","Qatari", "Syrian","Turkish","Yemenis","Lebanese","Druze Northern Israeli", "Sephardic Jewish")

# other_mixed <-list("Canadian Inuit","Native American", "Pacific Islander",
#                    "Lumbee Native Americans", "Trinidadian", "^Jewish$", "Other/ Mixed Caucasian")

unknown<-list("Unknown","Unknown- Adopted","Mediterranean",
              "Argentinian","Canary Islands","Caribbean","Guyanese","^Non-Jewish$","AS","TH","WH", "Canadian Inuit","Native American", "Pacific Islander",
              "Lumbee Native Americans", "Trinidadian", "^Jewish$", "Other/ Mixed Caucasian")

hispanic<-list("Brazilian","Colombian","Cuban","Dominican","Ecuadorian","Hispanic","Hispanic/Latino","Mexican","Puerto Rican")

t2 <- paste(unknown, collapse="|")
t3 <- paste(hispanic, collapse="|")
t4 <- paste(caucasian, collapse="|")
t5 <- paste(asian, collapse="|")
t6 <- paste(african_american, collapse="|")
t7 <- paste(middle_eastern, collapse="|")
# t8 <- paste(other_mixed, collapse="|")

____________

# biochem_copy <-biochem

biochem<-biochem_copy
# if race is cell is empty replace it with unknown
biochem$Race[biochem$Race==""] <- "Unknown"
# if ethnicity cell is empty replace it with unknown
biochem$Ethnicity[biochem$Ethnicity==""] <- "Unknown"
# if the cell race  is is.na ==TRUE  then replace it with unknown 
biochem$Race[is.na(biochem$Race)] <- "Unknown"
# if the cell ethnicity  is is.na ==TRUE  then replace it with unknown
biochem$Ethnicity[is.na(biochem$Ethnicity)] <- "Unknown"
# change Not Provided to unknown in race and ethnicity column
biochem <- within(biochem, Race[ grepl("^Not Provided$",Race) ] <- 'Unknown')
biochem <- within(biochem, Ethnicity[ grepl("^Not Provided$",Ethnicity) ] <- 'Unknown')
# # if ethnicity cell has a word ashkenazi then replace it with ashkenazi jewish and make the race ashkenazi jewish
biochem <- within(biochem, Ethnicity[ grepl("Ashkenazi",Ethnicity) ] <- 'Ashkenazi Jewish')
biochem[ Ethnicity == "Ashkenazi Jewish", Race := Ethnicity ]

# biochem$Race <- ifelse(biochem$Ethnicity == "Ashkenazi", Race := Ethnicity,  biochem$Race)

biochem$Race <- ifelse(biochem$Race == "Italian" & biochem$Ethnicity!="Unknown", biochem$Ethnicity, biochem$Race)
biochem$Race <- ifelse(biochem$Race == "Unknown", biochem$Ethnicity,  biochem$Race)



biochem$var2 <- ifelse(grepl(t2,biochem$Race),'unknown','')
biochem$var3 <- ifelse( grepl(t3,biochem$Race),'hispanic', '')
biochem$var4 <- ifelse( grepl(t4,biochem$Race),'caucasian','')
biochem$var5 <- ifelse( grepl(t5,biochem$Race),'asian', '')
biochem$var6 <- ifelse( grepl(t6,biochem$Race),'african_american', '')
biochem$var7 <- ifelse( grepl(t7,biochem$Race),'middle_eastern', '')
# biochem$var8 <- ifelse( grepl(t8,biochem$Race),'other_mixed', '')
biochem$var10 <- ifelse( grepl("Ashkenazi",biochem$Race),'ashkenazi', '')

biochem$race <- ifelse(biochem$var2 == "unknown" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "" , 'unknown', '')
biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "hispanic" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "" , 'hispanic', biochem$race)
biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "caucasian" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "" , 'caucasian', biochem$race)
biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "asian" & biochem$var6 == "" & biochem$var7 == ""  , 'asian', biochem$race)
biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "african_american" & biochem$var7 == "" , 'african_american', biochem$race)
biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "middle_eastern" , 'middle_eastern', biochem$race)
# biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "" & biochem$var8 == "other_mixed", 'other_mixed', biochem$race)


# only ashkenazi
biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == ""  & biochem$var10 == "ashkenazi", 'ashkenazi', biochem$race)
biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "caucasian" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == ""  & biochem$var10 == "ashkenazi", 'ashkenazi', biochem$race)

# # ashkenazi mix 
# biochem$race <- ifelse(biochem$var2 == "unknown" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "" & biochem$var8 == "" & biochem$var10 == "ashkenazi", 'ashkenazi_mix', biochem$race)
# biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "hispanic" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "" & biochem$var8 == "" & biochem$var10 == "ashkenazi", 'ashkenazi_mix', biochem$race)
# biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "asian" & biochem$var6 == "" & biochem$var7 == "" & biochem$var8 == "" & biochem$var10 == "ashkenazi", 'ashkenazi_mix', biochem$race)
# biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "african_american" & biochem$var7 == "" & biochem$var8 == "" & biochem$var10 == "ashkenazi", 'ashkenazi_mix', biochem$race)
# biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "middle_eastern" & biochem$var8 == "" & biochem$var10 == "ashkenazi", 'ashkenazi_mix', biochem$race)
# biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "" & biochem$var8 == "other_mixed" & biochem$var10 == "ashkenazi", 'ashkenazi_mix', biochem$race)
# biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "caucasian" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "" & biochem$var8 == "other_mixed" & biochem$var10 == "ashkenazi", 'ashkenazi_mix', biochem$race)
# biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "caucasian" & biochem$var5 == "" & biochem$var6 == "african_american" & biochem$var7 == "" & biochem$var8 == "" & biochem$var10 == "ashkenazi", 'ashkenazi_mix', biochem$race)
# biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "caucasian" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "middle_eastern" & biochem$var8 == "" & biochem$var10 == "ashkenazi", 'ashkenazi_mix', biochem$race)
# biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "hispanic" & biochem$var4 == "caucasian" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "" & biochem$var8 == "" & biochem$var10 == "ashkenazi", 'ashkenazi_mix', biochem$race)
# biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "hispanic" & biochem$var4 == "" & biochem$var5 == "asian" & biochem$var6 == "" & biochem$var7 == "" & biochem$var8 == "" & biochem$var10 == "ashkenazi", 'ashkenazi_mix', biochem$race)
# biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "caucasian" & biochem$var5 == "asian" & biochem$var6 == "" & biochem$var7 == "" & biochem$var8 == "" & biochem$var10 == "ashkenazi", 'ashkenazi_mix', biochem$race)
# biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "african_american" & biochem$var7 == "middle_eastern" & biochem$var8 == "" & biochem$var10 == "ashkenazi", 'ashkenazi_mix', biochem$race)
# biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "african_american" & biochem$var7 == "" & biochem$var8 == "other_mixed" & biochem$var10 == "ashkenazi", 'ashkenazi_mix', biochem$race)
# biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "hispanic" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "middle_eastern" & biochem$var8 == "" & biochem$var10 == "ashkenazi", 'ashkenazi_mix', biochem$race)
# biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "hispanic" & biochem$var4 == "caucasian" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "middle_eastern" & biochem$var8 == "" & biochem$var10 == "ashkenazi", 'ashkenazi_mix', biochem$race)



# ashkenazi mix
biochem$race <- ifelse(biochem$var2 == "unknown" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "" &  biochem$var10 == "ashkenazi", 'ashkenazi', biochem$race)
biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "hispanic" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "" & biochem$var10 == "ashkenazi", 'ashkenazi', biochem$race)
biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "asian" & biochem$var6 == "" & biochem$var7 == "" & biochem$var10 == "ashkenazi", 'ashkenazi', biochem$race)
biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "african_american" & biochem$var7 == ""  & biochem$var10 == "ashkenazi", 'ashkenazi', biochem$race)
biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "middle_eastern"  & biochem$var10 == "ashkenazi", 'ashkenazi', biochem$race)
# biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "" & biochem$var8 == "other_mixed" & biochem$var10 == "ashkenazi", 'ashkenazi', biochem$race)
# biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "caucasian" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "" & biochem$var8 == "other_mixed" & biochem$var10 == "ashkenazi", 'ashkenazi', biochem$race)
biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "caucasian" & biochem$var5 == "" & biochem$var6 == "african_american" & biochem$var7 == "" &  biochem$var10 == "ashkenazi", 'ashkenazi', biochem$race)
biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "caucasian" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "middle_eastern"  & biochem$var10 == "ashkenazi", 'ashkenazi', biochem$race)
biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "hispanic" & biochem$var4 == "caucasian" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == ""  & biochem$var10 == "ashkenazi", 'ashkenazi', biochem$race)
biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "hispanic" & biochem$var4 == "" & biochem$var5 == "asian" & biochem$var6 == "" & biochem$var7 == ""  & biochem$var10 == "ashkenazi", 'ashkenazi', biochem$race)
biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "caucasian" & biochem$var5 == "asian" & biochem$var6 == "" & biochem$var7 == ""  & biochem$var10 == "ashkenazi", 'ashkenazi', biochem$race)
biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "african_american" & biochem$var7 == "middle_eastern" & biochem$var10 == "ashkenazi", 'ashkenazi', biochem$race)
# biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "african_american" & biochem$var7 == "" & biochem$var8 == "other_mixed" & biochem$var10 == "ashkenazi", 'ashkenazi', biochem$race)
biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "hispanic" & biochem$var4 == "" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "middle_eastern" &  biochem$var10 == "ashkenazi", 'ashkenazi', biochem$race)
biochem$race <- ifelse(biochem$var2 == "" & biochem$var3 == "hispanic" & biochem$var4 == "caucasian" & biochem$var5 == "" & biochem$var6 == "" & biochem$var7 == "middle_eastern"  & biochem$var10 == "ashkenazi", 'ashkenazi', biochem$race)


# mixed <- subset(biochem, (race == "other_mixed" ))
# Mixed<-as.data.table(table(mixed$Race))
# write.csv(Mixed, file="mixed_race.csv")
# 
# ashkenazi_mix <- subset(biochem, (race == "ashkenazi_mix" ))
# ashkenazi_Mix<-as.data.table(table(ashkenazi_mix$Race))
# write.csv(ashkenazi_Mix, file="ashkenazi_mixed_race.csv")



# biochem$race <- ifelse(biochem$Race == "Ashkenazi Jewish" , 'Ashkenazi Jewish',biochem$race)
# biochem$race <- ifelse(biochem$Race == "Other/ Mixed Caucasian" , 'other_mixed',biochem$race)
# biochem$race <- ifelse(biochem$Race == "West Indian" , 'other_mixed',biochem$race)
# biochem$race <- ifelse(biochem$Race == "American Indian/Alaskan" , 'other_mixed',biochem$race)
# biochem$race <- ifelse(biochem$race == "" , 'other_mixed',biochem$race)

biochem$race <- ifelse(biochem$Race == "West Indian" , 'unknown',biochem$race)
biochem$race <- ifelse(biochem$Race == "American Indian/Alaskan" , 'unknown',biochem$race)
biochem$race <- ifelse(biochem$race == "" , 'unknown',biochem$race)

biochem$Result<-as.numeric(biochem$Result)
mean15 <- aggregate( Result ~  `race`, new, mean)

###
race1<-as.data.frame(table(biochem$race))
table(biochem$race)


#### Enzyme results ----

#Average Hex A%
biochem$Result<-as.numeric(biochem$Result)
mean15 <- aggregate( Result ~  `race`, biochem, mean)

# Tay-Sachs Carrier
tay_sachs_carrier_race <- subset(biochem, (biochem_groups == "Tay-Sachs Carrier" ))
tay_sachs_carrier_race1<-as.data.frame(table(tay_sachs_carrier_race$race))

# Inconclusive
inconclusive_race <- subset(biochem, (biochem_groups == "Inconclusive" ))
inconclusive_race_1<-as.data.frame(table(inconclusive_race $race))


# Normal(non-carrier)
non_carrier_race <- subset(biochem, (biochem_groups == "Non-Carrier" ))
non_carrier_race_1<-as.data.frame(table(non_carrier_race$race))

# Sandhoff Carrier
sandhoff_carrier_race <- subset(biochem, (biochem_groups == "Sandhoff-Carrier" ))
sandhoff_carrier_race_1<-as.data.frame(table(sandhoff_carrier_race $race))


#### NGS data summary table ----

merged_data_race= merge(biochem, combined_tables, by='Indiv',all = FALSE, sort = FALSE)
merged_data_race<-merged_data_race%>% separate(ChromPOS, c('Chrom', 'POS'), sep=":")
merged_data_race$Classification <- as.character(merged_data_race$Classification)

reclassified_vars <- within(merged_data_race, Classification[Classification == 'Likely pathogenic'] <- 'Pathogenic')
reclassified_vars <- within(reclassified_vars, Classification[Classification == 'pathogenic'] <- 'Pathogenic')
reclassified_vars <- within(reclassified_vars, Classification[Classification == 'Likely benign'] <- 'Benign')
reclassified_vars <- within(reclassified_vars, Classification[Classification == 'Likely Benign'] <- 'Benign')

# No HEXA/HEXB variants or only benign/likely benign variants found
not_benign<-subset(reclassified_vars, Classification=="Pathogenic"|Classification=="Pseudodeficienct"|Classification=="VOUS" )
# samples that are in non benign are excluded and we are left with samples that are only have benigh mutations
benign <- reclassified_vars[ ! reclassified_vars$Indiv %in% not_benign$Indiv, ]
kk<-subset(benign,   select=c(Indiv, race, Classification))
# unique the samples that are benign
only_benign<-kk[!duplicated(kk),]
table(only_benign$race)

# figure out which samples are duplicated in the non benign dataframe and these are the ones
# that have more then one mutation and exclude those from the list
bb<-subset(not_benign,   select=c(Indiv))
more_then_one_mutation<-bb[duplicated(bb),]


#exclude sample that have more then one mutation. these samples only 1 specific classification 
not_benign_with_one_mutation <- not_benign[ ! not_benign$Indiv %in% more_then_one_mutation$Indiv, ]

#15
Pathogenic_15<-subset(not_benign_with_one_mutation, Classification=="Pathogenic"  & Chrom=="15")
table(Pathogenic_15$race)


Pseudodeficienct_15<-subset(not_benign_with_one_mutation, Classification=="Pseudodeficienct"  & Chrom=="15")
table(Pseudodeficienct_15$race)

VOUS_15<-subset( not_benign_with_one_mutation, Classification=="VOUS"  & Chrom=="15")
table(VOUS_15$race)

#5
Pathogenic_5<-subset(not_benign_with_one_mutation, Classification=="Pathogenic"  & Chrom=="5")
table(Pathogenic_5$race)

Pseudodeficienct_5<-subset(not_benign_with_one_mutation, Classification=="Pseudodeficienct"  & Chrom=="5")
table(Pseudodeficienct_5$race)

VOUS_5<-subset(not_benign_with_one_mutation, Classification=="VOUS"  & Chrom=="5")
table(VOUS_5$race)

# samples that have at least 2 mutations 
multiple_mutations <- not_benign[not_benign$Indiv %in% more_then_one_mutation$Indiv, ]

ww<-subset(multiple_mutations,   select=c(Indiv, race))
ff<-ww[!duplicated(ww),]
table(ff$race)

aggregate(Result, FUN=mean, data=biochem)
mean(biochem$Result)


#### Summary table for each ethnicity ---- 
molecular_without_path

merged_data= merge(biochem, combined_tables, by='Indiv',all = FALSE, sort = FALSE)
merged_data<-merged_data%>% separate(ChromPOS, c('Chrom', 'POS'), sep=":")
merged_data$Classification <- as.character(merged_data$Classification)

reclassified_vars <- within(merged_data, Classification[Classification == 'Likely pathogenic'] <- 'Pathogenic')
reclassified_vars <- within(reclassified_vars, Classification[Classification == 'pathogenic'] <- 'Pathogenic')
reclassified_vars <- within(reclassified_vars, Classification[Classification == 'Likely benign'] <- 'Benign')
reclassified_vars <- within(reclassified_vars, Classification[Classification == 'Likely Benign'] <- 'Benign')


# race_string <- c("african_american" ,"ashkenazi","asian", "caucasian","hispanic","middle_eastern","unknown")
race_string <- ("caucasian")
# race_list<- vector("list", 70)
for ( i in race_string) {
  # x<-subset(reclassified_vars , grepl( i, race ))
  # x = subset(reclassified_vars, race == "asian")
  x = subset(reclassified_vars, race == i)
  
  not_benign<-subset(x, Classification=="Pathogenic"|Classification=="Pseudodeficienct"|Classification=="VOUS" )
  benign <- x[ !x$Indiv %in% not_benign$Indiv, ]
  kk<-subset(benign,   select=c(Indiv, biochem_groups, Classification))
  only_benign<-kk[!duplicated(kk),]
  #benign
  s<-as.data.frame(table(only_benign$biochem_groups))
  s<- rename.vars( s, c("Freq") , c("benign" ))
  print(s)
  
  
  
  # race_list[[i]] <- table(only_benign_cau$biochem_groups)
  bb<-subset(not_benign,   select=c(Indiv))
  more_then_one_mutation<-bb[duplicated(bb),]
  
  
  
  #exclude sample that have more then one mutation. these samples only 1 specific classification 
  not_benign_with_one_mutation <- not_benign[ ! not_benign$Indiv %in% more_then_one_mutation$Indiv, ]
  #15 # HEXA
  
  Pathogenic_15<-subset(not_benign_with_one_mutation, Classification=="Pathogenic"  & Chrom=="15")
  a<-as.data.frame(table(Pathogenic_15$biochem_groups))
  a<- rename.vars( a, c("Freq") , c("Pathogenic HEXA" ))
  print(a)
  
  
  
  Pseudodeficienct_15<-subset(not_benign_with_one_mutation, Classification=="Pseudodeficienct"  & Chrom=="15")
  b<-as.data.frame(table(Pseudodeficienct_15$biochem_groups))
  b<- rename.vars( b, c("Freq") , c("Pseudodeficienct HEXA" ))
  print(b)
  
  
  VOUS_15<-subset( not_benign_with_one_mutation, Classification=="VOUS"  & Chrom=="15")
  c<-as.data.frame(table(VOUS_15$biochem_groups))
  c<- rename.vars( c, c("Freq") , c(" VOUS HEXA" ))
  print(c)
  
  #5   # HEXB
  Pathogenic_5<-subset(not_benign_with_one_mutation, Classification=="Pathogenic"  & Chrom=="5")
  d<-as.data.frame(table(Pathogenic_5$biochem_groups))
  d<- rename.vars( d, c("Freq") , c(" Pathogenic HEXB" ))
  print(d)
  
  
  Pseudodeficienct_5<-subset(not_benign_with_one_mutation, Classification=="Pseudodeficienct"  & Chrom=="5")
  e<-as.data.frame(table(Pseudodeficienct_5$biochem_groups))
  e<- rename.vars( e, c("Freq") , c(" Pseudodeficienct HEXB" ))
  print(e)
  
  
  VOUS_5<-subset(not_benign_with_one_mutation, Classification=="VOUS"  & Chrom=="5")
  f<-as.data.frame(table(VOUS_5$biochem_groups))
  f<- rename.vars( f, c("Freq") , c(" VOUS HEXB" ))
  print(f)
  
  # samples that have at least 2 mutations 
  multiple_mutations <- not_benign[not_benign$Indiv %in% more_then_one_mutation$Indiv, ]
  
  multiple_mutations_HEXA<-subset(multiple_mutations, grepl("^15", Chrom))
  multiple_mutations_HEXB<-subset(multiple_mutations, grepl("^5", Chrom))
  
  # samples that have one hex a and hexb mutation 
  common_Indiv<-data.frame(intersect(multiple_mutations_HEXA$Indiv,multiple_mutations_HEXB$Indiv))
  mutations_in_HEX_A_and_B <- multiple_mutations [multiple_mutations$Indiv %in% common_Indiv$intersect.multiple_mutations_HEXA.Indiv..multiple_mutations_HEXB.Indiv., ]
  
  ww<-subset(mutations_in_HEX_A_and_B,   select=c(Indiv, biochem_groups))
  ff<-ww[!duplicated(ww),]
  g<-as.data.frame(table(ff$biochem_groups))
  g<- rename.vars( g, c("Freq") , c(" samples that have one hex a and hexb mutation " ))
  print(g)
  
  # mutations in either hex a or hex b
  multiple_mutations <- multiple_mutations [!multiple_mutations$Indiv %in% common_Indiv$intersect.multiple_mutations_HEXA.Indiv..multiple_mutations_HEXB.Indiv., ]
  
  multiple_mutations_HEXA<-subset(multiple_mutations, grepl("^15", Chrom))
  multiple_mutations_HEXA<-unique(multiple_mutations_HEXA, by=c(1:2))
  r<-as.data.frame(table(multiple_mutations_HEXA$biochem_groups))
  r<- rename.vars( r, c("Freq") , c(" hex a " ))
  print(r)
  
  multiple_mutations_HEXB<-subset(multiple_mutations, grepl("^5", Chrom))
  multiple_mutations_HEXB<-unique(multiple_mutations_HEXB, by=c(1:2))
  y<-as.data.frame(table(multiple_mutations_HEXB$biochem_groups))
  y<- rename.vars( y, c("Freq") , c(" hex b" ))
  print(y)
}
# race_list <- do.call(rbind, race_list)


#### variants that need ethnicla info and enzyme classification ----

cDNA=list("c.1274_1277dupTATC", "c.739C>T")

df_total = data.frame()

for (var in cDNA){
  name<-subset( reclassified_vars, cDNA==var)
  df <- data.frame(table(name$race))
  df_total <- rbind(df_total,df)
  
}
var72642925<-subset(reclassified_vars, cDNA=="c.*81_*82delTG" & gt_GT=="0/1")
table(var72642925$biochem_groups)
table(var72642925$race)

merged_biochem_molecular_data<-merge(biochem,molecular_without_path, by = "Indiv",all=FALSE, sort = FALSE)

var72642925<-subset(merged_biochem_molecular_data, ChromPOS=="5:74017080" & gt_GT=="1/1")
table(var72642925$race)
table(var72642925$biochem_groups)

#### Merging data ----
merged_biochem_molecular_data<-merge(biochem,molecular, by = "Indiv",all=FALSE, sort = FALSE)

#### Checking for missing samples in biochemical data ----

missing_in_biochem<-molecular%>%
  anti_join(biochem, by=c("Indiv"))
missing_samples_in_biochem<-as.data.table(unique(missing_in_biochem$Indiv))
export(list(missing_samples_in_biochem=missing_samples_in_biochem), "missing_samples_in_biochem.xlsx")


#### Checking for missing samples in molecular data ----

missing_in_molecular<-biochem%>%
  anti_join(molecular, by=c("Indiv"))

# Ordering missing samples

missing_samples_in_molecular<-as.data.table(unique(missing_in_molecular$Indiv))
missing_samples_in_molecular$V1<-as.numeric(missing_samples_in_molecular$V1)
missing_samples_in_molecular<-arrange(missing_samples_in_molecular,desc(V1))

export(list(missing_samples_in_molecular=missing_samples_in_molecular), "missing_samples_in_molecular.xlsx")

#### Excluding samples from biochem data set ----

biochem <- biochem[!(biochem$Indiv %in% missing_samples_in_molecular$V1),]
# get rid of rows that have empty indiv
biochem <- biochem[-which(biochem$Indiv == ""), ]

save(biochem,file="data.biochem")



#### Classified variants table ----

# All the classified variants for HEXA and HEXB have been exported from dashboard. There have been changes made to the original file 
# where the locus has been changed by one base pair location in order to accurately merge this file to the vcf file. 
c.p. <-read.table("/Users/olgalukatskaya/Desktop/work_projects/SR-2738/DataFiles/c.p.refaltchrompos.txt",header=TRUE, sep="\t")

# get rid of extra transcript 
c.p.<- subset(c.p., Transcript!= "NM_001292004.1")
# drop extra level "NM_001292004.1" in transcript column
c.p.$Transcript<-droplevels(c.p.$Transcript)

c.p.<- subset(c.p., select = c(ChromPOS,Classification,cDNA,REF,ALT,Protein) )
c.p.<- rename.vars( c.p., "Protein" , "AA")


# import new table with new variant classifications. These variants were classified by Jun. 
newdata<-read.xlsx("/Users/olgalukatskaya/Desktop/work_projects/SR-2738/DataFiles/variant_list_HEXA_HEXB_Unclassified_Variants_refaltchrompos.xlsx")
# divide a singel column into two columns where one columns contains  classifications and second columns has extra info that i dont need 
new_classifiedvar<-newdata %>% separate(Labels, into=c("Classification", "other"),sep=",")

#rename columns in new_classifiedvar
new_classifiedvar<- rename.vars( new_classifiedvar, c("Ref2","Alt2","cDNA.(cNomen)", "Protein.(pNomen)") , c("REF","ALT","cDNA","AA"))

# subset the data table and select only specific columns
new_classifiedvar<-subset(new_classifiedvar,select=c( ChromPOS,REF, ALT,cDNA,Classification,AA))

# find variants that are common to both tables and remove them from the c.p. data set
common_variant<-intersect(new_classifiedvar$ChromPOS,c.p.$ChromPOS)
common_variant<-common_variant[-(4)]
c.p.1 <- c.p.[ ! c.p.$ChromPOS %in% common_variant, ]

# bind two tables by adding one on top of the other by common columns 
variant_classification<-rbind(c.p.1, new_classifiedvar)
#creating a extra column genotype by combining REF and ALT columns togethers 
variant_classification_1<- variant_classification %>% unite(Genotype, REF,ALT, sep="/", remove=FALSE)


# import new table with new variant classifications. These variants were classified by Jun on 111617. 
unclassified_variants_111617<-read.xlsx("/Users/olgalukatskaya/Desktop/work_projects/SR-2738/DataFiles/extrs_variant_classification/unclassified_variants_111617.xlsx")

unclassified_variants_2<-unclassified_variants_111617 %>% separate(cDNA, into=c("other", "cDNA"),sep=":")

unclassified_variants_2<-subset(unclassified_variants_2,select=c( ChromPOS,REF, ALT,cDNA,`Jun's.Classification`, AA))

unclassified_variants_2a<-unite(unclassified_variants_2, "Genotype", REF, ALT, sep="/", remove = FALSE)
unclassified_variants_2a<- rename.vars( unclassified_variants_2a, c("Jun's.Classification") , c("Classification" ))

complete_list_of_classified_variants<-rbind( variant_classification_1, unclassified_variants_2a)

# get rid of specific var from the table     
complete_list_of_classified_variants<-complete_list_of_classified_variants[!(complete_list_of_classified_variants$ChromPOS=="5:74011541" & complete_list_of_classified_variants$cDNA=="c.1082+27delA"),]
complete_list_of_classified_variants<-complete_list_of_classified_variants[!(complete_list_of_classified_variants$ChromPOS=="5:74011541" & complete_list_of_classified_variants$cDNA=="c.1082+26_1082+27insA"),]
complete_list_of_classified_variants<-complete_list_of_classified_variants[!(complete_list_of_classified_variants$ChromPOS=="5:74011541" & complete_list_of_classified_variants$cDNA=="c.1082+26_1082+27insAAAA"),]

export(list(complete_list_of_classified_variants=complete_list_of_classified_variants), "complete_list_of_classified_variants.xlsx")


#### Excluding non-target regions from the variant classification table ----

#editing the table with the hexa/b targe regions
target_regions<-read.xlsx("/Users/olgalukatskaya/Desktop/work_projects/SR-2738/Jun_analysis_details /HEXA_HEXB_Target_Regions.xlsx",colNames=FALSE,1)
target_regions<-target_regions %>% separate(X1, c('h', 'Chrom'), sep="r")
target_regions<- subset(target_regions, select=-c(h))

#separating the molecular table and target region  table based on 2 chromosomes 
regions_5 <- subset(target_regions, grepl("^5", Chrom))
regions_15 <- subset(target_regions, grepl("^15", Chrom))
complete_list_of_classified_variants<-complete_list_of_classified_variants %>% separate(ChromPOS, c('Chrom', 'POS'), sep=":")

complete_list_of_classified_variants_15 <- subset(complete_list_of_classified_variants, grepl("^15", Chrom))
complete_list_of_classified_variants_5 <- subset(complete_list_of_classified_variants, grepl("^5", Chrom))


# getting all the values that are not in the target regions
#listing all the target sequences from start(X2) to end(X3) and stuffing them into a vector and  excluding regions from the 
# molecular table so that now we are left with non target regions that we want to remove
#15
complete_list_of_classified_variants_non_target_15<-as.data.frame(setdiff(complete_list_of_classified_variants_15$POS,unlist(with(regions_15, Map( `:`, X2, X3)))))
# rename vars 
complete_list_of_classified_variants_non_target_15<- rename.vars( complete_list_of_classified_variants_non_target_15, from="setdiff(complete_list_of_classified_variants_15$POS, unlist(with(regions_15, Map(`:`, X2, X3))))", to = "var1")
# removing non target reigons from the molecular table   
mo <- complete_list_of_classified_variants_15[ !complete_list_of_classified_variants_15$POS %in% complete_list_of_classified_variants_non_target_15$var1, ]

#5
complete_list_of_classified_variants_non_target_5<-as.data.frame(setdiff(complete_list_of_classified_variants_5$POS,unlist(with(regions_5, Map( `:`, X2, X3)))))
# rename vars 
complete_list_of_classified_variants_non_target_5<- rename.vars( complete_list_of_classified_variants_non_target_5, from="setdiff(complete_list_of_classified_variants_5$POS, unlist(with(regions_5, Map(`:`, X2, X3))))", to = "var1")
# removing non target reigons from the molecular table   
ma <- complete_list_of_classified_variants_5[ ! complete_list_of_classified_variants_5$POS %in% complete_list_of_classified_variants_non_target_5$var1, ]

bb <- rbind(mo, ma)
complete_list_of_classified_variants_within_target<- bb %>% unite(ChromPOS, Chrom,POS, sep=":", remove=TRUE)

# change the classification to make it more consistent 
complete_list_of_classified_variants_within_target <- within(complete_list_of_classified_variants_within_target, Classification[Classification == 'Likely benign'] <- 'Likely Benign')
complete_list_of_classified_variants_within_target$Classification<-droplevels(complete_list_of_classified_variants_within_target$Classification)
table(complete_list_of_classified_variants_within_target$Classification)

levels(complete_list_of_classified_variants_within_target$Classification) <- c(levels(complete_list_of_classified_variants_within_target$Classification),"Pseudodeficienct")

#### Make changes in the classified variant table ---- 
complete_list_of_classified_variants_within_target<- within(complete_list_of_classified_variants_within_target,
                                                            Classification[Classification == 'VOUS' & ChromPOS == '15:72636487' ] <- 'Likely Benign')

complete_list_of_classified_variants_within_target<- within(complete_list_of_classified_variants_within_target, 
                                                            Classification[Classification == 'VOUS' & ChromPOS == '15:72637878' ] <- 'Likely Benign')

complete_list_of_classified_variants_within_target<- within(complete_list_of_classified_variants_within_target,
                                                            Classification[Classification == 'Likely Benign' & ChromPOS == '15:72638659' & cDNA == 'c.1338T>C' ] <- 'VOUS')

complete_list_of_classified_variants_within_target<- within(complete_list_of_classified_variants_within_target,
                                                            Classification[Classification == 'Likely Benign' & ChromPOS == '5:73981270' ] <- 'Benign')

complete_list_of_classified_variants_within_target<- within(complete_list_of_classified_variants_within_target,
                                                            Classification[Classification == 'Pathogenic' & ChromPOS == '15:72642925' ] <- 'Pseudodeficienct')

complete_list_of_classified_variants_within_target<- within(complete_list_of_classified_variants_within_target,
                                                            Classification[Classification == 'VOUS' & ChromPOS == '5:74016956' ] <- 'Pseudodeficienct')

complete_list_of_classified_variants_within_target<- within(complete_list_of_classified_variants_within_target, 
                                                            Classification[Classification == 'Pathogenic' & ChromPOS == '15:72642919' ] <- 'Pseudodeficienct')

save(complete_list_of_classified_variants_within_target,file="data.complete_list_of_classified_variants_within_target")


