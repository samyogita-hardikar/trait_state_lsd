.libPaths(c("/data/tu_hardikar/R/x86_64-pc-linux-gnu-library/4.0", "/afs/cbs.mpg.de/software/r/4.0.2/ubuntu-bionic-amd64/lib/R/library"))

#setwd("W:/MPI/Daydreaming/LSD3/")
setwd("/data/p_02260/all_lsd/gradients/")

# load necessary packages -------------------------------------------------

library(FactoMineR)
library(ggplot2)
library(factoextra)
library(reshape2)
library(missMDA)
library(VIM)
library(dplyr)
library(car)
library(ggiraph)
library(ggiraphExtra)
#library(MVN)
library(gplots)
library(ComplexHeatmap)
library(dendextend)
library(corrplot)
#library(MVLM)
library(heplots)
library(viridis)
library(ggwordcloud)
library(circlize)
library(RColorBrewer)
library(lm.beta)
library(psych)
library(tibble)
library(PerformanceAnalytics)



# data import and cleaning ------------------------------------------------

folder <- "dataverse_files/behavioral_data_MPILMBB/phenotype/"   # path to folder that holds multiple .tsv files 

spss_pca_scores<-read.csv("trait_LSD_spss_pca_scores.csv", sep= ",", header = TRUE)
spss_pca_rotated_var_scores<-read.csv("trait_LSD_spss_var_scores_rotated.csv", sep= ";", header = TRUE)
spss_pca_unrotated_var_scores<-read.csv("trait_LSD_spss_var_scores_unrotated.csv", sep= ";", header = TRUE)

thought_spss_pca_scores<-read.csv("thought_LSD_spss_pca_scores.csv", sep= ",", header = TRUE)
thought_spss_pca_rotated_var_scores<-read.csv("thought_LSD_spss_var_scores_rotated.csv", sep= ";", header = TRUE)
thought_spss_pca_unrotated_var_scores<-read.csv("thought_LSD_spss_var_scores_unrotated.csv", sep= ";", header = TRUE)



# create list of all .tsv files in folder
file_list <- list.files(path=folder, pattern="*.tsv") 
# create list df names (without.tsv)
df_list<-gsub(".tsv", "", file_list)

# read in each .tsv file in file_list and create a data frame with the same name
for (i in 1:length(file_list)){
  df<-assign(df_list[i],
             read.csv(file = file.path(folder, file_list[i]), sep='\t', na.strings= "n/a"))
  df[[1]]<-gsub("sub-0", "", df[[1]])         #remove "sub-0" 
  colnames(df)[1] <- "participant_id"         # change first colum to "participant_id"
  assign(df_list[i], df)
  remove(df)
}

#minor adjustments 
myBDI<-BDI[c("participant_id", "BDI_summary_sum")]      #get the BDI summary score
myNEO<-NEO[c("participant_id", "NEO_N", "NEO_E", "NEO_O", "NEO_A", "NEO_C")]         #Keep only the NEO-5 

colnames(ASR)<-c("participant_id", "Adaptive Func Friends (ASR)", 
                 "Adaptive Func Spouse (ASR)", "Adaptive Func Family (ASR)",           
                 "Adaptive Func Job (ASR)", "Adaptive Func Edu (ASR)",              
                 "Tobacco (ASR)", "Alcohol(ASR)", "Drugs(ASR)", 
                 "Critical Items (ASR)","Ancious-depressed (ASR)",
                 "Withdrawn (ASR)", "Somatic Complaints (ASR)",
                 "Thought Problems (ASR)", "Attention Problems (ASR)",
                 "Aggressive Beh (ASR)", "Rule-breaking Beh (ASR)",
                 "Intrusive (ASR)", "Internalizing (ASR)", 
                 "Externalizing (ASR)")




#load raw motion parametres
mean_FD<-read.csv("mean_FD.csv", sep = ",",  header = TRUE)
mean_FD$participant_id<-gsub("sub-0", "", mean_FD$participant_id )    #remove "sub-0" 

#load gender
motion_age_gender<- read.csv("motion_age_gender.csv", sep = ",", header = TRUE)
gender<- motion_age_gender[, c("participant_id", "gender")]
gender$gender[which(gender$gender=="F")]<-1
gender$gender[which(gender$gender=="M")]<-0

#load age csv
age<-read.csv("age.csv", header = TRUE, sep = ",")

my_mag<-Reduce(function(x, y) merge(x, y, all=TRUE), 
               list(age, mean_FD,gender, by= "participant_id" ))
my_mag<-subset(my_mag[, 1:7])


# SNYCQ rearrage SNYCQ to cbind all instances properly
mySNYCQ<-as.data.frame(na.pass(SNYCQ[1:73]))
keepvars_mySNYCQ<-c("participant_id", colnames(mySNYCQ[grep("post.ses.02.run", colnames(mySNYCQ))]))
mySNYCQ<-SNYCQ[keepvars_mySNYCQ]

# Keep participants with only one set of SNYCQ missing
#mySNYCQ<-as.data.frame(mySNYCQ[which(rowSums(is.na(mySNYCQ))<13),])

# SNYCQ colname corrections
colnames(mySNYCQ)<-gsub("surrpundings", "surroundings", colnames(mySNYCQ))
colnames(mySNYCQ)<-gsub("vague", "specific", colnames(mySNYCQ))



#Only the psychometric questionnaires that I'm interested in
myLSD <- Reduce(function(x, y) merge(x, y, all=TRUE), 
                list(myBDI, BISBAS, BPS, ESS, GoldMSI, HADS, IAT, IMIS, 
                     MMI,  myNEO, PSSI,SCS, SD3, SDMW, SDS, 
                     SE, STAXI, TPS, UPPS, ASR, ACS, my_mag, by= "participant_id" ))

drop_y<-colnames(myLSD)%in% c("y")
myLSD<-myLSD[!drop_y] 

all_LSD<- Reduce(function(x, y) merge(x, y, all=TRUE), 
                 list(myLSD, mySNYCQ,  by= "participant_id" ))

drop_y<-colnames(all_LSD)%in% c("y")
all_LSD<-all_LSD[!drop_y] 

#all the cleaning steps: 
#remove older, remove missing, replace bad scan thoughts with NAs

#participant with 3 scans with fewer timepoints, for some reason
remove_participants<-all_LSD$participant_id %in%c( "10070", "10156", "10184", "10198", "10214")
all_LSD_clean<-all_LSD[-which(remove_participants),]

all_LSD_clean<-all_LSD_clean[which(rowSums(is.na(
  all_LSD_clean[,grep("SNYCQ", colnames(all_LSD_clean))]))<14),]

remove_observations_rows<- c("10027","10029","10030","10038","10044","10045","10083","10090",
                             "10099", "10101","10103","10133","10164","10201","10211")

remove_observations_columns<-c(grep("SNYCQ_post.ses.02.run.02.acq.PA", colnames(all_LSD_clean)),
                               grep("acq.PA_run.02", colnames(all_LSD_clean)))
all_LSD_clean[which(all_LSD_clean$participant_id%in%remove_observations_rows), 
          remove_observations_columns]<-NA

all_LSD_clean_young<-all_LSD_clean[which(all_LSD_clean$age_years<50),]
all_LSD_clean_young<-all_LSD_clean_young[which(rowSums(is.na(all_LSD_clean_young))<126),]

motion<-as.numeric(as.vector(rowMeans(all_LSD_clean_young[c(
  grep("_run", colnames(all_LSD_clean_young)))], na.rm = TRUE)))   

all_LSD_clean_young<-subset(all_LSD_clean_young[-c(grep("_run", colnames(all_LSD_clean)))])

all_LSD_clean_young<-add_column(all_LSD_clean_young, motion, .after = "age_years") 

#write.table(all_LSD_clean_young, file = "all_lsd_clean_young.csv", sep= ",", na="NA", col.names= TRUE)

all_LSD_clean_young_participants<-as.list(paste("sub-0", all_LSD_clean_young$participant_id, sep=""))


### dataframe for trait PCA 
trait_LSD<-as.data.frame(all_LSD_clean_young[colnames(myLSD[1:73])])
#Renames Variables as you want them for eventual plots
colnames(trait_LSD)[2]<-"Depression (BDI)"
colnames(trait_LSD)[3]<-"Inhibition Sum (BIS/BAS)"
colnames(trait_LSD)[4]<-"Activationn Sum (BIS/BAS)"
colnames(trait_LSD)[5]<-"Boredom Proneness (BPS)"
colnames(trait_LSD)[6]<-"Sleepiness (Epworth)"
colnames(trait_LSD)[7]<-"Active Engagement in Music (Gold MSI)"
colnames(trait_LSD)[8]<-"Musical Training (Gold MSI)"
colnames(trait_LSD)[9]<-"Anxiety Sum (HADS)"
colnames(trait_LSD)[10]<-"Depression Sum (HADS)"
colnames(trait_LSD)[11]<-"Internet Addiction (IAT)"
colnames(trait_LSD)[12]<-"Negative Valence (IMIS)"
colnames(trait_LSD)[13]<-"Help (IMIS)"
colnames(trait_LSD)[14]<-"Movement (IMIS)"
colnames(trait_LSD)[15]<-"Personal Reflections (IMIS)"
colnames(trait_LSD)[16]<-"Multimedia Multitasking"
colnames(trait_LSD)[17:21]<- c("Neuroticism (NEO PI-R)" ,"Extraversion (NEO PI-R)", "Openness (NEO PI-R)", 
                           "Agreeableness (NEO PI-R)", "Conscientiousness (NEO PI-R)")
colnames(trait_LSD)[22:35]<-c("Paranoid (PSSI)", "Schizoid (PSSI)", "Schizotypal (PSSI)", 
                          "Borderline (PSSI)", "Histrionic  (PSSI)", "Narcissistic (PSSI)", 
                          "Avoidant (PSSI)", "Dependent  (PSSI)", "Obsessive-compulsive (PSSI)",
                          "Negativistic (PSSI)", "Depressive (PSSI)", "Altruistic (PSSI)",
                          "Rhapsodic (PSSI)",  "Antisocial (PSSI)")
colnames(trait_LSD)[36]<-"Self Control Summary (SCS)"
colnames(trait_LSD)[37:39]<-c("Machievellianism Sum (SD3)", "Narcissism Sum (SD3)", 
                          "Psychopathy Sum (SD3)")
colnames(trait_LSD)[40:41]<- c("Deliberate (SDMW)", "Spontaneous (SDMW")
colnames(trait_LSD)[42]<-"Social Desirability Sum (SDS)"
colnames(trait_LSD)[43]<-"Global Self Worth (SES)"
colnames(trait_LSD)[44:47]<-c("Anger Trait (STAXI)", "Anger Inward (STAXI)", "Anger Outward (STAXI)",
                          "Anger Control (STAXI)")
colnames(trait_LSD)[48]<-"Procrastination Sum (TPS)"
colnames(trait_LSD)[49:53]<-c("Negative Urgency (UPPS)", "Lack of Premeditation (UPPS)", 
                          "Lack of Perseverance (UPPS)", "Sensation Seeking (UPPS)",
                          "Postitive Urgency (UPPS)")


#write.table(trait_LSD, file = "trait_LSD.csv", sep= ",", na="NA", col.names= TRUE)


#KMO(na.exclude(trait_LSD[,2:73]))

# trait_LSDnPCs <- estim_ncpPCA(trait_LSD[,2:73])
# PCA_trait_LSD_impute<- imputePCA(trait_LSD [,2:73 ], ncp=5, graph= FALSE)
# KMO(PCA_trait_LSD_impute$completeObs)  # Does PCA make sense for this data? (sampling adequacy)
# sampling_adequacy<-KMO(PCA_trait_LSD_impute$completeObs)  
# cortest.bartlett(PCA_trait_LSD_impute$completeObs, n=144)

PCA_trait_LSD<-subset(trait_LSD, select= -c(`Anger Control (STAXI)`))
PCA_trait_LSD_impute<- imputePCA(PCA_trait_LSD [,2:72 ], ncp=5, graph= FALSE) 
#write.table(PCA_trait_LSD_impute$completeObs, file = "trait_LSD_imputed.csv", sep= ",", na="", col.names= TRUE, row.names = FALSE)
res.pca = PCA(PCA_trait_LSD_impute$completeObs, ncp=5, graph= FALSE, scale.unit = TRUE )

#fviz_screeplot(res.pca, ncp=10, addlabels = TRUE, main = "Traits PCA") 
fviz_screeplot(res.pca, ncp=10, addlabels = TRUE, main = "Traits PCA", labelsize = 8, repel = TRUE) +
  theme(text = element_text(size = 22),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 22))
#export at 700 x 700




#cortest.bartlett(PCA_trait_LSD_impute$completeObs)




#Heatmap of component loadings -  export at  788 x 980
# 
# col_fun = colorRamp2(c(-1, 0, 1), c(viridis(3)))
# 
# row_dend = hclust(dist(res.pca$var$coord)) # row clustering
# col_dend = hclust(dist(t(res.pca$var$coord))) # column clustering
# Heatmap(res.pca$var$cor, col=col_fun,show_row_dend = FALSE, show_column_dend = FALSE,
#         name = "Loadings", row_names_side = "left",
#         row_names_gp = gpar(fontsize=7),
#         width = unit(4, "cm"), height = unit(15, "cm"),
#         cluster_rows = color_branches(row_dend), cluster_columns = color_branches(col_dend),
#         column_labels = c("C1", "C2", "C3", "C4", "C5"),
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.1f", res.pca$var$cor[i, j]), x, y, gp = gpar(fontsize = 6))
#         })
#
# 
# 
# 

## VARIMAX rotated

rotated_loadings<-varimax(res.pca$var$coord)

col_fun = colorRamp2(c(-1, 0, 1), c(viridis(3)))

row_dend = hclust(dist(rotated_loadings$loadings)) # row clustering
col_dend = hclust(dist(t(rotated_loadings$loadings))) # column clustering
Heatmap(as.matrix(rotated_loadings$loadings[1:71,1:5]),
        col=col_fun,show_row_dend = FALSE, show_column_dend = FALSE,
        name = "Loadings", row_names_side = "left", 
        row_names_gp = gpar(fontsize=11),
        column_names_gp = gpar(fontsize=11),
        width = unit(8, "cm"), height = unit(30, "cm"),
        cluster_rows = color_branches(row_dend), 
        cluster_columns = FALSE,
        column_labels = c("Dim_1", "Dim_2", "Dim_3",
                          "Dim_4", "Dim_5"),
        column_names_rot = 45,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", rotated_loadings$loadings[i, j]), x, y, gp = gpar(fontsize = 11))
        }
)




# #SPSS rotated loadings 

spss_rotated_loadings<-as.matrix(spss_pca_rotated_var_scores[1:71,2:6]) 

col_fun = colorRamp2(c(-1, 0, 1), c(viridis(3)))

dimnames(spss_rotated_loadings)[[1]]<- spss_pca_rotated_var_scores[1:71,1]
dimnames(spss_rotated_loadings)[[2]]<- colnames(spss_pca_rotated_var_scores)[2:6]

row_dend = hclust(dist(spss_rotated_loadings)) # row clustering
col_dend = hclust(dist(t(spss_rotated_loadings))) # column clustering
Heatmap(as.matrix(spss_rotated_loadings),
        col=col_fun,show_row_dend = FALSE, show_column_dend = FALSE,
        name = "Loadings", row_names_side = "left",
        row_names_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize=12),
        width = unit(5, "cm"), height = unit(18, "cm"),
        cluster_rows = color_branches(row_dend),
        cluster_columns = FALSE,
        column_labels = c("Trait_1", "Trait_2", "Trait_3",
                          "Trait_4", "Trait_5"),
        column_names_rot = 45,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", spss_rotated_loadings[i, j]), x, y, gp = gpar(fontsize = 8))
        }
)
#exported at w 600 h 800






### SNYCQ PCA-----------------------------------------

mynewSNYCQ<-as.data.frame(all_LSD_clean_young[,c(1,grep("SNYCQ_post.ses.02", colnames(all_LSD_clean_young)))])


SNYCQ_varlist<-list("positive", "negative", "future", "past", "myself", 
                    "people", "surroundings", "vigilance", "images", "words", "specific", "intrusive")



for (i in 1:length(SNYCQ_varlist)){
  newdf<-as.data.frame(cbind(mynewSNYCQ$participant_id, 
                             mynewSNYCQ[grep(SNYCQ_varlist[i],colnames(mynewSNYCQ), ignore.case = TRUE)]), na.rm=TRUE)
  colnames(newdf)<- c("participant_id", paste0("newSNYCQ1_", SNYCQ_varlist[[i]]),paste0("newSNYCQ2_", SNYCQ_varlist[[i]]),
                      paste0("newSNYCQ3_", SNYCQ_varlist[[i]]),paste0("newSNYCQ4_", SNYCQ_varlist[[i]]))
  newdf<-melt(newdf, id.vars="participant_id")
  colnames(newdf)<- c("participant_id", "var", SNYCQ_varlist[[i]])
  newdf<- newdf[c(1,3)]
  assign(paste0("newSNYCQ_", SNYCQ_varlist[[i]]), newdf)
}



newSNYCQ<-as.data.frame(cbind(newSNYCQ_future, newSNYCQ_images[,2], newSNYCQ_intrusive[,2], newSNYCQ_myself[,2], newSNYCQ_negative[,2], 
                              newSNYCQ_past[,2], newSNYCQ_people[,2], newSNYCQ_positive[,2], newSNYCQ_surroundings[,2], 
                              newSNYCQ_specific[,2], newSNYCQ_vigilance[,2], newSNYCQ_words[,2]))


newcolnames<-gsub("newSNYCQ_", "", colnames(newSNYCQ))
newcolnames<-gsub("[, 2]", "", newcolnames)
newcolnames<-gsub("\\[\\]", "", newcolnames)
colnames(newSNYCQ)<-newcolnames

#write.table(newSNYCQ, file = "thoughts_LSD.csv", sep= ",", na="NA", col.names= TRUE)

#KMO(na.exclude(trait_LSD[,2:73]))
thoughts_LSD<-as.data.frame(newSNYCQ)
thoughts_LSDnPCs <- estim_ncpPCA(thoughts_LSD[,2:13])
thoughts_LSD_impute<- imputePCA(thoughts_LSD[,2:13 ], ncp=5, graph= FALSE) 
#write.table(thoughts_LSD_impute$completeObs, file = "thought_LSD_R_imputed.csv", sep= ",", na="", col.names= TRUE, row.names = FALSE)

#PCA_thoughts_LSD_impute<- imputePCA(thoughts_LSD [,2:13 ], ncp=5, graph= FALSE) 
KMO(thoughts_LSD_impute$completeObs)  # Does PCA make sense for this data? (sampling adequacy)
cortest.bartlett(thoughts_LSD_impute$completeObs)

res.pca.thoughts = PCA(thoughts_LSD_impute$completeObs, ncp=5,graph= FALSE) 
#fviz_screeplot(res.pca.thoughts, ncp=10, addlabels = TRUE, main = "Thoughts Patterns PCA") 

fviz_screeplot(res.pca.thoughts, ncp=10, addlabels = TRUE, main = "Thoughts PCA", labelsize = 8, repel = TRUE) +
  theme(text = element_text(size = 22),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 22))
#export at 700 x 700


row_dend = hclust(dist(res.pca.thoughts$var$coord)) # row clustering
col_dend = hclust(dist(t(res.pca.thoughts$var$coord))) # column clustering
Heatmap(as.matrix(res.pca.thoughts$var$coord),
        col=col_fun, show_row_dend = FALSE, show_column_dend = FALSE,
        name = "Loadings", row_names_side = "left", 
        row_names_gp = gpar(fontsize=9),
        width = unit(4, "cm"), height = unit(15, "cm"),
        cluster_rows = color_branches(row_dend), 
        cluster_columns = FALSE,
        column_labels = c("Thought_1", "Thought_2", "Thought_3",
                          "Thought_4", "Thought_5"),
        column_names_rot = 45,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", res.pca.thoughts$var$coord[i, j]), x, y, gp = gpar(fontsize = 7))
        }
)



##rotated loadings

rotated_loadings_thoughts<-varimax(res.pca.thoughts$var$coord)

col_fun = colorRamp2(c(-1, 0, 1), c(viridis(3)))

row_dend = hclust(dist(rotated_loadings_thoughts$loadings)) # row clustering
col_dend = hclust(dist(t(rotated_loadings_thoughts$loadings))) # column clustering
Heatmap(as.matrix(rotated_loadings_thoughts$loadings[1:12,1:5]),
        col=col_fun, show_row_dend = FALSE, show_column_dend = FALSE,
        name = "Loadings", row_names_side = "left", 
        row_names_gp = gpar(fontsize=9),
        width = unit(5, "cm"), height = unit(18, "cm"),
        cluster_rows = color_branches(row_dend), 
        column_names_gp = gpar(fontsize=12),
        cluster_columns = FALSE,
        column_labels = c("Thought_1", "Thought_2", "Thought_3",
                          "Thought_4", "Thought_5"),
        column_names_rot = 45,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", rotated_loadings_thoughts$loadings[i, j]), x, y, gp = gpar(fontsize = 8))
        }
)


## 

spss_rotated_loadings_thoughts<-as.matrix(thought_spss_pca_rotated_var_scores[1:12,2:6])

dimnames(spss_rotated_loadings_thoughts)[[1]]<- thought_spss_pca_rotated_var_scores[1:12,1]
dimnames(spss_rotated_loadings_thoughts)[[2]]<- colnames(thought_spss_pca_rotated_var_scores)[2:6]

col_fun = colorRamp2(c(-1, 0, 1), c(viridis(3)))

row_dend = hclust(dist(spss_rotated_loadings_thoughts)) # row clustering
col_dend = hclust(dist(t(spss_rotated_loadings_thoughts))) # column clustering
Heatmap(as.matrix(spss_rotated_loadings_thoughts),
        col=col_fun,show_row_dend = FALSE, show_column_dend = FALSE,
        name = "Loadings", row_names_side = "left",
        row_names_gp = gpar(fontsize=10),
        column_names_gp = gpar(fontsize=12),
        width = unit(5, "cm"), height = unit(18, "cm"),
        cluster_rows = color_branches(row_dend),
        cluster_columns = FALSE,
        column_labels = c("Thought_1", "Thought_2", "Thought_3",
                          "Thought_4", "Thought_5"),
        column_names_rot = 45,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", spss_rotated_loadings_thoughts[i, j]), x, y, gp = gpar(fontsize = 8))
        }
)








for (i in 1:length(SNYCQ_varlist)){
  df<-as.data.frame(cbind(mynewSNYCQ$participant_id, (rowMeans(mynewSNYCQ[grep(SNYCQ_varlist[i],
                                                                         colnames(mynewSNYCQ),  ignore.case = TRUE)], na.rm=TRUE))))
  colnames(df)<- c("participant_id", paste0("SNYCQ_", SNYCQ_varlist[[i]]))
  assign(paste0("SNYCQ_", SNYCQ_varlist[[i]]), df)
}



mynewSNYCQ_avg<- Reduce(function(x, y) merge(x, y, all=TRUE, na.strings= "NaN"), 
                     list(SNYCQ_future, SNYCQ_images, SNYCQ_intrusive, SNYCQ_myself, SNYCQ_negative, SNYCQ_past,
                          SNYCQ_people, SNYCQ_positive, SNYCQ_surroundings, SNYCQ_specific, SNYCQ_vigilance, SNYCQ_words, by= "participant_id"))

mynewSNYCQ_avg<-type.convert(mynewSNYCQ_avg[1:13], na.strings= "NaN")
colnames(mynewSNYCQ_avg)<-c(gsub("SNYCQ_", "", colnames(mynewSNYCQ_avg)))

#mything5<-as.numeric(res.pca.thoughts$ind$coord[,5])

thoughts_LSD_coord<-as.data.frame(cbind(newSNYCQ$participant_id, thought_spss_pca_scores$FAC1_2,
                                        thought_spss_pca_scores$FAC2_2, thought_spss_pca_scores$FAC3_2,
                                        thought_spss_pca_scores$FAC4_2,thought_spss_pca_scores$FAC5_2)) 

thoughts_LSD_coord[] <- lapply(thoughts_LSD_coord, function(x) as.numeric(as.character(x))) 

thoughts_LSD_coord_mean<-aggregate(thoughts_LSD_coord[, 2:6], list(thoughts_LSD_coord$V1 ), mean)
colnames(thoughts_LSD_coord_mean)<-c("participant_id", "Thought_C1", "Thought_C2", "Thought_C3", "Thought_C4", "Thought_C5")



LSD<-as.data.frame(cbind(all_LSD_clean_young[1], spss_pca_scores[77:81],
                             thoughts_LSD_coord_mean[2:6], mynewSNYCQ_avg[2:13],
                         all_LSD_clean_young[c(grep("age_years", colnames(all_LSD_clean_young)))],
                         all_LSD_clean_young[c(grep("motion", colnames(all_LSD_clean_young)))],
                         all_LSD_clean_young[c(grep("gender", colnames(all_LSD_clean_young)))]))

LSD[2:26] <- lapply(LSD[2:26], function(x) as.numeric(as.character(x))) 

LSD_std<-as.data.frame(cbind(LSD[1], scale(LSD[2:23]), LSD[24:26]))



#write.table(LSD_std, file = "LSD_std_spss.csv", sep= ",", na="NA", col.names= TRUE)



## save imputed, unstandardized PCA rotated scored from spss as CSV
newmat<-cbind(LSD$participant_id,  spss_pca_scores[77:81],   thoughts_LSD_coord_mean[2:6] )

my_outlier_impute <- function(x) {
  x[x < quantile(x,0.25) - 1.5 * IQR(x) | x > quantile(x,0.75) + 1.5 * IQR(x)] <- 
    mean(x[!(x < quantile(x,0.25) - 1.5 * IQR(x) | x > quantile(x,0.75) + 1.5 * IQR(x))])
  x
}

newmat_new<-cbind(newmat[1], lapply(newmat[2:11] , my_outlier_impute))
write.table(newmat_new, file = "LSD_spss_rotated_imputed_unstd.csv", sep= ",", na="NA", col.names= TRUE, row.names = FALSE)


outlier <- function(x) {
  x[x < quantile(x,0.25) - 1.5 * IQR(x) | x > quantile(x,0.75) + 1.5 * IQR(x)] <- mean(x)
  x
}

imputed_vars<-lapply(LSD_std[2:23], outlier)
LSD_std_imputed<-cbind(LSD_std[1], lapply(LSD_std[2:23], outlier), LSD_std[24:26])

write.table(LSD_std_imputed, file = "LSD_std_imputed.csv", sep= ",", na="NA", col.names= TRUE)
