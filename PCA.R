# install.packages('...')
library("factoextra")
library("dplyr")
library("tidyr")
library("RColorBrewer")
library("pheatmap")
library("caret")

options(stringsAsFactors = FALSE)

#input tables
HPV_I=read.csv2('/Users/bichmann/Projects/LenaHPVCluster/OPSCC-Cluster/TemplateUltimate_Proteins_ClassI.csv')
HPV_II=read.csv2('/Users/bichmann/Projects/LenaHPVCluster/OPSCC-Cluster/TemplateUltimate_Proteins_ClassII.csv')

#specify which patients are HPV positive and HPV negative
nPats<-c('OPSCC_1003','OPSCC_1012','OPSCC_1014','OPSCC_1021','OPSCC_1024','OPSCC_1008','OPSCC_1382','OPSCC_1404','OPSCC_1393','OPSCC_1401','OPSCC_1403', 'OPSCC_04-001', 'OPSCC_1427')
pPats<-c('OPSCC_007','OPSCC_020','OPSCC_022','OPSCC_1118','OPSCC_1209','OPSCC_1349','OPSCC_1359','OPSCC_1345','OPSCC_1361','OPSCC_1372','OPSCC_1365','OPSCC_1385','OPSCC_1407','OPSCC_1416','OPSCC_1420')


#HPV_I[c('Patient','X...Protein.Group.Accessions')]
#HPV_II[c('Patient','X...Accession')]

#proteins=c('P23470','IL17REL','Q14993','Q5IJ48','Q96T17','GRIN2C',
#           'Q6PF15','Q9UPR6','Q8WWU5','P19012','O75346','Q9BXR5',
#          'Q9H0D2','RNF212','Q9UJ83','Q9UJ98','Q5TIA1','Q0VDD7',
#           'Q9UKT9','Q9BX26','A0A075B6H9','O43761','P02788','P42771',
#           'P42773','P33076','P31785','Q99758','P43360','P43364',
#           'O60663','Q9Y5F6','P01138','P08151','P78396','Q07283',
#           'P23141','P29353','Q96J84','P28290','Q5XXA6','O95678',
#           'P21589','Q86SJ2','P04637','O15391','P12004','Q9UJY4',
#           'Q01892','P41217','Q9UKQ2','P10415','P19320','O75144',
#           'P31749','P31751','Q9Y243','P42345','Q92569','O00459',
#           'P27986','P02795','P14780','P17948','P15692','Q01860',
#           'O95715','P29590','P35222','Q9UP38','Q14332','P04426',
#           'P09544','Q15465','Q14623','O43323','P10070','P10071')

#read and store tables in data frame HPVmatrix
HPVmatrix=as.data.frame(matrix(0,length(unique(HPV_I$Patient)),length(unique(c(HPV_I$Protein.Group.Accessions, HPV_II$Accession)))))
row.names(HPVmatrix)=unique(HPV_I$Patient)
colnames(HPVmatrix)=unique(c(HPV_I$Protein.Group.Accessions, HPV_II$Accession))

r1 <- HPV_I %>% group_by(Patient) %>% summarise(protein = list(unique(Protein.Group.Accessions)))
r1 %>% ungroup %>% unnest()

r2 <- HPV_II %>% group_by(Patient) %>% summarise(protein = list(unique(Accession)))
r2 %>% ungroup %>% unnest()

for (pat in seq(1,length(unique(HPV_I$Patient)))) {
  print(pat)
  plist1<-r1[pat,]$protein
  HPVmatrix[r1[pat,]$Patient,unlist(plist1)]<-1
  plist2<-r2[pat,]$protein
  HPVmatrix[r2[pat,]$Patient,unlist(plist2)]<-1
}

#Throw out proteins that have variance 0 across all patients
HPVmatrix<-HPVmatrix[ - as.numeric(which(apply(HPVmatrix, 2, var) <= 0))]
dim(HPVmatrix)

#Add column describing infection state
HPVmatrix$Infection<-rep('unkown',dim(HPVmatrix)[1])
HPVmatrix$Infection[which(rownames(HPVmatrix) %in% nPats)]<-'other'
HPVmatrix$Infection[which(rownames(HPVmatrix) %in% pPats)]<-'HPV_infected'


# Machine learning script to define import proteins for distinction of HPV negative and positive
HPVmatrixML<-HPVmatrix[which(HPVmatrix$Infection!='unkown'),]
HPVstat<-HPVmatrixML$Infection
HPVmatrixML<-HPVmatrixML[,-dim(HPVmatrix)[2]]

param <- c(100) 
tunegrid <- expand.grid(.mtry=param) 

fitControl <- trainControl(method = "cv",
                           ## 10-fold CV...
                           number = 10,
                           classProbs = TRUE,
                           ## repeated ten times
                           #repeats = 10,
                           summaryFunction = twoClassSummary)

Training <- train(x=HPVmatrixML, y=as.factor(HPVstat), 
                  method = "rf", 
                  importance = TRUE,
                  #preProc = c("center"),
                  tuneGrid = tunegrid,
                  trControl = fitControl, metric = "ROC")

importance<-as.data.frame(varImp(Training)[1])

#Importance threshold set to 60
#Check Importance curve using plot(varImp(Training))
important_proteins<-rownames(importance)[which(importance$importance.HPV_infected>=60)]

#If one wants to check all proteins set important proteins to all
#important_proteins<-colnames(HPVmatrix)

#Visualize results using PCA
res.pca <- prcomp(HPVmatrix[,which(colnames(HPVmatrix) %in% important_proteins)], scale = FALSE) #[,which(colnames(HPVmatrix) %in% proteins)],-dim(HPVmatrix)[2]]

#Visualize variance explained by PCA axis
fviz_eig(res.pca)

#Visualize PCA
fviz_pca_ind(res.pca,
             #col.ind = "cos2", # Color by the quality of representation
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             #repel = TRUE     # Avoid text overlapping
             habillage = HPVmatrix$Infection, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses,
             ellipse.level = 0.7
)

#Visualize all contributions of all proteins in PCA
fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)

#Visualize contributions of individual proteins
fviz_contrib(res.pca, choice = "var", axes = 1, top = 50)
fviz_contrib(res.pca, choice = "var", axes = 2, top = 50)
col.pal <- brewer.pal(9,"Blues")

#Visualize all important proteins as heatmap
annotation_row <- data.frame(Infection = HPVmatrix$Infection)
rownames(annotation_row)=rownames(HPVmatrix)
pheatmap(HPVmatrix[,which(colnames(HPVmatrix) %in% important_proteins)], show_colnames = TRUE, color = col.pal, annotation_row = annotation_row )

