##### Librarying packages #####
library(readxl)
library(tidyverse)
library(paco)
library(devtools)
library(treeman)
library(ape)
library(phangorn)
library(seqinr)
library(ips)
library(ape)
library(ips)
library(MCMCtreeR)
library(treeio)
library("seqinr")
library(phytools)
library(data.table)
library(vegan)
library(ggplot2)
library(brms)
library(tidytree)
setwd("~/Lab Poulin/PhD/PhD/Dados/Capítulo 4/Tabelas")

##### Assembling the data #####

dados <- read_excel("Table and Sites 2021 Malavi.xlsx")
dados <- dados[-c(9:10, 14:18, 20:22, 24:25)]
dados <- filter(dados, Host_Environment == "Wild")
Birdtraits <- read_excel("MigTraits.xlsx")

MyTraits <- Birdtraits %>%
  filter(Sp.Scien.jetz %in% dados$species)

MyTraits <- MyTraits[-c(2:6)]
names(MyTraits)[names(MyTraits) == "Sp.Scien.jetz"] <- "species"

dados1 <- inner_join(dados, MyTraits, by = "species")

write.csv2(file = "Cap4.csv", dados1)

##### Filtering the Data #####

dados2 <- with(dados1, names(table(species)[table(species) > 4]))
dados2 <- dados1[dados1$species%in% dados2, ]
dados2$species <- factor(dados2$species)

dados3 <- with(dados2, names(table(Lineage)[table(Lineage) > 4]))
dados3 <- dados2[dados2$Lineage%in% dados3, ]
dados3$Lineage <- factor(dados3$Lineage)

dados3 <- filter(dados3, parasiteGenus != "Leucocytozoon")

dados3$strategy_3 <- as.factor(dados3$strategy_3)
length(unique(dados3$Lineage))# [dados3$parasiteGenus == "Haemoproteus"]))
length(unique(dados3$species))

dados3$Loc <-  paste(dados3$Latitude,dados3$Longitude, sep = "_")
length(unique(dados3$Loc))


##### Creating Fasta File for Parasite Phylogeny #####

teste <- read.FASTA("Malavi_Long.fasta", type = "DNA")
class(teste)
print(teste)

teste1 <- teste[names(teste) %in% dados3$Lineage]

write.FASTA(teste1, "HaemSeqCap4.fasta")

write.dna(teste1, "HaemSeqCap4.NEX", format = "sequential", append = FALSE,
          nbcol = 3, colsep = "", colw = 500, indent = NULL,
          blocksep = 1)

##### Adding Parasite Phylogeny #####

phylo <- read.mrbayes("HaemSeqCap4.tre")
class(phylo)
phylo <- as.phylo(phylo)
class(phylo)
phylo$edge.length
is.rooted(phylo)

phy3 <- root(phylo, "L_GALLUS05_Leucocytozoon_caulleryi", resolve.root=TRUE)
length(unique(phy3$tip.label))
is.rooted(phy3)
class(phy3)

matches2 <-match(dados3$Lineage, phy3$tip.label)
matches2 <-na.omit(matches2)
phylotree <-drop.tip(phy3, phy3$tip.label[-matches2])

sptest=as.factor(phylotree$tip.label)
sp=as.data.frame(sptest)

dados3=dados3 %>%
  filter(Lineage %in% sp$sptest)

##### Removing duplicated sequences #####

Seq <- read_excel("MalAvi_Sequences.xlsx")
Seq <- Seq[c(1,24)]
Seq <- Seq %>%
  filter(Lineage_Name %in% dados3$Lineage_Name)
length(unique(Seq$sequence))
Seq$sequence[duplicated(Seq$sequence)]

dados3$Lineage_Name[dados3$Lineage_Name == "SETAUD08"] <- "PASILI01"
dados3$Lineage_Name[dados3$Lineage_Name == "TULEU06"] <- "TULEU01"

dados3$Lineage[dados3$Lineage == "H_SETAUD08"] <- "H_PASILI01"
dados3$Lineage[dados3$Lineage == "P_TULEU06"] <- "P_TULEU01"

length(unique(phylotree$tip.label))
length(unique(dados3$Lineage))

matches2 <-match(dados3$Lineage, phy3$tip.label)
matches2 <-na.omit(matches2)
phylotree <-drop.tip(phy3, phy3$tip.label[-matches2])

length(unique(phylotree$tip.label))
length(unique(dados3$Lineage))

plot.phylo(phylotree, type = "fan", show.tip.label = FALSE)

##### Bird Phylogeny #####

load("random_trees2.RData")
mytree <- random_trees2[[1]]

matchesBirds <-match(dados3$species, mytree$tip.label)
matchesBirds <-na.omit(matchesBirds)
mytree1 <-drop.tip(mytree, mytree$tip.label[-matchesBirds])

length(unique(mytree1$tip.label))
length(unique(dados3$species))

dados3$strategy_3 <- as.factor(dados3$strategy_3)
summary(dados3$strategy_3)
length(dados3$continent[dados3$continent == "South_America"])

dados4 <- filter(dados3, strategy_3 == "strict_mig")
length(unique(dados4$species))

dados4 <- filter(dados3, strategy_3 == "partial_mig")
length(unique(dados4$species))

dados4 <- filter(dados3, strategy_3 == "resident")
length(unique(dados4$species))

plot.phylo(mytree1, type = "fan", show.tip.label = FALSE)

Fig <- filter(dados3, continent == "South_America")
Fig$link <- paste(Fig$Lineage_Name, Fig$species)

CommonLinks <- with(Fig, names(table(link)[table(link) > 4]))
CommonLinks <- Fig[Fig$link%in% CommonLinks, ]
CommonLinks$link <- factor(CommonLinks$link)
CommonLinks <- CommonLinks[,19]
CommonLinks <- distinct(CommonLinks)

###### Analysing cophylogeny #####

HostMatrix <- cophenetic(mytree1)
ParMatrix <- cophenetic(phylotree)

dados3 <- as.data.table(dados3)
HPMatrix <- dcast.data.table(dados3, species ~ Lineage)
HPMatrix <- column_to_rownames(HPMatrix, var = "species")
HPMatrix <- replace(HPMatrix, HPMatrix >= 1, 1)

HostMatrix <- HostMatrix[rownames(HPMatrix),rownames(HPMatrix)]
ParMatrix <- ParMatrix[colnames(HPMatrix),colnames(HPMatrix)]

PACo <- function (H.dist, P.dist, HP.bin)
{HP.bin <- which(HP.bin > 0, arr.in=TRUE)
H.PCo <- pcoa(H.dist, correction="cailliez")$vectors
P.PCo <- pcoa(P.dist, correction="cailliez")$vectors
H.PCo <- H.PCo[HP.bin[,1],]
P.PCo <- P.PCo[HP.bin[,2],]
list (H.PCo = H.PCo, P.PCo = P.PCo)}

PACo.fit <- PACo(HostMatrix, ParMatrix, HPMatrix)
HP.proc <- procrustes(PACo.fit$H.PCo, PACo.fit$P.PCo)

Hostx <- HP.proc$X
ParY <- HP.proc$Yrot
plot(Hostx)
points(ParY, pch=3)
arrows(ParY[,1], ParY[,2], Hostx[,1], Hostx[,2], length=0.12, angle=15,
       xpd=FALSE)

plot(Hostx)
points(ParY, pch=3)
arrows(ParY[,1], ParY[,2], Hostx[,1], Hostx[,2], length=0.12, angle=15,
       xpd=FALSE)


m2.obs <- HP.proc$ss
N.perm = 100000
P.value = 0

for (n in c(1:N.perm))
{if (NLinks <= nrow(HPMatrix) | NLinks <= ncol(HPMatrix))
{flag2 <- TRUE
while (flag2 == TRUE)
{HP.perm <- t(apply(HPMatrix,1,sample))
if(any(colSums(HP.perm) == NLinks)) flag2 <- TRUE
else flag2 <- FALSE}
} else { HP.perm <- t(apply(HPMatrix,1,sample))}
  PACo.perm <- PACo(HostMatrix, ParMatrix, HP.perm)
  m2.perm <- procrustes(PACo.perm$H.PCo, PACo.perm$P.PCo)$ss
  write(m2.perm, file="m2_perm.txt", sep="\t", append=TRUE)
  if (m2.perm <= m2.obs){P.value = P.value + 1}
}

P.value <- P.value/N.perm

NLinks = sum(HPMatrix)


HP.ones <- which(HPMatrix > 0, arr.in=TRUE)
SQres.jackn <- matrix(rep(NA, NLinks**2), NLinks)
colnames (SQres.jackn) <- paste(rownames(Hostx),rownames(ParY), sep="-")
t.critical = qt(0.975,NLinks-1)
for(i in c(1:NLinks))
{HP.ind <- HPMatrix
HP.ind[HP.ones[i,1],HP.ones[i,2]]=0
PACo.ind <- PACo(HostMatrix, ParMatrix, HP.ind)
Proc.ind <- procrustes(PACo.ind$H.PCo, PACo.ind$P.PCo)
res.Proc.ind <- c(residuals(Proc.ind))
res.Proc.ind <- append (res.Proc.ind, NA, after= i-1)
SQres.jackn [i, ] <- res.Proc.ind}
SQres.jackn <- SQres.jackn**2
SQres <- residuals (HP.proc)**2
          SQres.jackn <- SQres.jackn*(-(NLinks-1))
          SQres <- SQres*NLinks
          SQres.jackn <- t(apply(SQres.jackn, 1, "+", SQres))
          phi.mean <- apply(SQres.jackn, 2, mean, na.rm = TRUE)
          phi.UCI <- apply(SQres.jackn, 2, sd, na.rm = TRUE)
          phi.UCI <- phi.mean + t.critical * phi.UCI/sqrt(NLinks)
          

cophyvalues <- cbind.data.frame(cophy = phi.mean, sdcophy = phi.UCI, interaction = colnames(SQres.jackn))
cophyvalues <- separate(cophyvalues, interaction, sep = "-", into = c ("species", "parasite"))
cophyvalues <- rownames_to_column(cophyvalues, var = "interaction")
cophyvalues <- inner_join(cophyvalues, select(dados3, species, family, parasiteGenus, strategy_3, distance_4, distance_quanti_RES0,
                                              distance_quanti_ALL), by = "species")
cophyvalues <- distinct(cophyvalues)

cophyvalues$distance_4 <- factor(cophyvalues$distance_4, levels = c("resident", "short", "variable", "long"))

summary(cophyvalues$distance_4)

pat.bar <- barplot(phi.mean, names.arg = " ", space = 0.25, col="white", xlab=
                     "Host-parasite link", ylab= "Squared residuals",ylim=c(0, max(phi.UCI)),
                   cex.lab=1.2)
text(pat.bar, par("usr")[3] - 0.001, srt = 330, adj = 0, labels =
       colnames(SQres.jackn), xpd = TRUE, font = 1, cex=0.3)
arrows(pat.bar, phi.mean, pat.bar, phi.UCI, length= 0.05, angle=90)
abline(a=median(phi.mean), b=0, lty=2)

ggbar <- ggplot(cophyvalues, aes(x=factor(interaction), y = cophy)) + labs(x = "Host-parasite link", y = "Squared residuals") +
  geom_col(width=0.7, color = cophyvalues$distance_4) + theme_classic()

ggbar 


ggbar1 <- ggplot(cophyvalues, aes(x=factor(interaction), y = cophy)) + 
  geom_col(width=0.7, color = cophyvalues$distance_4) + theme_classic() +
  labs(x = "Host-parasite link", y = "Squared residuals") 
  
ggbar1

levels(cophyvalues$distance_4)

#write.csv2(cophyvalues, file = "cap5.csv")
cophyvalues[,c(1, 4:9)] <- sapply(cophyvalues[,c(1, 4:9)], as.factor)

boxplot(cophy ~ distance_4, data = cophyvalues)
plot(cophy ~ distance_quanti_ALL, data = cophyvalues)

cophyvalues[,c(2:3, 10:11)] <- sapply(cophyvalues[,c(2:3, 10:11)], as.numeric)
basicTrendline::trendline(cophyvalues$distance_quanti_ALL, cophyvalues$cophy)

#save.image("Cap4NoModels.RData")

##### Running PAco Model #####

D <- prepare_paco_data(HostMatrix, ParMatrix, HPMatrix)
D <- add_pcoord(D, correction = "none")

Paco <- PACo(D, nperm = 1000,
             seed = NA,
             method = "r0",
             symmetric = FALSE,
             proc.warnings = TRUE,
             shuffled = TRUE)
print(Paco$gof)

###### Baysian Models ######

inv.phylobird <- MCMCglmm::inverseA(mytree1, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylobird$Ainv)
rownames(A) <- rownames(inv.phylobird$Ainv)
class(A)

inv.phylo <- MCMCglmm::inverseA(phylotree, nodes = "TIPS", scale = FALSE)
B <- solve(inv.phylo$Ainv)
rownames(B) <- rownames(inv.phylo$Ainv)
class(B)

hist(cophyvalues$cophy, breaks = 100)
ggpubr::ggdensity(cophyvalues$cophy)
ggpubr::ggdensity(cophyvalues$distance_quanti_ALL)

variaveis <- bf(cophy ~ distance_4 + log1p(distance_quanti_ALL) + (1|interaction) + (1|gr(species, cov = A)),
                family = Gamma(link = "log"))

prior <- get_prior(variaveis, data = cophyvalues) #getting priors
prior


model <- brm(cophy ~ distance_4 + (1|interaction) + (1|gr(species, cov = A)),
             data = cophyvalues,
             family = Gamma(link = "log"), chains = 4,
             iter = 4000,
             data2 = list(A = A),
             prior = c(
               prior(student_t(3, 0, 2.5), "Intercept"),                                             
               prior(student_t(3, 0, 2.5), "sd"),
               prior(gamma(0.01, 0.01),"shape")
             ))
summary(model)

weights <- loo_model_weights(model, model1, model2, method = c("stacking", "pseudo-BMA"))
weights

summary(model)

a <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))

plot1 <- plot(conditional_effects(modelTesteDisSA), points = FALSE, theme = a)

plot2 <- plot1$distance_4 + labs(x = "Migratory distance category", y = "Cophylogenetic squared residuals ") #caption = "Plasmodium")
plot2


cophyPlas <- filter(cophyvalues, parasiteGenus == "Plasmodium")

modelPlas <- brm(cophy ~ distance_4 + (1|interaction) + (1|gr(species, cov = A)),
                      data = cophyPlas,
                      family = Gamma(link = "log"), chains = 4,
                      iter = 4000,
                      data2 = list(A = A),
                      prior = c(
                        prior(student_t(3, 8, 2.5), "Intercept"),                                             
                        prior(student_t(3, 0, 2.5), "sd"),
                        prior(gamma(0.01, 0.01),"shape")
                      ))
summary(modelPlas)

cophyHaem <- filter(cophyvalues, parasiteGenus == "Haemoproteus")

modelHaem <- brm(cophy ~ distance_4 + (1|interaction) + (1|gr(species, cov = A)),
                      data = cophyHaem,
                      family = Gamma(link = "log"), chains = 4,
                      iter = 4000,
                      data2 = list(A = A),
                      prior = c(
                        prior(student_t(3, 8, 2.5), "Intercept"),                                             
                        prior(student_t(3, 0, 2.5), "sd"),
                        prior(gamma(0.01, 0.01),"shape")
                      ))
summary(modelHaem)

##### South America Analyses #####

dados3SA <- filter(dados3, continent == "South_America")

dados4SA <- filter(dados3SA, strategy_3 == "strict_mig")
length(unique(dados4SA$species))

dados4SA <- filter(dados3SA, strategy_3 == "partial_mig")
length(unique(dados4SA$species))

dados4SA <- filter(dados3SA, strategy_3 == "resident")
length(unique(dados4SA$species))

matches2SA <-match(dados3SA$Lineage, phy3$tip.label)
matches2SA <-na.omit(matches2SA)
phylotreeSA <-drop.tip(phy3, phy3$tip.label[-matches2SA])

matches1SA <-match(dados3SA$Lineage, mytree1$tip.label)
matches1SA <-na.omit(matches1SA)
mytree1SA <-drop.tip(mytree1, mytree1$tip.label[-matches1SA])

HostMatrixSA <- cophenetic(mytree1SA)
ParMatrixSA <- cophenetic(phylotreeSA)

dados3SA <- as.data.table(dados3SA)
HPMatrixSA <- dcast.data.table(dados3SA, species ~ Lineage)
HPMatrixSA <- column_to_rownames(HPMatrixSA, var = "species")
HPMatrixSA <- replace(HPMatrixSA, HPMatrixSA >= 1, 1)

HostMatrixSA <- HostMatrixSA[rownames(HPMatrixSA),rownames(HPMatrixSA)]
ParMatrixSA <- ParMatrixSA[colnames(HPMatrixSA),colnames(HPMatrixSA)]

PACoSA <- function (H.dist, P.dist, HP.bin)
{HP.bin <- which(HP.bin > 0, arr.in=TRUE)
H.PCoSA <- pcoa(H.dist, correction="cailliez")$vectors
P.PCoSA <- pcoa(P.dist, correction="cailliez")$vectors
H.PCoSA <- H.PCoSA[HP.bin[,1],]
P.PCoSA <- P.PCoSA[HP.bin[,2],]
list (H.PCoSA = H.PCoSA, P.PCoSA = P.PCoSA)}

PACo.fitSA <- PACo(HostMatrixSA, ParMatrixSA, HPMatrixSA)
HP.procSA <- procrustes(PACo.fitSA$H.PCo, PACo.fitSA$P.PCo)

HostxSA <- HP.procSA$X
ParYSA <- HP.procSA$Yrot

m2.obsSA <- HP.procSA$ss
N.perm = 100000
NLinksSA = sum(HPMatrixSA)


HP.onesSA <- which(HPMatrixSA > 0, arr.in=TRUE)
SQres.jacknSA <- matrix(rep(NA, NLinksSA**2), NLinksSA)
colnames (SQres.jacknSA) <- paste(rownames(HostxSA),rownames(ParYSA), sep="-")
t.critical = qt(0.975,NLinksSA-1)
for(i in c(1:NLinksSA))
{HP.indSA <- HPMatrixSA
HP.indSA[HP.onesSA[i,1],HP.onesSA[i,2]]=0
PACo.indSA <- PACo(HostMatrixSA, ParMatrixSA, HP.indSA)
Proc.indSA <- procrustes(PACo.indSA$H.PCo, PACo.indSA$P.PCo)
res.Proc.indSA <- c(residuals(Proc.indSA))
res.Proc.indSA <- append (res.Proc.indSA, NA, after= i-1)
SQres.jacknSA [i, ] <- res.Proc.indSA}
SQres.jacknSA <- SQres.jacknSA**2
SQresSA <- residuals (HP.procSA)**2
SQres.jacknSA <- SQres.jacknSA*(-(NLinksSA-1))
SQresSA <- SQresSA*NLinksSA
SQres.jacknSA <- t(apply(SQres.jacknSA, 1, "+", SQresSA))
phi.meanSA <- apply(SQres.jacknSA, 2, mean, na.rm = TRUE)
phi.UCISA <- apply(SQres.jacknSA, 2, sd, na.rm = TRUE)
phi.UCISA <- phi.meanSA + t.critical * phi.UCISA/sqrt(NLinksSA)


cophyvaluesSA <- cbind.data.frame(cophy = phi.meanSA, sdcophy = phi.UCISA, interaction = colnames(SQres.jacknSA))
cophyvaluesSA <- separate(cophyvaluesSA, interaction, sep = "-", into = c ("species", "parasite"))
cophyvaluesSA <- rownames_to_column(cophyvaluesSA, var = "interaction")
cophyvaluesSA <- inner_join(cophyvaluesSA, select(dados3SA, species, family, parasiteGenus, strategy_3, distance_4, distance_quanti_RES0,
                                              distance_quanti_ALL), by = "species")
cophyvaluesSA <- distinct(cophyvaluesSA)


summary(cophyvaluesSA$strategy_3)

pat.barSA <- barplot(phi.meanSA, names.arg = " ", space = 0.25, col="white", xlab=
                     "Host-parasite link", ylab= "Squared residuals",ylim=c(0, max(phi.UCI)),
                   cex.lab=1.2)
text(pat.barSA, par("usr")[3] - 0.001, srt = 330, adj = 0, labels =
       colnames(SQres.jacknSA), xpd = TRUE, font = 1, cex=0.3)
arrows(pat.barSA, phi.meanSA, pat.barSA, phi.UCISA, length= 0.05, angle=90)
abline(a=median(phi.meanSA), b=0, lty=2)

ggbarSA <- ggplot(cophyvaluesSA, aes(x=factor(interaction), y = cophy)) + labs(x = "Host-parasite link", y = "Squared residuals") +
  geom_col(width=0.7, color = cophyvaluesSA$strategy_3) + theme_classic()

ggbarSA 


levels(cophyvaluesSA$strategy_3)

#write.csv2(cophyvalues, file = "cap5.csv")
cophyvaluesSA[,c(1, 4:7)] <- sapply(cophyvaluesSA[,c(1, 4:7)], as.factor)

boxplot(cophy ~ distance_4, data = cophyvaluesSA)
plot(cophy ~ distance_quanti_ALL, data = cophyvaluesSA)

cophyvaluesSA[,c(2:3, 10:11)] <- sapply(cophyvaluesSA[,c(2:3, 10:11)], as.numeric)
basicTrendline::trendline(cophyvaluesSA$distance_quanti_ALL, cophyvaluesSA$cophy)

###### Baysian Models ######

inv.phylobirdSA <- MCMCglmm::inverseA(mytree1SA, nodes = "TIPS", scale = TRUE)
ASA <- solve(inv.phylobirdSA$Ainv)
rownames(ASA) <- rownames(inv.phylobirdSA$Ainv)

hist(cophyvalues$cophy, breaks = 100)
ggpubr::ggdensity(cophyvalues$cophy)
ggpubr::ggdensity(cophyvalues$distance_quanti_ALL)

variaveis <- bf(cophy ~  distance_quanti_ALL + (1|interaction) + (1|gr(species, cov = A)),
                family = Gamma(link = "log"))

prior <- get_prior(variaveis, data = cophyvaluesSA) #getting priors
prior

cophyvaluesSA$distance_4 <- relevel(cophyvaluesSA$distance_4, ref = "resident")


modelSA <- brm(cophy ~ distance_4 + (1|interaction) + (1|gr(species, cov = ASA)),
             data = cophyvaluesSA,
             family = Gamma(link = "log"), chains = 4,
             iter = 4000,
             data2 = list(ASA = ASA),
             prior = c(
               prior(student_t(3, 8.6, 2.5), "Intercept"),                                             
               prior(student_t(3, 0, 2.5), "sd"),
               prior(gamma(0.01, 0.01),"shape")
             ))
summary(modelSA)

plot1SA <- plot(conditional_effects(modelSA), points = FALSE, theme = a)

plot2SA <- plot1SA$distance_4 + labs(x = "Migratory host category", y = "cophylogeny distance") #caption = "Plasmodium")
plot2SA


cophyPlasSA <- filter(cophyvaluesSA, parasiteGenus == "Plasmodium")

modelPlasSA <- brm(cophy ~ distance_4 + (1|interaction) + (1|gr(species, cov = ASA)),
                      data = cophyPlasSA,
                      family = Gamma(link = "log"), chains = 4,
                      iter = 4000,
                      data2 = list(ASA = ASA),
                      prior = c(
                        prior(student_t(3, 8.6, 2.5), "Intercept"),                                             
                        prior(student_t(3, 0, 2.5), "sd"),
                        prior(gamma(0.01, 0.01),"shape")
                      ))
summary(modelPlasSA)

cophyHaemSA <- filter(cophyvaluesSA, parasiteGenus == "Haemoproteus")

modelHaemSA <- brm(cophy ~ distance_4 + (1|interaction) + (1|gr(species, cov = A)),
                      data = cophyHaem,
                      family = Gamma(link = "log"), chains = 4,
                      iter = 4000,
                      data2 = list(A = A),
                      prior = c(
                        prior(student_t(3, 8, 2.5), "Intercept"),                                             
                        prior(student_t(3, 0, 2.5), "sd"),
                        prior(gamma(0.01, 0.01),"shape")
                      ))
summary(modelHaemSA)


##### Plotting Cophylogeny ######

library(ape)
library(phytools)
library(tidyverse)

load("C:/Users/danid/Documents/Lab Poulin/PhD/PhD/Dados/Capítulo 4/Tabelas/PlotCophy.RData")

cophyvalues[,c(1, 4:9)] <- sapply(cophyvalues[,c(1, 4:9)], as.factor)
cophyvaluesSA[,c(1, 4:9)] <- sapply(cophyvaluesSA[,c(1, 4:9)], as.factor)
cophyvalues[,c(2:3, 10:11)] <- sapply(cophyvalues[,c(2:3, 10:11)], as.numeric)
cophyvaluesSA[,c(2:3, 10:11)] <- sapply(cophyvaluesSA[,c(2:3, 10:11)], as.numeric)

Intdata <- cbind.data.frame(birds = cophyvalues$species, parasites <- cophyvalues$parasite)
class(Intdata)
Intdata <- distinct(Intdata)

cophylo <- cophylo(mytree1, phylotree, assoc = Intdata)

plot(mytree1, cex = 0.2)

plot(cophylo)
plot(cophylo,link.lwd=1, fsize = 0.5,
     link.lty="solid",link.col=make.transparent("blue",0.5))

cophyloplot(mytree1, phylotree, assoc = Intdata, show.tip.label = FALSE, type = "cladogram")

IntdataSA <- cbind.data.frame(birds = cophyvaluesSA$species, parasites <- cophyvaluesSA$parasite)
class(IntdataSA)
IntdataSA <- distinct(IntdataSA)

cophyloSA <- cophylo(mytree1SA, phylotreeSA, assoc = IntdataSA)

plot(cophyloSA,link.lwd=1, fsize = 0.30, #link.type="curved",
     link.lty="solid",link.col=make.transparent("red",
                                                0.5))
cophyloplot(mytree1SA, phylotreeSA, assoc = IntdataSA, show.tip.label = FALSE, type = "cladogram")


###### Plotting Map ######

newmap <- rworldmap::getMap(resolution = "high")
plot(newmap, asp = 1, col = "white")#, xlim = c(-60, -45), ylim = c(-40, 10))
dados3$continent <- as.factor(dados3$continent)
points(dados3$Longitude, dados3$Latitude, col= dados3$continent, cex = 1.0, pch = 20)
legend(x = "bottomleft", title = "Biomes", pch = 18, cex = 1.0, pt.cex = 1.6, bty = "n",
       #legend = c("Amazonia", "Andean Forest", "Caatinga", "Cerrado",
                  #"Grassland", "Atlantic Rain Forest", "Pantanal"),
       col = 1:7)
dev.off()


print(unique(Table1$biomes))

