## creating social anxiety disorder common metric
library(data.table)
library(mirt)
library(mirtCAT)
library(plyr)
library(BlandAltmanLeh)
library(weights)
library(foreign)
library(moments)
library(plink)
library(equate)
library(psych)
library(car)
library(ltm)
library(SNSequate)

## load dataset
comdata <- read.csv("H:/CRE/NHMRC Early Career Fellowship/Common metrics/papers/social anxiety equating/r-data2.csv")

totscores<-as.matrix(subset(comdata, select=c(SOCtot, LSAStot, SPINtot, SIAStot, SPStot, DSMSADtot, SIAS6tot, SPS6tot,GENSAD12tot)))

socbank <- subset(comdata, select=c(SADBNK1:SADBNK26))
lsas <- subset(comdata, select=c(LSASr1:LSASr24))
spin <- subset(comdata, select=c(SPIN1:SPIN17))
sias <- subset(comdata, select=c(SIAS1:SIAS20))
sias6<- subset(comdata, select=c(SIAS2, SIAS4, SIAS6, SIAS8, SIAS10, SIAS13))
sps <- subset(comdata, select=c(SPS1:SPS20))
sps6<- subset(comdata, select=c(SPS4, SPS7, SPS8, SPS15, SPS16, SPS17))
siassps12<-subset(comdata, select=c(SIAS2, SIAS4, SIAS6, SIAS8, SIAS10, SIAS13,SPS4, SPS7, SPS8, SPS15, SPS16, SPS17))
dsmsad <- subset(comdata, select=c(DSMSAD1:DSMSAD10))
minispin<-subset(comdata, select=c(SPIN6,SPIN9,SPIN15))

comb <- subset(comdata, select=c(SADBNK1:SADBNK26,LSASr1:LSASr24,SPIN1:SPIN17,SIAS1:SIAS20,SPS1:SPS20,DSMSAD1:DSMSAD10))
comb1 <- subset(comdata, select=c(SADBNK1:SADBNK26,LSASr1:LSASr24))
comb2 <- subset(comdata, select=c(SADBNK1:SADBNK26,SPIN1:SPIN17))
comb3 <- subset(comdata, select=c(SADBNK1:SADBNK26,SIAS1:SIAS20))
comb4 <- subset(comdata, select=c(SADBNK1:SADBNK26,SPS1:SPS20))
comb5 <- subset(comdata, select=c(SADBNK1:SADBNK26,DSMSAD1:DSMSAD10))
comb6<- subset(comdata, select=c(SADBNK1:SADBNK26,SIAS2, SIAS4, SIAS6, SIAS8, SIAS10, SIAS13))
comb7<- subset(comdata, select=c(SADBNK1:SADBNK26,SPS4, SPS7, SPS8, SPS15, SPS16, SPS17))
comb8<- subset(comdata, select=c(SADBNK1:SADBNK26,SIAS2, SIAS4, SIAS6, SIAS8, SIAS10, SIAS13,SPS4, SPS7, SPS8, SPS15, SPS16, SPS17))
comb9<- subset(comdata, select=c(SADBNK1:SADBNK26, SPIN6,SPIN9,SPIN15))

combuni<-mirt(comb,1, itemtype="graded", technical=list(NCYCLES=2000))
M2combuni<-M2(combuni)

##IRT item fit and residual matrices to check for local dependence
comb1uni<-mirt(comb1,1, itemtype="graded")
coef(comb1uni, simplify=TRUE)
ifit1<-itemfit(comb1uni)
residmat1<-M2(comb1uni,residmat=TRUE)


comb2uni<-mirt(comb2,1, itemtype="graded")
coef(comb2uni, simplify=TRUE)
ifit2<-itemfit(comb2uni)
residmat2<-M2(comb2uni,residmat=TRUE)

comb3uni<-mirt(comb3,1, itemtype="graded")
coef(comb3uni, simplify=TRUE)
ifit3<-itemfit(comb3uni)
residmat3<-M2(comb3uni,residmat=TRUE)

comb4uni<-mirt(comb4,1, itemtype="graded")
coef(comb4uni, simplify=TRUE)
ifit4<-itemfit(comb4uni)
residmat4<-M2(comb4uni,residmat=TRUE)

comb5uni<-mirt(comb5,1, itemtype="graded")
coef(comb5uni, simplify=TRUE)
ifit5<-itemfit(comb5uni)
residmat5<-M2(comb5uni,residmat=TRUE)

comb6uni<-mirt(comb6,1, itemtype="graded")
coef(comb6uni, simplify=TRUE)
ifit6<-itemfit(comb6uni)
residmat6<-M2(comb6uni,residmat=TRUE)



##fixed-anchor IRT calibration
sadpars <- read.csv("H:/CRE/NHMRC Early Career Fellowship/Common metrics/papers/social anxiety equating/sadpars.csv")
rownames(sadpars)<-sadpars[,1]
sadpars[,1]<- NULL
sadpars<-as.matrix(sadpars)

##SAD item bank and LSAS
## generating LSAS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb1,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parlsas<-mirt(comb1,1, itemtype="graded", pars=sv2)
coef(parlsas, IRTpars=TRUE, simplify=TRUE)
coeflsas<-coef(parlsas, as.data.frame=TRUE)
eap<-fscores(parlsas, method="EAP")
sum(eap[1:1052,])/1052

coef(parlsas, printSE=FALSE)

##calculate lsas cross-walk tables
svlsas <-mirt(lsas,1, itemtype="graded", pars="values")
sv2lsas<-svlsas
sv2lsas$value<-c(coeflsas[131:228,])

sv2lsas$est<- FALSE

lsaswalk<-mirt(lsas,1, itemtype="graded", pars=sv2lsas)
coef(lsaswalk)

infolsasfixed<-plot(lsaswalk, type="info")
infofixedlsas<-as.data.table(cbind(infolsasfixed$panel.args[[1]][[1]],infolsasfixed$panel.args[[1]][[2]]))

tcclsasfixed<-plot(lsaswalk, type="score")
scorefixedlsas<-as.data.table(cbind(tcclsasfixed$panel.args[[1]][[1]],tcclsasfixed$panel.args[[1]][[2]]))

lsaseap<-fscores(lsaswalk, method="EAP", full.scores.SE=TRUE)
sum(lsaseap[1:1052,1])/1052
lsaseapsum<-fscores(lsaswalk, method="EAPsum", full.scores = FALSE)

##SAD item bank and SPIN
## generating SPIN irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb2,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parspin<-mirt(comb2,1, itemtype="graded", pars=sv2)
coef(parspin, IRTpars=TRUE, simplify=TRUE)
coefspin<-coef(parspin, as.data.frame=TRUE)
eap<-fscores(parspin, method="EAP")
sum(eap[1:1052,])/1052

coef(parspin, printSE=TRUE)

##calculate spin cross-walk tables
svspin <-mirt(spin,1, itemtype="graded", pars="values")
sv2spin<-svspin
sv2spin$value<-c(coefspin[131:217,])

sv2spin$est<- FALSE

spinwalk<-mirt(spin,1, itemtype="graded", pars=sv2spin)
coef(spinwalk, simplify=TRUE)

infospinfixed<-plot(spinwalk, type="info")
infofixedspin<-as.data.table(cbind(infospinfixed$panel.args[[1]][[1]],infospinfixed$panel.args[[1]][[2]]))

tccspinfixed<-plot(spinwalk, type="score")
scorefixedspin<-as.data.table(cbind(tccspinfixed$panel.args[[1]][[1]],tccspinfixed$panel.args[[1]][[2]]))


spineap<-fscores(spinwalk, method="EAP", full.scores.SE=TRUE)
sum(spineap[1:1052,1])/1052
spineapsum<-fscores(spinwalk, method="EAPsum", full.scores = FALSE)


##SAD item bank and mini-SPIN
## generating mini-SPIN irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb9,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parmspin<-mirt(comb9,1, itemtype="graded", pars=sv2)
coef(parmspin, IRTpars=TRUE, simplify=TRUE)
coefmspin<-coef(parmspin, as.data.frame=TRUE)
eap<-fscores(parmspin, method="EAP")
sum(eap[1:1052,])/1052

##calculate spin cross-walk tables
svmspin <-mirt(minispin,1, itemtype="graded", pars="values")
sv2mspin<-svmspin
sv2mspin$value<-c(coefmspin[131:147,])

sv2mspin$est<- FALSE

mspinwalk<-mirt(minispin,1, itemtype="graded", pars=sv2mspin)
coef(mspinwalk, simplify=TRUE)

infomspinfixed<-plot(mspinwalk, type="info")
infofixedmspin<-as.data.table(cbind(infomspinfixed$panel.args[[1]][[1]],infomspinfixed$panel.args[[1]][[2]]))

tccmspinfixed<-plot(mspinwalk, type="score")
scorefixedmspin<-as.data.table(cbind(tccmspinfixed$panel.args[[1]][[1]],tccmspinfixed$panel.args[[1]][[2]]))


mspineap<-fscores(mspinwalk, method="EAP", full.scores.SE=TRUE)
sum(mspineap[1:1052,1])/1052
mspineapsum<-fscores(mspinwalk, method="EAPsum", full.scores = FALSE)


##SAD item bank and SIAS
## generating SIAS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb3,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parsias<-mirt(comb3,1, itemtype="graded", pars=sv2)
coef(parsias, IRTpars=TRUE, simplify=TRUE)
coefsias<-coef(parsias, as.data.frame=TRUE)
eap<-fscores(parsias, method="EAP")
sum(eap[1:1052,])/1052

coef(parsias, printSE=TRUE)

##calculate sias cross-walk tables
svsias <-mirt(sias,1, itemtype="graded", pars="values")
sv2sias<-svsias
sv2sias$value<-c(coefsias[131:232,])

sv2sias$est<- FALSE

siaswalk<-mirt(sias,1, itemtype="graded", pars=sv2sias)
coef(siaswalk, simplify=TRUE)

infosiasfixed<-plot(siaswalk, type="info")
infofixedsias<-as.data.table(cbind(infosiasfixed$panel.args[[1]][[1]],infosiasfixed$panel.args[[1]][[2]]))

tccsiasfixed<-plot(siaswalk, type="score")
scorefixedsias<-as.data.table(cbind(tccsiasfixed$panel.args[[1]][[1]],tccsiasfixed$panel.args[[1]][[2]]))


siaseap<-fscores(siaswalk, method="EAP", full.scores.SE=TRUE)
sum(siaseap[1:1052,])/1052
siaseapsum<-fscores(siaswalk, method="EAPsum", full.scores = FALSE)

##SAD item bank and SPS
## generating SPS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb4,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parsps<-mirt(comb4,1, itemtype="graded", pars=sv2)
coef(parsps, IRTpars=TRUE, simplify=TRUE)
coefsps<-coef(parsps, as.data.frame=TRUE)
eap<-fscores(parsps, method="EAP")
sum(eap[1:1052,])/1052

coef(parsps, printSE=TRUE)

##calculate sps cross-walk tables
svsps <-mirt(sps,1, itemtype="graded", pars="values")
sv2sps<-svsps
sv2sps$value<-c(coefsps[131:232,])

sv2sps$est<- FALSE

spswalk<-mirt(sps,1, itemtype="graded", pars=sv2sps)
coef(spswalk, simplify=TRUE)

infospsfixed<-plot(spswalk, type="info")
infofixedsps<-as.data.table(cbind(infospsfixed$panel.args[[1]][[1]],infospsfixed$panel.args[[1]][[2]]))

tccspsfixed<-plot(spswalk, type="score")
scorefixedsps<-as.data.table(cbind(tccspsfixed$panel.args[[1]][[1]],tccspsfixed$panel.args[[1]][[2]]))


spseap<-fscores(spswalk, method="EAP", full.scores.SE=TRUE)
sum(spseap[1:1052,])/1052
spseapsum<-fscores(spswalk, method="EAPsum", full.scores = FALSE)

##SAD item bank and dsmsad
## generating dsmsad irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb5,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

pardsm<-mirt(comb5,1, itemtype="graded", pars=sv2)
coef(pardsm, IRTpars=TRUE, simplify=TRUE)
coefdsm<-coef(pardsm, as.data.frame=TRUE)
eap<-fscores(pardsm, method="EAP")
sum(eap[1:1052,])/1052

coef(pardsm, printSE=TRUE)

##calculate dsmsad cross-walk tables
svdsm <-mirt(dsmsad,1, itemtype="graded", pars="values")
sv2dsm<-svdsm
sv2dsm$value<-c(coefdsm[131:182,])

sv2dsm$est<- FALSE

dsmwalk<-mirt(dsmsad,1, itemtype="graded", pars=sv2dsm)
coef(dsmwalk, simplify=TRUE)

infodsmfixed<-plot(dsmwalk, type="info")
infofixeddsm<-as.data.table(cbind(infodsmfixed$panel.args[[1]][[1]],infodsmfixed$panel.args[[1]][[2]]))

tccdsmfixed<-plot(dsmwalk, type="score")
scorefixeddsm<-as.data.table(cbind(tccdsmfixed$panel.args[[1]][[1]],tccdsmfixed$panel.args[[1]][[2]]))


dsmeap<-fscores(dsmwalk, method="EAP", full.scores.SE=TRUE)
sum(dsmeap[1:1052,])/1052
dsmeapsum<-fscores(dsmwalk, method="EAPsum", full.scores = FALSE)

##SAD item bank and SIAS-6
## generating SIAS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb6,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parsias6<-mirt(comb6,1, itemtype="graded", pars=sv2)
coef(parsias6, IRTpars=TRUE, simplify=TRUE)
coefsias6<-coef(parsias6, as.data.frame=TRUE)
eap<-fscores(parsias6, method="EAP")
sum(eap[1:1052,])/1052

##calculate sias6 cross-walk tables
svsias6 <-mirt(sias6,1, itemtype="graded", pars="values")
sv2sias6<-svsias6
sv2sias6$value<-c(coefsias6[131:162,])

#testing if taking parameters from full scale improves fit results
sv2sias6$value<-c(1.75184519,
                  -0.47341958,
                  -1.91630507,
                  -3.00843537,
                  -4.27000768,
                  2.18088004,
                  0.12079298,
                  -1.83200034,
                  -3.29545988,
                  -4.73504618,
                  2.14457292,
                  -0.2106793,
                  -2.13581268,
                  -3.482685,
                  -5.07212058,
                  1.70280061,
                  -0.80462497,
                  -2.35599203,
                  -3.43462696,
                  -4.85423246,
                  2.56628114,
                  -0.12004936,
                  -2.45057471,
                  -4.01169328,
                  -5.69183403,
                  1.0422526,
                  0.2079037,
                  -1.08527414,
                  -2.26428883,
                  -3.58159414,
                  0,
                  1)

sv2sias6$est<- FALSE

sias6walk<-mirt(sias6,1, itemtype="graded", pars=sv2sias6)
coef(sias6walk, simplify=TRUE)

infosias6fixed<-plot(sias6walk, type="info")
infofixedsias6<-as.data.table(cbind(infosias6fixed$panel.args[[1]][[1]],infosias6fixed$panel.args[[1]][[2]]))

tccsias6fixed<-plot(sias6walk, type="score")
scorefixedsias6<-as.data.table(cbind(tccsias6fixed$panel.args[[1]][[1]],tccsias6fixed$panel.args[[1]][[2]]))


sias6eap<-fscores(sias6walk, method="EAP", full.scores.SE=TRUE)
sum(sias6eap[1:1052,1])/1052
sias6eapsum<-fscores(sias6walk, method="EAPsum", full.scores = FALSE)

##SAD item bank and SPS-6
## generating sps6 irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb7,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parsps6<-mirt(comb7,1, itemtype="graded", pars=sv2)
coef(parsps6, IRTpars=TRUE, simplify=TRUE)
coefsps6<-coef(parsps6, as.data.frame=TRUE)
eap<-fscores(parsps6, method="EAP")
sum(eap[1:1052,])/1052

##calculate sps6 cross-walk tables
svsps6 <-mirt(sps6,1, itemtype="graded", pars="values")
sv2sps6<-svsps6
sv2sps6$value<-c(coefsps6[131:162,])

sv2sps6$est<- FALSE

sps6walk<-mirt(sps6,1, itemtype="graded", pars=sv2sps6)
coef(sps6walk, simplify=TRUE)

infosps6fixed<-plot(sps6walk, type="info")
infofixedsps6<-as.data.table(cbind(infosps6fixed$panel.args[[1]][[1]],infosps6fixed$panel.args[[1]][[2]]))

tccsps6fixed<-plot(sps6walk, type="score")
scorefixedsps6<-as.data.table(cbind(tccsps6fixed$panel.args[[1]][[1]],tccsps6fixed$panel.args[[1]][[2]]))


sps6eap<-fscores(sps6walk, method="EAP", full.scores.SE=TRUE)
sum(sps6eap[1:1052,1])/1052
sps6eapsum<-fscores(sps6walk, method="EAPsum", full.scores = FALSE)

##SAD item bank and SIAS-SPS-12
## generating SIAS-SPS-12 irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb8,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parsps12<-mirt(comb8,1, itemtype="graded", pars=sv2)
coef(parsps12, IRTpars=TRUE, simplify=TRUE)
coefsps12<-coef(parsps12, as.data.frame=TRUE)
eap<-fscores(parsps12, method="EAP")
sum(eap[1:1052,])/1052

##calculate sps6 cross-walk tables
svsps12 <-mirt(siassps12,1, itemtype="graded", pars="values")
sv2sps12<-svsps12
sv2sps12$value<-c(coefsps12[131:192,])

sv2sps12$est<- FALSE

sps12walk<-mirt(siassps12,1, itemtype="graded", pars=sv2sps12)
coef(sps12walk, simplify=TRUE)

infosps12fixed<-plot(sps12walk, type="info")
infofixedsps12<-as.data.table(cbind(infosps12fixed$panel.args[[1]][[1]],infosps12fixed$panel.args[[1]][[2]]))

tccsps12fixed<-plot(sps12walk, type="score")
scorefixedsps12<-as.data.table(cbind(tccsps12fixed$panel.args[[1]][[1]],tccsps12fixed$panel.args[[1]][[2]]))


sps12eap<-fscores(sps12walk, method="EAP", full.scores.SE=TRUE)
sum(sps12eap[1:1052,1])/1052
sps12eapsum<-fscores(sps12walk, method="EAPsum", full.scores = FALSE)

##calculate socbank cross-walk tables
svsocb <-mirt(socbank,1, itemtype="graded", pars="values")
sv2socb<-svsocb
sv2socb$value[1:length(sadvals)] = sadvals

sv2socb$est<- FALSE

socbwalk<-mirt(socbank,1, itemtype="graded", pars=sv2socb)
coef(socbwalk, simplify=TRUE)

infosbankfixed<-plot(socbwalk, type="info")
infofixedsbank<-as.data.table(cbind(infosbankfixed$panel.args[[1]][[1]],infosbankfixed$panel.args[[1]][[2]]))


socbeap<-fscores(socbwalk, method="EAP", full.scores.SE=TRUE)
sum(socbeap[1:1052,1])/1052
socbeapsum<-fscores(socbwalk, method="EAPsum", full.scores = FALSE)


##fixed parameter crosswalk for all items included in pool
sv <-mirt(comb,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parall<-mirt(comb,1, itemtype="graded", pars=sv2)
coef(parall, IRTpars=TRUE, simplify=TRUE)
coefall<-coef(parall, as.data.frame=TRUE)
eap<-fscores(parall, method="EAP")
sum(eap[1:1052,])/1052

id<-1:nrow(coefall)
coefall<-cbind(coefall,id)

svlsas <-mirt(lsas,1, itemtype="graded", pars="values")
sv2lsas<-svlsas
sv2lsas$value<-c(coefall[131:226,1],coefall[562:563,1])

sv2lsas$est<- FALSE

lsaswalk<-mirt(lsas,1, itemtype="graded", pars=sv2lsas)
coef(lsaswalk, simplify=TRUE)

tcclsasfixed<-plot(lsaswalk, type="score")
scorefixedlsas<-as.data.table(cbind(tcclsasfixed$panel.args[[1]][[1]],tcclsasfixed$panel.args[[1]][[2]]))

lsaseap<-fscores(lsaswalk, method="EAP", full.scores.SE=TRUE)
sum(lsaseap[1:1052,1])/1052
lsaseapsum<-fscores(lsaswalk, method="EAPsum", full.scores = FALSE)

cor(socbeap[,1],lsaseap[,1])
mean(socbeap[,1]-lsaseap[,1])
SD(socbeap[,1]-lsaseap[,1])
MSD <- mean(((socbeap[,1] - lsaseap[,1])^2))
(sqrt(MSD))
ICC(as.matrix(cbind(socbeap[,1],lsaseap[,1])))

svspin <-mirt(spin,1, itemtype="graded", pars="values")
sv2spin<-svspin
sv2spin$value<-c(coefall[227:311,1],coefall[562:563,1])

sv2spin$est<- FALSE

spinwalk<-mirt(spin,1, itemtype="graded", pars=sv2spin)
coef(spinwalk, simplify=TRUE)

tccspinfixed<-plot(spinwalk, type="score")
scorefixedspin<-as.data.table(cbind(tccspinfixed$panel.args[[1]][[1]],tccspinfixed$panel.args[[1]][[2]]))

spineap<-fscores(spinwalk, method="EAP", full.scores.SE=TRUE)
sum(spineap[1:1052,1])/1052
spineapsum<-fscores(spinwalk, method="EAPsum", full.scores = FALSE)

cor(socbeap[,1],spineap[,1])
mean(socbeap[,1]-spineap[,1])
SD(socbeap[,1]-spineap[,1])
MSD <- mean(((socbeap[,1] - spineap[,1])^2))
(sqrt(MSD))
ICC(as.matrix(cbind(socbeap[,1],spineap[,1])))


##plink IRT separate calibration equating
##SAD item bank and LSAS
##generating sadbank coefficients
sv <-mirt(comb1,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE
mirtlsas<-mirt(comb1,1, itemtype="graded", pars=sv2)
mirtceof<-coef(mirtlsas, IRTpars=TRUE, simplify=TRUE)
sadparsirt<-mirtceof$items[1:26,]

## generating separeately calibrated items for sadbank and lsas
comb1uni<-mirt(comb1,1, itemtype="graded")
sadparlsas<-coef(comb1uni, simplify=TRUE, IRTpars=TRUE)
sadparlsas<-as.matrix(sadparlsas$items)
##specific number of categories for both item parameters to be linked
catsad<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
catlsas<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)
##specific the model to be used for both sets of item paramters
pmsad<-as.poly.mod(26,c("grm"),list(c(1:26)))
pmlsas<-as.poly.mod(50,c("grm"),list(c(1:50)))
##indicate the number of common items between the two sets of item parameters
commonlsas<-matrix(c(1:26, 1:26),26,2)
##combine all info into a single object
parslsas<-as.irt.pars(x=list(sadparsirt,sadparlsas), common=list(commonlsas), cat=list(catsad, catlsas), poly.mod=list(pmsad, pmlsas),grp.names=c("sadbank","lsas"))

##stocking lord
plinklsas<-plink(parslsas, rescale="SL")
summary(plinklsas)
reparlsasSL<-link.pars(plinklsas)
reparlsasSL<-reparlsasSL$lsas[27:50,1:4]
d1=-1*(reparlsasSL[,1]*reparlsasSL[,2])
d2=-1*(reparlsasSL[,1]*reparlsasSL[,3])
d3=-1*(reparlsasSL[,1]*reparlsasSL[,4])
reparlsasSL2<-cbind(reparlsasSL[,1],d1,d2,d3)
colnames(reparlsasSL2)<-c("a1","d1","d2","d3")
object_reparlsasSL<-generate.mirt_object(reparlsasSL2, "graded")
tcclsasSF<-plot(object_reparlsasSL, type="score")
tcclsasSF
scoreSFlsas<-as.data.table(cbind(tcclsasSF$panel.args[[1]][[1]],tcclsasSF$panel.args[[1]][[2]]))
#mean/mean
plinklsas<-plink(parslsas, rescale="MM")
summary(plinklsas)
reparlsasMM<-link.pars(plinklsas)
reparlsasMM<-reparlsasMM$lsas[27:50,1:4]
d1=-1*(reparlsasMM[,1]*reparlsasMM[,2])
d2=-1*(reparlsasMM[,1]*reparlsasMM[,3])
d3=-1*(reparlsasMM[,1]*reparlsasMM[,4])
reparlsasMM2<-cbind(reparlsasMM[,1],d1,d2,d3)
colnames(reparlsasMM2)<-c("a1","d1","d2","d3")
object_reparlsasMM<-generate.mirt_object(reparlsasMM2, "graded")
tcclsasMM<-plot(object_reparlsasMM, type="score")
scoreMMlsas<-as.data.table(cbind(tcclsasMM$panel.args[[1]][[1]],tcclsasMM$panel.args[[1]][[2]]))
#mean/sigma
plinklsas<-plink(parslsas, rescale="MS")
summary(plinklsas)
reparlsasMS<-link.pars(plinklsas)
reparlsasMS<-reparlsasMS$lsas[27:50,1:4]
d1=-1*(reparlsasMS[,1]*reparlsasMS[,2])
d2=-1*(reparlsasMS[,1]*reparlsasMS[,3])
d3=-1*(reparlsasMS[,1]*reparlsasMS[,4])
reparlsasMS2<-cbind(reparlsasMS[,1],d1,d2,d3)
colnames(reparlsasMS2)<-c("a1","d1","d2","d3")
object_reparlsasMS<-generate.mirt_object(reparlsasMS2, "graded")
tcclsasMS<-plot(object_reparlsasMS, type="score")
scoreMSlsas<-as.data.table(cbind(tcclsasMS$panel.args[[1]][[1]],tcclsasMS$panel.args[[1]][[2]]))
#haebara
plinklsas<-plink(parslsas, rescale="HB")
summary(plinklsas)
reparlsasHB<-link.pars(plinklsas)
reparlsasHB<-reparlsasHB$lsas[27:50,1:4]
d1=-1*(reparlsasHB[,1]*reparlsasHB[,2])
d2=-1*(reparlsasHB[,1]*reparlsasHB[,3])
d3=-1*(reparlsasHB[,1]*reparlsasHB[,4])
reparlsasHB2<-cbind(reparlsasHB[,1],d1,d2,d3)
colnames(reparlsasHB2)<-c("a1","d1","d2","d3")
object_reparlsasHB<-generate.mirt_object(reparlsasHB2, "graded")
tcclsasHB<-plot(object_reparlsasHB, type="score")
scoreHBlsas<-as.data.table(cbind(tcclsasHB$panel.args[[1]][[1]],tcclsasHB$panel.args[[1]][[2]]))


## generating separeately calibrated items for sadbank and SPIN
comb2uni<-mirt(comb2,1, itemtype="graded")
sadparspin<-coef(comb2uni, simplify=TRUE, IRTpars=TRUE)
sadparspin<-as.matrix(sadparspin$items)
##specific number of categories for both item parameters to be linked
catsad<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
catspin<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
##specific the model to be used for both sets of item paramters
pmsad<-as.poly.mod(26,c("grm"),list(c(1:26)))
pmspin<-as.poly.mod(43,c("grm"),list(c(1:43)))
##indicate the number of common items between the two sets of item parameters
commonspin<-matrix(c(1:26, 1:26),26,2)
##combine all info into a single object
parsspin<-as.irt.pars(x=list(sadparsirt,sadparspin), common=list(commonspin), cat=list(catsad, catspin), poly.mod=list(pmsad, pmspin),grp.names=c("sadbank","spin"))

##stocking lord
plinkspin<-plink(parsspin, rescale="SL")
summary(plinkspin)
reparspinSL<-link.pars(plinkspin)
reparspinSL<-reparspinSL$spin[27:43,1:5]
d1=-1*(reparspinSL[,1]*reparspinSL[,2])
d2=-1*(reparspinSL[,1]*reparspinSL[,3])
d3=-1*(reparspinSL[,1]*reparspinSL[,4])
d4=-1*(reparspinSL[,1]*reparspinSL[,5])
reparspinSL2<-cbind(reparspinSL[,1],d1,d2,d3,d4)
colnames(reparspinSL2)<-c("a1","d1","d2","d3","d4")
object_reparspinSL<-generate.mirt_object(reparspinSL2, "graded")
tccspinSF<-plot(object_reparspinSL, type="score")
scoreSFspin<-as.data.table(cbind(tccspinSF$panel.args[[1]][[1]],tccspinSF$panel.args[[1]][[2]]))
#mean/mean
plinkspin<-plink(parsspin, rescale="MM")
summary(plinkspin)
reparspinMM<-link.pars(plinkspin)
reparspinMM<-reparspinMM$spin[27:43,1:5]
d1=-1*(reparspinMM[,1]*reparspinMM[,2])
d2=-1*(reparspinMM[,1]*reparspinMM[,3])
d3=-1*(reparspinMM[,1]*reparspinMM[,4])
d4=-1*(reparspinMM[,1]*reparspinMM[,5])
reparspinMM2<-cbind(reparspinMM[,1],d1,d2,d3,d4)
colnames(reparspinMM2)<-c("a1","d1","d2","d3","d4")
object_reparspinMM<-generate.mirt_object(reparspinMM2, "graded")
tccspinMM<-plot(object_reparspinMM, type="score")
scoreMMspin<-as.data.table(cbind(tccspinMM$panel.args[[1]][[1]],tccspinMM$panel.args[[1]][[2]]))
#mean/sigma
plinkspin<-plink(parsspin, rescale="MS")
summary(plinkspin)
reparspinMS<-link.pars(plinkspin)
reparspinMS<-reparspinMS$spin[27:43,1:5]
d1=-1*(reparspinMS[,1]*reparspinMS[,2])
d2=-1*(reparspinMS[,1]*reparspinMS[,3])
d3=-1*(reparspinMS[,1]*reparspinMS[,4])
d4=-1*(reparspinMS[,1]*reparspinMS[,5])
reparspinMS2<-cbind(reparspinMS[,1],d1,d2,d3,d4)
colnames(reparspinMS2)<-c("a1","d1","d2","d3","d4")
object_reparspinMS<-generate.mirt_object(reparspinMS2, "graded")
tccspinMS<-plot(object_reparspinMS, type="score")
scoreMSspin<-as.data.table(cbind(tccspinMS$panel.args[[1]][[1]],tccspinMS$panel.args[[1]][[2]]))
#haebara
plinkspin<-plink(parsspin, rescale="HB")
summary(plinkspin)
reparspinHB<-link.pars(plinkspin)
reparspinHB<-reparspinHB$spin[27:43,1:5]
d1=-1*(reparspinHB[,1]*reparspinHB[,2])
d2=-1*(reparspinHB[,1]*reparspinHB[,3])
d3=-1*(reparspinHB[,1]*reparspinHB[,4])
d4=-1*(reparspinHB[,1]*reparspinHB[,5])
reparspinHB2<-cbind(reparspinHB[,1],d1,d2,d3,d4)
colnames(reparspinHB2)<-c("a1","d1","d2","d3","d4")
object_reparspinHB<-generate.mirt_object(reparspinHB2, "graded")
tccspinHB<-plot(object_reparspinHB, type="score")
scoreHBspin<-as.data.table(cbind(tccspinHB$panel.args[[1]][[1]],tccspinHB$panel.args[[1]][[2]]))

## generating separeately calibrated items for sadbank and SIAS
comb3uni<-mirt(comb3,1, itemtype="graded")
sadparsias<-coef(comb3uni, simplify=TRUE, IRTpars=TRUE)
sadparsias<-as.matrix(sadparsias$items)
##specific number of categories for both item parameters to be linked
catsad<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
catsias<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
##specific the model to be used for both sets of item paramters
pmsad<-as.poly.mod(26,c("grm"),list(c(1:26)))
pmsias<-as.poly.mod(46,c("grm"),list(c(1:46)))
##indicate the number of common items between the two sets of item parameters
commonsias<-matrix(c(1:26, 1:26),26,2)
##combine all info into a single object
parssias<-as.irt.pars(x=list(sadparsirt,sadparsias), common=list(commonsias), cat=list(catsad, catsias), poly.mod=list(pmsad, pmsias),grp.names=c("sadbank","sias"))

##stocking lord
plinksias<-plink(parssias, rescale="SL")
summary(plinksias)
reparsiasSL<-link.pars(plinksias)
reparsiasSL<-reparsiasSL$sias[27:46,1:5]
d1=-1*(reparsiasSL[,1]*reparsiasSL[,2])
d2=-1*(reparsiasSL[,1]*reparsiasSL[,3])
d3=-1*(reparsiasSL[,1]*reparsiasSL[,4])
d4=-1*(reparsiasSL[,1]*reparsiasSL[,5])
reparsiasSL2<-cbind(reparsiasSL[,1],d1,d2,d3,d4)
colnames(reparsiasSL2)<-c("a1","d1","d2","d3","d4")
object_reparsiasSL<-generate.mirt_object(reparsiasSL2, "graded")
tccsiasSL<-plot(object_reparsiasSL, type="score")
scoreSLsias<-as.data.table(cbind(tccsiasSL$panel.args[[1]][[1]],tccsiasSL$panel.args[[1]][[2]]))
##mean/mean
plinksias<-plink(parssias, rescale="MM")
summary(plinksias)
reparsiasMM<-link.pars(plinksias)
reparsiasMM<-reparsiasMM$sias[27:46,1:5]
d1=-1*(reparsiasMM[,1]*reparsiasMM[,2])
d2=-1*(reparsiasMM[,1]*reparsiasMM[,3])
d3=-1*(reparsiasMM[,1]*reparsiasMM[,4])
d4=-1*(reparsiasMM[,1]*reparsiasMM[,5])
reparsiasMM2<-cbind(reparsiasMM[,1],d1,d2,d3,d4)
colnames(reparsiasMM2)<-c("a1","d1","d2","d3","d4")
object_reparsiasMM<-generate.mirt_object(reparsiasMM2, "graded")
tccsiasMM<-plot(object_reparsiasMM, type="score")
scoreMMsias<-as.data.table(cbind(tccsiasMM$panel.args[[1]][[1]],tccsiasMM$panel.args[[1]][[2]]))
##mean/sigma
plinksias<-plink(parssias, rescale="MS")
summary(plinksias)
reparsiasMS<-link.pars(plinksias)
reparsiasMS<-reparsiasMS$sias[27:46,1:5]
d1=-1*(reparsiasMS[,1]*reparsiasMS[,2])
d2=-1*(reparsiasMS[,1]*reparsiasMS[,3])
d3=-1*(reparsiasMS[,1]*reparsiasMS[,4])
d4=-1*(reparsiasMS[,1]*reparsiasMS[,5])
reparsiasMS2<-cbind(reparsiasMS[,1],d1,d2,d3,d4)
colnames(reparsiasMS2)<-c("a1","d1","d2","d3","d4")
object_reparsiasMS<-generate.mirt_object(reparsiasMS2, "graded")
tccsiasMS<-plot(object_reparsiasMS, type="score")
scoreMSsias<-as.data.table(cbind(tccsiasMS$panel.args[[1]][[1]],tccsiasMS$panel.args[[1]][[2]]))
##haebara
plinksias<-plink(parssias, rescale="HB")
summary(plinksias)
reparsiasHB<-link.pars(plinksias)
reparsiasHB<-reparsiasHB$sias[27:46,1:5]
d1=-1*(reparsiasHB[,1]*reparsiasHB[,2])
d2=-1*(reparsiasHB[,1]*reparsiasHB[,3])
d3=-1*(reparsiasHB[,1]*reparsiasHB[,4])
d4=-1*(reparsiasHB[,1]*reparsiasHB[,5])
reparsiasHB2<-cbind(reparsiasHB[,1],d1,d2,d3,d4)
colnames(reparsiasHB2)<-c("a1","d1","d2","d3","d4")
object_reparsiasHB<-generate.mirt_object(reparsiasHB2, "graded")
tccsiasHB<-plot(object_reparsiasHB, type="score")
scoreHBsias<-as.data.table(cbind(tccsiasHB$panel.args[[1]][[1]],tccsiasHB$panel.args[[1]][[2]]))


## generating separeately calibrated items for sadbank and SPS
comb4uni<-mirt(comb4,1, itemtype="graded", technical=list(NCYCLES=1000))
sadparsps<-coef(comb4uni, simplify=TRUE, IRTpars=TRUE)
sadparsps<-as.matrix(sadparsps$items)
##specific number of categories for both item parameters to be linked
catsad<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
catsps<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
##specific the model to be used for both sets of item paramters
pmsad<-as.poly.mod(26,c("grm"),list(c(1:26)))
pmsps<-as.poly.mod(46,c("grm"),list(c(1:46)))
##indicate the number of common items between the two sets of item parameters
commonsps<-matrix(c(1:26, 1:26),26,2)
##combine all info into a single object
parssps<-as.irt.pars(x=list(sadparsirt,sadparsps), common=list(commonsps), cat=list(catsad, catsps), poly.mod=list(pmsad, pmsps),grp.names=c("sadbank","sps"))

##stocking lord
plinksps<-plink(parssps, rescale="SL")
summary(plinksps)
reparspsSL<-link.pars(plinksps)
reparspsSL<-reparspsSL$sps[27:46,1:5]
d1=-1*(reparspsSL[,1]*reparspsSL[,2])
d2=-1*(reparspsSL[,1]*reparspsSL[,3])
d3=-1*(reparspsSL[,1]*reparspsSL[,4])
d4=-1*(reparspsSL[,1]*reparspsSL[,5])
reparspsSL2<-cbind(reparspsSL[,1],d1,d2,d3,d4)
colnames(reparspsSL2)<-c("a1","d1","d2","d3","d4")
object_reparspsSL<-generate.mirt_object(reparspsSL2, "graded")
tccspsSL<-plot(object_reparspsSL, type="score")
scoreSLsps<-as.data.table(cbind(tccspsSL$panel.args[[1]][[1]],tccspsSL$panel.args[[1]][[2]]))
##mean mean
plinksps<-plink(parssps, rescale="MM")
summary(plinksps)
reparspsMM<-link.pars(plinksps)
reparspsMM<-reparspsMM$sps[27:46,1:5]
d1=-1*(reparspsMM[,1]*reparspsMM[,2])
d2=-1*(reparspsMM[,1]*reparspsMM[,3])
d3=-1*(reparspsMM[,1]*reparspsMM[,4])
d4=-1*(reparspsMM[,1]*reparspsMM[,5])
reparspsMM2<-cbind(reparspsMM[,1],d1,d2,d3,d4)
colnames(reparspsMM2)<-c("a1","d1","d2","d3","d4")
object_reparspsMM<-generate.mirt_object(reparspsMM2, "graded")
tccspsMM<-plot(object_reparspsMM, type="score")
scoreMMsps<-as.data.table(cbind(tccspsMM$panel.args[[1]][[1]],tccspsMM$panel.args[[1]][[2]]))
spsmmfscores<-fscores(object_reparspsMM, response.pattern=sps,method="EAP")
##mean sigma
plinksps<-plink(parssps, rescale="MS")
summary(plinksps)
reparspsMS<-link.pars(plinksps)
reparspsMS<-reparspsMS$sps[27:46,1:5]
d1=-1*(reparspsMS[,1]*reparspsMS[,2])
d2=-1*(reparspsMS[,1]*reparspsMS[,3])
d3=-1*(reparspsMS[,1]*reparspsMS[,4])
d4=-1*(reparspsMS[,1]*reparspsMS[,5])
reparspsMS2<-cbind(reparspsMS[,1],d1,d2,d3,d4)
colnames(reparspsMS2)<-c("a1","d1","d2","d3","d4")
object_reparspsMS<-generate.mirt_object(reparspsMS2, "graded")
tccspsMS<-plot(object_reparspsMS, type="score")
scoreMSsps<-as.data.table(cbind(tccspsMS$panel.args[[1]][[1]],tccspsMS$panel.args[[1]][[2]]))
##haebara
plinksps<-plink(parssps, rescale="HB")
summary(plinksps)
reparspsHB<-link.pars(plinksps)
reparspsHB<-reparspsHB$sps[27:46,1:5]
d1=-1*(reparspsHB[,1]*reparspsHB[,2])
d2=-1*(reparspsHB[,1]*reparspsHB[,3])
d3=-1*(reparspsHB[,1]*reparspsHB[,4])
d4=-1*(reparspsHB[,1]*reparspsHB[,5])
reparspsHB2<-cbind(reparspsHB[,1],d1,d2,d3,d4)
colnames(reparspsHB2)<-c("a1","d1","d2","d3","d4")
object_reparspsHB<-generate.mirt_object(reparspsHB2, "graded")
tccspsHB<-plot(object_reparspsHB, type="score")
scoreHBsps<-as.data.table(cbind(tccspsHB$panel.args[[1]][[1]],tccspsHB$panel.args[[1]][[2]]))

## generating separeately calibrated items for sadbank and dsmsad
comb5uni<-mirt(comb5,1, itemtype="graded", technical=list(NCYCLES=1000))
sadpardsmsad<-coef(comb5uni, simplify=TRUE, IRTpars=TRUE)
sadpardsmsad<-as.matrix(sadpardsmsad$items)
##specific number of categories for both item parameters to be linked
catsad<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
catdsmsad<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
##specific the model to be used for both sets of item paramters
pmsad<-as.poly.mod(26,c("grm"),list(c(1:26)))
pmdsmsad<-as.poly.mod(36,c("grm"),list(c(1:36)))
##indicate the number of common items between the two sets of item parameters
commondsmsad<-matrix(c(1:26, 1:26),26,2)
##combine all info into a single object
parsdsmsad<-as.irt.pars(x=list(sadparsirt,sadpardsmsad), common=list(commondsmsad), cat=list(catsad, catdsmsad), poly.mod=list(pmsad, pmdsmsad),grp.names=c("sadbank","dsmsad"))

##stocking lord
plinkdsmsad<-plink(parsdsmsad, rescale="SL")
summary(plinkdsmsad)
repardsmsadSL<-link.pars(plinkdsmsad)
repardsmsadSL<-repardsmsadSL$dsmsad[27:36,1:5]
d1=-1*(repardsmsadSL[,1]*repardsmsadSL[,2])
d2=-1*(repardsmsadSL[,1]*repardsmsadSL[,3])
d3=-1*(repardsmsadSL[,1]*repardsmsadSL[,4])
d4=-1*(repardsmsadSL[,1]*repardsmsadSL[,5])
repardsmsadSL2<-cbind(repardsmsadSL[,1],d1,d2,d3,d4)
colnames(repardsmsadSL2)<-c("a1","d1","d2","d3","d4")
object_repardsmsadSL<-generate.mirt_object(repardsmsadSL2, "graded")
tccdsmsadSL<-plot(object_repardsmsadSL, type="score")
scoreSLdsmsad<-as.data.table(cbind(tccdsmsadSL$panel.args[[1]][[1]],tccdsmsadSL$panel.args[[1]][[2]]))
##mean mean
plinkdsmsad<-plink(parsdsmsad, rescale="MM")
summary(plinkdsmsad)
repardsmsadMM<-link.pars(plinkdsmsad)
repardsmsadMM<-repardsmsadMM$dsmsad[27:36,1:5]
d1=-1*(repardsmsadMM[,1]*repardsmsadMM[,2])
d2=-1*(repardsmsadMM[,1]*repardsmsadMM[,3])
d3=-1*(repardsmsadMM[,1]*repardsmsadMM[,4])
d4=-1*(repardsmsadMM[,1]*repardsmsadMM[,5])
repardsmsadMM2<-cbind(repardsmsadMM[,1],d1,d2,d3,d4)
colnames(repardsmsadMM2)<-c("a1","d1","d2","d3","d4")
object_repardsmsadMM<-generate.mirt_object(repardsmsadMM2, "graded")
tccdsmsadMM<-plot(object_repardsmsadMM, type="score")
scoreMMdsmsad<-as.data.table(cbind(tccdsmsadMM$panel.args[[1]][[1]],tccdsmsadMM$panel.args[[1]][[2]]))
##mean sigma
plinkdsmsad<-plink(parsdsmsad, rescale="MS")
summary(plinkdsmsad)
repardsmsadMS<-link.pars(plinkdsmsad)
repardsmsadMS<-repardsmsadMS$dsmsad[27:36,1:5]
d1=-1*(repardsmsadMS[,1]*repardsmsadMS[,2])
d2=-1*(repardsmsadMS[,1]*repardsmsadMS[,3])
d3=-1*(repardsmsadMS[,1]*repardsmsadMS[,4])
d4=-1*(repardsmsadMS[,1]*repardsmsadMS[,5])
repardsmsadMS2<-cbind(repardsmsadMS[,1],d1,d2,d3,d4)
colnames(repardsmsadMS2)<-c("a1","d1","d2","d3","d4")
object_repardsmsadMS<-generate.mirt_object(repardsmsadMS2, "graded")
tccdsmsadMS<-plot(object_repardsmsadMS, type="score")
scoreMSdsmsad<-as.data.table(cbind(tccdsmsadMS$panel.args[[1]][[1]],tccdsmsadMS$panel.args[[1]][[2]]))
##haebara
plinkdsmsad<-plink(parsdsmsad, rescale="HB")
summary(plinkdsmsad)
repardsmsadHB<-link.pars(plinkdsmsad)
repardsmsadHB<-repardsmsadHB$dsmsad[27:36,1:5]
d1=-1*(repardsmsadHB[,1]*repardsmsadHB[,2])
d2=-1*(repardsmsadHB[,1]*repardsmsadHB[,3])
d3=-1*(repardsmsadHB[,1]*repardsmsadHB[,4])
d4=-1*(repardsmsadHB[,1]*repardsmsadHB[,5])
repardsmsadHB2<-cbind(repardsmsadHB[,1],d1,d2,d3,d4)
colnames(repardsmsadHB2)<-c("a1","d1","d2","d3","d4")
object_repardsmsadHB<-generate.mirt_object(repardsmsadHB2, "graded")
tccdsmsadHB<-plot(object_repardsmsadHB, type="score")
scoreHBdsmsad<-as.data.table(cbind(tccdsmsadHB$panel.args[[1]][[1]],tccdsmsadHB$panel.args[[1]][[2]]))


## generating separeately calibrated items for sadbank and SIAS-6
comb6uni<-mirt(comb6,1, itemtype="graded")
sadparsias6<-coef(comb6uni, simplify=TRUE, IRTpars=TRUE)
sadparsias6<-as.matrix(sadparsias6$items)
##specific number of categories for both item parameters to be linked
catsad<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
catsias6<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
##specific the model to be used for both sets of item paramters
pmsad<-as.poly.mod(26,c("grm"),list(c(1:26)))
pmsias6<-as.poly.mod(32,c("grm"),list(c(1:32)))
##indicate the number of common items between the two sets of item parameters
commonsias6<-matrix(c(1:26, 1:26),26,2)
##combine all info into a single object
parssias6<-as.irt.pars(x=list(sadparsirt,sadparsias6), common=list(commonsias6), cat=list(catsad, catsias6), poly.mod=list(pmsad, pmsias6),grp.names=c("sadbank","sias6"))

##stocking lord
plinksias6<-plink(parssias6, rescale="SL")
summary(plinksias6)
plinksias6SL<-link.pars(plinksias6)
reparsias6SL<-plinksias6SL$sias6[27:32,1:5]
d1=-1*(reparsias6SL[,1]*reparsias6SL[,2])
d2=-1*(reparsias6SL[,1]*reparsias6SL[,3])
d3=-1*(reparsias6SL[,1]*reparsias6SL[,4])
d4=-1*(reparsias6SL[,1]*reparsias6SL[,5])
reparsias6SL2<-cbind(reparsias6SL[,1],d1,d2,d3,d4)
colnames(reparsias6SL2)<-c("a1","d1","d2","d3","d4")
object_reparsias6SL<-generate.mirt_object(reparsias6SL2, "graded")
tccsias6SL<-plot(object_reparsias6SL, type="score")
scoreSLsias6<-as.data.table(cbind(tccsias6SL$panel.args[[1]][[1]],tccsias6SL$panel.args[[1]][[2]]))
##mean/mean
plinksias6<-plink(parssias6, rescale="MM")
summary(plinksias6)
reparsias6MM<-link.pars(plinksias6)
reparsias6MM<-reparsias6MM$sias6[27:32,1:5]
d1=-1*(reparsias6MM[,1]*reparsias6MM[,2])
d2=-1*(reparsias6MM[,1]*reparsias6MM[,3])
d3=-1*(reparsias6MM[,1]*reparsias6MM[,4])
d4=-1*(reparsias6MM[,1]*reparsias6MM[,5])
reparsias6MM2<-cbind(reparsias6MM[,1],d1,d2,d3,d4)
colnames(reparsias6MM2)<-c("a1","d1","d2","d3","d4")
object_reparsias6MM<-generate.mirt_object(reparsias6MM2, "graded")
tccsias6MM<-plot(object_reparsias6MM, type="score")
scoreMMsias6<-as.data.table(cbind(tccsias6MM$panel.args[[1]][[1]],tccsias6MM$panel.args[[1]][[2]]))
##mean/sigma
plinksias6<-plink(parssias6, rescale="MS")
summary(plinksias6)
reparsias6MS<-link.pars(plinksias6)
reparsias6MS<-reparsias6MS$sias6[27:32,1:5]
d1=-1*(reparsias6MS[,1]*reparsias6MS[,2])
d2=-1*(reparsias6MS[,1]*reparsias6MS[,3])
d3=-1*(reparsias6MS[,1]*reparsias6MS[,4])
d4=-1*(reparsias6MS[,1]*reparsias6MS[,5])
reparsias6MS2<-cbind(reparsias6MS[,1],d1,d2,d3,d4)
colnames(reparsias6MS2)<-c("a1","d1","d2","d3","d4")
object_reparsias6MS<-generate.mirt_object(reparsias6MS2, "graded")
tccsias6MS<-plot(object_reparsias6MS, type="score")
scoreMSsias6<-as.data.table(cbind(tccsias6MS$panel.args[[1]][[1]],tccsias6MS$panel.args[[1]][[2]]))
##haebara
plinksias6<-plink(parssias6, rescale="HB")
summary(plinksias6)
reparsias6HB<-link.pars(plinksias6)
reparsias6HB<-reparsias6HB$sias6[27:32,1:5]
d1=-1*(reparsias6HB[,1]*reparsias6HB[,2])
d2=-1*(reparsias6HB[,1]*reparsias6HB[,3])
d3=-1*(reparsias6HB[,1]*reparsias6HB[,4])
d4=-1*(reparsias6HB[,1]*reparsias6HB[,5])
reparsias6HB2<-cbind(reparsias6HB[,1],d1,d2,d3,d4)
colnames(reparsias6HB2)<-c("a1","d1","d2","d3","d4")
object_reparsias6HB<-generate.mirt_object(reparsias6HB2, "graded")
tccsias6HB<-plot(object_reparsias6HB, type="score")
scoreHBsias6<-as.data.table(cbind(tccsias6HB$panel.args[[1]][[1]],tccsias6HB$panel.args[[1]][[2]]))


## generating separeately calibrated items for sadbank and SPS-6
comb7uni<-mirt(comb7,1, itemtype="graded")
sadparsps6<-coef(comb7uni, simplify=TRUE, IRTpars=TRUE)
sadparsps6<-as.matrix(sadparsps6$items)
##specific number of categories for both item parameters to be linked
catsad<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
catsps6<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
##specific the model to be used for both sets of item paramters
pmsad<-as.poly.mod(26,c("grm"),list(c(1:26)))
pmsps6<-as.poly.mod(32,c("grm"),list(c(1:32)))
##indicate the number of common items between the two sets of item parameters
commonsps6<-matrix(c(1:26, 1:26),26,2)
##combine all info into a single object
parssps6<-as.irt.pars(x=list(sadparsirt,sadparsps6), common=list(commonsps6), cat=list(catsad, catsps6), poly.mod=list(pmsad, pmsps6),grp.names=c("sadbank","sps6"))

##stocking lord
plinksps6<-plink(parssps6, rescale="SL")
summary(plinksps6)
plinksps6SL<-link.pars(plinksps6)
reparsps6SL<-plinksps6SL$sps6[27:32,1:5]
d1=-1*(reparsps6SL[,1]*reparsps6SL[,2])
d2=-1*(reparsps6SL[,1]*reparsps6SL[,3])
d3=-1*(reparsps6SL[,1]*reparsps6SL[,4])
d4=-1*(reparsps6SL[,1]*reparsps6SL[,5])
reparsps6SL2<-cbind(reparsps6SL[,1],d1,d2,d3,d4)
colnames(reparsps6SL2)<-c("a1","d1","d2","d3","d4")
object_reparsps6SL<-generate.mirt_object(reparsps6SL2, "graded")
tccsps6SL<-plot(object_reparsps6SL, type="score")
scoreSLsps6<-as.data.table(cbind(tccsps6SL$panel.args[[1]][[1]],tccsps6SL$panel.args[[1]][[2]]))
##mean/mean
plinksps6<-plink(parssps6, rescale="MM")
summary(plinksps6)
reparsps6MM<-link.pars(plinksps6)
reparsps6MM<-reparsps6MM$sps6[27:32,1:5]
d1=-1*(reparsps6MM[,1]*reparsps6MM[,2])
d2=-1*(reparsps6MM[,1]*reparsps6MM[,3])
d3=-1*(reparsps6MM[,1]*reparsps6MM[,4])
d4=-1*(reparsps6MM[,1]*reparsps6MM[,5])
reparsps6MM2<-cbind(reparsps6MM[,1],d1,d2,d3,d4)
colnames(reparsps6MM2)<-c("a1","d1","d2","d3","d4")
object_reparsps6MM<-generate.mirt_object(reparsps6MM2, "graded")
tccsps6MM<-plot(object_reparsps6MM, type="score")
scoreMMsps6<-as.data.table(cbind(tccsps6MM$panel.args[[1]][[1]],tccsps6MM$panel.args[[1]][[2]]))
##mean/sigma
plinksps6<-plink(parssps6, rescale="MS")
summary(plinksps6)
reparsps6MS<-link.pars(plinksps6)
reparsps6MS<-reparsps6MS$sps6[27:32,1:5]
d1=-1*(reparsps6MS[,1]*reparsps6MS[,2])
d2=-1*(reparsps6MS[,1]*reparsps6MS[,3])
d3=-1*(reparsps6MS[,1]*reparsps6MS[,4])
d4=-1*(reparsps6MS[,1]*reparsps6MS[,5])
reparsps6MS2<-cbind(reparsps6MS[,1],d1,d2,d3,d4)
colnames(reparsps6MS2)<-c("a1","d1","d2","d3","d4")
object_reparsps6MS<-generate.mirt_object(reparsps6MS2, "graded")
tccsps6MS<-plot(object_reparsps6MS, type="score")
scoreMSsps6<-as.data.table(cbind(tccsps6MS$panel.args[[1]][[1]],tccsps6MS$panel.args[[1]][[2]]))
##haebara
plinksps6<-plink(parssps6, rescale="HB")
summary(plinksps6)
reparsps6HB<-link.pars(plinksps6)
reparsps6HB<-reparsps6HB$sps6[27:32,1:5]
d1=-1*(reparsps6HB[,1]*reparsps6HB[,2])
d2=-1*(reparsps6HB[,1]*reparsps6HB[,3])
d3=-1*(reparsps6HB[,1]*reparsps6HB[,4])
d4=-1*(reparsps6HB[,1]*reparsps6HB[,5])
reparsps6HB2<-cbind(reparsps6HB[,1],d1,d2,d3,d4)
colnames(reparsps6HB2)<-c("a1","d1","d2","d3","d4")
object_reparsps6HB<-generate.mirt_object(reparsps6HB2, "graded")
tccsps6HB<-plot(object_reparsps6HB, type="score")
scoreHBsps6<-as.data.table(cbind(tccsps6HB$panel.args[[1]][[1]],tccsps6HB$panel.args[[1]][[2]]))


## generating separeately calibrated items for sadbank and GENSAD-12
comb8uni<-mirt(comb8,1, itemtype="graded", technical=list(NCYCLES=1000))
sadpargen12<-coef(comb8uni, simplify=TRUE, IRTpars=TRUE)
sadpargen12<-as.matrix(sadpargen12$items)
##specific number of categories for both item parameters to be linked
catsad<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
catgen12<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
##specific the model to be used for both sets of item paramters
pmsad<-as.poly.mod(26,c("grm"),list(c(1:26)))
pmgen12<-as.poly.mod(38,c("grm"),list(c(1:38)))
##indicate the number of common items between the two sets of item parameters
commongen12<-matrix(c(1:26, 1:26),26,2)
##combine all info into a single object
parsgen12<-as.irt.pars(x=list(sadparsirt,sadpargen12), common=list(commongen12), cat=list(catsad, catgen12), poly.mod=list(pmsad, pmgen12),grp.names=c("sadbank","gen12"))

##stocking lord
plinkgen12<-plink(parsgen12, rescale="SL")
summary(plinkgen12)
plinkgen12SL<-link.pars(plinkgen12)
repargen12SL<-plinkgen12SL$gen12[27:38,1:5]
d1=-1*(repargen12SL[,1]*repargen12SL[,2])
d2=-1*(repargen12SL[,1]*repargen12SL[,3])
d3=-1*(repargen12SL[,1]*repargen12SL[,4])
d4=-1*(repargen12SL[,1]*repargen12SL[,5])
repargen12SL2<-cbind(repargen12SL[,1],d1,d2,d3,d4)
colnames(repargen12SL2)<-c("a1","d1","d2","d3","d4")
object_repargen12SL<-generate.mirt_object(repargen12SL2, "graded")
tccgen12SL<-plot(object_repargen12SL, type="score")
scoreSLgen12<-as.data.table(cbind(tccgen12SL$panel.args[[1]][[1]],tccgen12SL$panel.args[[1]][[2]]))
##mean/mean
plinkgen12<-plink(parsgen12, rescale="MM")
summary(plinkgen12)
repargen12MM<-link.pars(plinkgen12)
repargen12MM<-repargen12MM$gen12[27:38,1:5]
d1=-1*(repargen12MM[,1]*repargen12MM[,2])
d2=-1*(repargen12MM[,1]*repargen12MM[,3])
d3=-1*(repargen12MM[,1]*repargen12MM[,4])
d4=-1*(repargen12MM[,1]*repargen12MM[,5])
repargen12MM2<-cbind(repargen12MM[,1],d1,d2,d3,d4)
colnames(repargen12MM2)<-c("a1","d1","d2","d3","d4")
object_repargen12MM<-generate.mirt_object(repargen12MM2, "graded")
tccgen12MM<-plot(object_repargen12MM, type="score")
scoreMMgen12<-as.data.table(cbind(tccgen12MM$panel.args[[1]][[1]],tccgen12MM$panel.args[[1]][[2]]))
##mean/sigma
plinkgen12<-plink(parsgen12, rescale="MS")
summary(plinkgen12)
repargen12MS<-link.pars(plinkgen12)
repargen12MS<-repargen12MS$gen12[27:38,1:5]
d1=-1*(repargen12MS[,1]*repargen12MS[,2])
d2=-1*(repargen12MS[,1]*repargen12MS[,3])
d3=-1*(repargen12MS[,1]*repargen12MS[,4])
d4=-1*(repargen12MS[,1]*repargen12MS[,5])
repargen12MS2<-cbind(repargen12MS[,1],d1,d2,d3,d4)
colnames(repargen12MS2)<-c("a1","d1","d2","d3","d4")
object_repargen12MS<-generate.mirt_object(repargen12MS2, "graded")
tccgen12MS<-plot(object_repargen12MS, type="score")
scoreMSgen12<-as.data.table(cbind(tccgen12MS$panel.args[[1]][[1]],tccgen12MS$panel.args[[1]][[2]]))
##haebara
plinkgen12<-plink(parsgen12, rescale="HB")
summary(plinkgen12)
repargen12HB<-link.pars(plinkgen12)
repargen12HB<-repargen12HB$gen12[27:38,1:5]
d1=-1*(repargen12HB[,1]*repargen12HB[,2])
d2=-1*(repargen12HB[,1]*repargen12HB[,3])
d3=-1*(repargen12HB[,1]*repargen12HB[,4])
d4=-1*(repargen12HB[,1]*repargen12HB[,5])
repargen12HB2<-cbind(repargen12HB[,1],d1,d2,d3,d4)
colnames(repargen12HB2)<-c("a1","d1","d2","d3","d4")
object_repargen12HB<-generate.mirt_object(repargen12HB2, "graded")
tccgen12HB<-plot(object_repargen12HB, type="score")
scoreHBgen12<-as.data.table(cbind(tccgen12HB$panel.args[[1]][[1]],tccgen12HB$panel.args[[1]][[2]]))


## generating separeately calibrated items for sadbank and mini-spin
comb9uni<-mirt(comb9,1, itemtype="graded", technical=list(NCYCLES=1000))
sadpargen12<-coef(comb9uni, simplify=TRUE, IRTpars=TRUE)
sadpargen12<-as.matrix(sadpargen12$items)
##specific number of categories for both item parameters to be linked
catsad<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
catgen12<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
##specific the model to be used for both sets of item paramters
pmsad<-as.poly.mod(26,c("grm"),list(c(1:26)))
pmgen12<-as.poly.mod(38,c("grm"),list(c(1:38)))
##indicate the number of common items between the two sets of item parameters
commongen12<-matrix(c(1:26, 1:26),26,2)
##combine all info into a single object
parsgen12<-as.irt.pars(x=list(sadparsirt,sadpargen12), common=list(commongen12), cat=list(catsad, catgen12), poly.mod=list(pmsad, pmgen12),grp.names=c("sadbank","gen12"))

##stocking lord
plinkgen12<-plink(parsgen12, rescale="SL")
summary(plinkgen12)
plinkgen12SL<-link.pars(plinkgen12)
repargen12SL<-plinkgen12SL$gen12[27:38,1:5]
d1=-1*(repargen12SL[,1]*repargen12SL[,2])
d2=-1*(repargen12SL[,1]*repargen12SL[,3])
d3=-1*(repargen12SL[,1]*repargen12SL[,4])
d4=-1*(repargen12SL[,1]*repargen12SL[,5])
repargen12SL2<-cbind(repargen12SL[,1],d1,d2,d3,d4)
colnames(repargen12SL2)<-c("a1","d1","d2","d3","d4")
object_repargen12SL<-generate.mirt_object(repargen12SL2, "graded")
tccgen12SL<-plot(object_repargen12SL, type="score")
scoreSLgen12<-as.data.table(cbind(tccgen12SL$panel.args[[1]][[1]],tccgen12SL$panel.args[[1]][[2]]))
##mean/mean
plinkgen12<-plink(parsgen12, rescale="MM")
summary(plinkgen12)
repargen12MM<-link.pars(plinkgen12)
repargen12MM<-repargen12MM$gen12[27:38,1:5]
d1=-1*(repargen12MM[,1]*repargen12MM[,2])
d2=-1*(repargen12MM[,1]*repargen12MM[,3])
d3=-1*(repargen12MM[,1]*repargen12MM[,4])
d4=-1*(repargen12MM[,1]*repargen12MM[,5])
repargen12MM2<-cbind(repargen12MM[,1],d1,d2,d3,d4)
colnames(repargen12MM2)<-c("a1","d1","d2","d3","d4")
object_repargen12MM<-generate.mirt_object(repargen12MM2, "graded")
tccgen12MM<-plot(object_repargen12MM, type="score")
scoreMMgen12<-as.data.table(cbind(tccgen12MM$panel.args[[1]][[1]],tccgen12MM$panel.args[[1]][[2]]))
##mean/sigma
plinkgen12<-plink(parsgen12, rescale="MS")
summary(plinkgen12)
repargen12MS<-link.pars(plinkgen12)
repargen12MS<-repargen12MS$gen12[27:38,1:5]
d1=-1*(repargen12MS[,1]*repargen12MS[,2])
d2=-1*(repargen12MS[,1]*repargen12MS[,3])
d3=-1*(repargen12MS[,1]*repargen12MS[,4])
d4=-1*(repargen12MS[,1]*repargen12MS[,5])
repargen12MS2<-cbind(repargen12MS[,1],d1,d2,d3,d4)
colnames(repargen12MS2)<-c("a1","d1","d2","d3","d4")
object_repargen12MS<-generate.mirt_object(repargen12MS2, "graded")
tccgen12MS<-plot(object_repargen12MS, type="score")
scoreMSgen12<-as.data.table(cbind(tccgen12MS$panel.args[[1]][[1]],tccgen12MS$panel.args[[1]][[2]]))
##haebara
plinkgen12<-plink(parsgen12, rescale="HB")
summary(plinkgen12)
repargen12HB<-link.pars(plinkgen12)
repargen12HB<-repargen12HB$gen12[27:38,1:5]
d1=-1*(repargen12HB[,1]*repargen12HB[,2])
d2=-1*(repargen12HB[,1]*repargen12HB[,3])
d3=-1*(repargen12HB[,1]*repargen12HB[,4])
d4=-1*(repargen12HB[,1]*repargen12HB[,5])
repargen12HB2<-cbind(repargen12HB[,1],d1,d2,d3,d4)
colnames(repargen12HB2)<-c("a1","d1","d2","d3","d4")
object_repargen12HB<-generate.mirt_object(repargen12HB2, "graded")
tccgen12HB<-plot(object_repargen12HB, type="score")
scoreHBgen12<-as.data.table(cbind(tccgen12HB$panel.args[[1]][[1]],tccgen12HB$panel.args[[1]][[2]]))


### agreement between factor scores generated by IRT methods
eapcomb<-cbind(socbeap, lsaseap, spineap, siaseap, spseap, dsmeap)
colnames(eapcomb)<-c("Socbank","Socbank SE", "LSAS","LSAS SE","SPIN","SPIN SE","SIAS", "SIAS SE", "SPS","SPS SE", "DSM SAD","DSM SAD SE")

mean(eapcomb[,c("Socbank")])
SD(eapcomb[,c("Socbank")])
mean(eapcomb[,c("LSAS")])
SD(eapcomb[,c("LSAS")])
mean(eapcomb[,c("SPIN")])
SD(eapcomb[,c("SPIN")])
mean(eapcomb[,c("SIAS")])
SD(eapcomb[,c("SIAS")])
mean(eapcomb[,c("SPS")])
SD(eapcomb[,c("SPS")])
mean(eapcomb[,c("DSM SAD")])
SD(eapcomb[,c("DSM SAD")])
mean(sias6eap[ ,1])
SD(sias6eap[ ,1])
mean(sps6eap[ ,1])
SD(sps6eap[ ,1])
mean(sps12eap[ ,1])
SD(sps12eap[ ,1])
mean(mspineap[ ,1])
SD(mspineap[ ,1])


cor(eapcomb[,c("Socbank")], eapcomb[,c("LSAS")])
cor(eapcomb[,c("Socbank")], eapcomb[,c("SPIN")])
cor(eapcomb[,c("Socbank")], eapcomb[,c("SIAS")])
cor(eapcomb[,c("Socbank")], eapcomb[,c("SPS")])
cor(eapcomb[,c("Socbank")], eapcomb[,c("DSM SAD")])
cor(socbeap[,1], sias6eap[,1])
cor(socbeap[,1], sps6eap[,1])
cor(socbeap[,1], sps12eap[,1])
cor(socbeap[,1], mspineap[,1])

cor(eapcomb[,c("Socbank")], spsmmfscores[,c("F1")])

mean(eapcomb[,c("Socbank")]-eapcomb[,c("LSAS")])
mean(eapcomb[,c("Socbank")]-eapcomb[,c("SPIN")])
mean(eapcomb[,c("Socbank")]-eapcomb[,c("SIAS")])
mean(eapcomb[,c("Socbank")]-eapcomb[,c("SPS")])
mean(eapcomb[,c("Socbank")]-eapcomb[,c("DSM SAD")])
mean(socbeap[,1]-sias6eap[,1])
mean(socbeap[,1]-sps6eap[,1])
mean(socbeap[,1]-sps12eap[,1])
mean(socbeap[,1]-mspineap[,1])

mean(eapcomb[,c("Socbank")]-spsmmfscores[,c("F1")])

SD(eapcomb[,c("Socbank")]-eapcomb[,c("LSAS")])
SD(eapcomb[,c("Socbank")]-eapcomb[,c("SPIN")])
SD(eapcomb[,c("Socbank")]-eapcomb[,c("SIAS")])
SD(eapcomb[,c("Socbank")]-eapcomb[,c("SPS")])
SD(eapcomb[,c("Socbank")]-eapcomb[,c("DSM SAD")])
SD(socbeap[,1]-sias6eap[,1])
SD(socbeap[,1]-sps6eap[,1])
SD(socbeap[,1]-sps12eap[,1])
SD(socbeap[,1]-mspineap[,1])

MSD <- mean(((eapcomb[,c("Socbank")] - eapcomb[,c("LSAS")])^2))
(sqrt(MSD))
MSD <- mean(((eapcomb[,c("Socbank")] - eapcomb[,c("SPIN")])^2))
(sqrt(MSD))
MSD <- mean(((eapcomb[,c("Socbank")] - eapcomb[,c("SIAS")])^2))
(sqrt(MSD))
MSD <- mean(((eapcomb[,c("Socbank")] - eapcomb[,c("SPS")])^2))
(sqrt(MSD))
MSD <- mean(((eapcomb[,c("Socbank")] - eapcomb[,c("DSM SAD")])^2))
(sqrt(MSD))

MSD <- mean(((eapcomb[,c("Socbank")] - spsmmfscores[,c("F1")])^2))
(sqrt(MSD))

ICC(as.matrix(cbind(eapcomb[,c("Socbank")],eapcomb[,c("LSAS")])))
ICC(as.matrix(cbind(eapcomb[,c("Socbank")],eapcomb[,c("SPIN")])))
ICC(as.matrix(cbind(eapcomb[,c("Socbank")],eapcomb[,c("SIAS")])))
ICC(as.matrix(cbind(eapcomb[,c("Socbank")],eapcomb[,c("SPS")])))
ICC(as.matrix(cbind(eapcomb[,c("Socbank")],eapcomb[,c("DSM SAD")])))
ICC(as.matrix(cbind(socbeap[,1],sias6eap[,1])))
ICC(as.matrix(cbind(socbeap[,1],sps6eap[,1])))
ICC(as.matrix(cbind(socbeap[,1],sps12eap[,1])))
ICC(as.matrix(cbind(socbeap[,1],mspineap[,1])))

bland.altman.plot(eapcomb[,c("Socbank")],eapcomb[,c("LSAS")], main="SAD-bank and LSAS", xlab="Mean", ylab="Differences")
bland.altman.plot(eapcomb[,c("Socbank")],eapcomb[,c("SPIN")], main="SAD-bank and SPIN", xlab="Mean", ylab="Differences")
bland.altman.plot(eapcomb[,c("Socbank")],eapcomb[,c("SIAS")], main="SAD-bank and SIAS-20", xlab="Mean", ylab="Differences")
bland.altman.plot(eapcomb[,c("Socbank")],eapcomb[,c("SPS")], main="SAD-bank and SPS-20", xlab="Mean", ylab="Differences")
bland.altman.plot(eapcomb[,c("Socbank")],eapcomb[,c("DSM SAD")], main="SAD-bank and DSM SAD", xlab="Mean", ylab="Differences")
bland.altman.plot(socbeap[,1],sias6eap[,1], main="SAD-bank and SIAS-6", xlab="Mean", ylab="Differences")
bland.altman.plot(socbeap[,1],sps6eap[,1], main="SAD-bank and SPS-6", xlab="Mean", ylab="Differences")
bland.altman.plot(socbeap[,1],sps12eap[,1], main="SAD-bank and GSAS-12", xlab="Mean", ylab="Differences")
bland.altman.plot(socbeap[,1],mspineap[,1], main="SAD-bank and Mini-SPIN", xlab="Mean", ylab="Differences")


##cross walked scores from fixed parameter linking
##lsas
lsastot<-subset(comdata, select=c(LSAStot))
rawlsas<-lsaseapsum$Sum.Scores
lsastheta<-lsaseapsum$Theta
recodemat<-data.frame(rawlsas,lsastheta)

idx<-match(lsastot$LSAStot, recodemat$rawlsas)
lsas_cross<-as.matrix(cbind(lsastot, recodemat[idx,-1]))

mean(lsas_cross$`recodemat[idx, -1]`)
SD(lsas_cross$`recodemat[idx, -1]`)

cor(socbeap[,1],lsas_cross[,2])
ICC(as.matrix(cbind(socbeap[,1],lsas_cross[,2])))
mean(socbeap[,1]-lsas_cross[,2])
SD(socbeap[,1]-lsas_cross[,2])

#spin
spintot<-subset(comdata, select=c(SPINtot))
rawspin<-spineapsum$Sum.Scores
spintheta<-spineapsum$Theta
recodemat<-data.frame(rawspin,spintheta)

idx<-match(spintot$SPINtot, recodemat$rawspin)
spin_cross<-as.matrix(cbind(spintot, recodemat[idx,-1]))

mean(spin_cross[,2])
SD(spin_cross[,2])

cor(socbeap[,1],spin_cross[,2])
ICC(as.matrix(cbind(socbeap[,1],spin_cross[,2])))
mean(socbeap[,1]-spin_cross[,2])
SD(socbeap[,1]-spin_cross[,2])

#sias
siastot<-subset(comdata, select=c(SIAStot))
rawsias<-siaseapsum$Sum.Scores
siastheta<-siaseapsum$Theta
recodemat<-data.frame(rawsias,siastheta)

idx<-match(siastot$SIAStot, recodemat$rawsias)
sias_cross<-as.matrix(cbind(siastot, recodemat[idx,-1]))

mean(sias_cross[,2])
SD(sias_cross[,2])

cor(socbeap[,1],sias_cross[,2])
ICC(as.matrix(cbind(socbeap[,1],sias_cross[,2])))
mean(socbeap[,1]-sias_cross[,2])
SD(socbeap[,1]-sias_cross[,2])

#sps
spstot<-subset(comdata, select=c(SPStot))
rawsps<-spseapsum$Sum.Scores
spstheta<-spseapsum$Theta
recodemat<-data.frame(rawsps,spstheta)

idx<-match(spstot$SPStot, recodemat$rawsps)
sps_cross<-as.matrix(cbind(spstot, recodemat[idx,-1]))

mean(sps_cross[,2])
SD(sps_cross[,2])

cor(socbeap[,1],sps_cross[,2])
ICC(as.matrix(cbind(socbeap[,1],sps_cross[,2])))
mean(socbeap[,1]-sps_cross[,2])
SD(socbeap[,1]-sps_cross[,2])

#dsmsad
dsmtot<-subset(comdata, select=c(DSMSADtot))
rawdsm<-dsmeapsum$Sum.Scores
dsmtheta<-dsmeapsum$Theta
recodemat<-data.frame(rawdsm,dsmtheta)

idx<-match(dsmtot$DSMSADtot, recodemat$rawdsm)
dsm_cross<-as.matrix(cbind(dsmtot, recodemat[idx,-1]))

mean(dsm_cross[,2])
SD(dsm_cross[,2])

cor(socbeap[,1],dsm_cross[,2])
ICC(as.matrix(cbind(socbeap[,1],dsm_cross[,2])))
mean(socbeap[,1]-dsm_cross[,2])
SD(socbeap[,1]-dsm_cross[,2])

#minispin
mspintot<-subset(comdata, select=c(SPIN3tot))
rawmspin<-mspineapsum$Sum.Scores
mspintheta<-mspineapsum$Theta
recodemat<-data.frame(rawmspin,mspintheta)

idx<-match(mspintot$SPIN3tot, recodemat$rawmspin)
mspin_cross<-as.matrix(cbind(mspintot, recodemat[idx,-1]))

mean(mspin_cross[,2])
SD(mspin_cross[,2])

cor(socbeap[,1],mspin_cross[,2])
ICC(as.matrix(cbind(socbeap[,1],mspin_cross[,2])))
mean(socbeap[,1]-mspin_cross[,2])
SD(socbeap[,1]-mspin_cross[,2])

#sias6
sias6tot<-subset(comdata, select=c(SIAS6tot))
rawsias6<-sias6eapsum$Sum.Scores
sias6theta<-sias6eapsum$Theta
recodemat<-data.frame(rawsias6,sias6theta)

idx<-match(sias6tot$SIAS6tot, recodemat$rawsias6)
sias6_cross<-as.matrix(cbind(sias6tot, recodemat[idx,-1]))

mean(sias6_cross[,2])
SD(sias6_cross[,2])

cor(socbeap[,1],sias6_cross[,2])
ICC(as.matrix(cbind(socbeap[,1],sias6_cross[,2])))
mean(socbeap[,1]-sias6_cross[,2])
SD(socbeap[,1]-sias6_cross[,2])

#sps6
sps6tot<-subset(comdata, select=c(SPS6tot))
rawsps6<-sps6eapsum$Sum.Scores
sps6theta<-sps6eapsum$Theta
recodemat<-data.frame(rawsps6,sps6theta)

idx<-match(sps6tot$SPS6tot, recodemat$rawsps6)
sps6_cross<-as.matrix(cbind(sps6tot, recodemat[idx,-1]))

mean(sps6_cross[,2])
SD(sps6_cross[,2])

cor(socbeap[,1],sps6_cross[,2])
ICC(as.matrix(cbind(socbeap[,1],sps6_cross[,2])))
mean(socbeap[,1]-sps6_cross[,2])
SD(socbeap[,1]-sps6_cross[,2])


#gensad12
gensad12tot<-subset(comdata, select=c(GENSAD12tot))
rawgensad12<-sps12eapsum$Sum.Scores
gensad12theta<-sps12eapsum$Theta
recodemat<-data.frame(rawgensad12,gensad12theta)

idx<-match(gensad12tot$GENSAD12tot, recodemat$rawgensad12)
gensad12_cross<-as.matrix(cbind(gensad12tot, recodemat[idx,-1]))

mean(gensad12_cross[,2])
SD(gensad12_cross[,2])

cor(socbeap[,1],gensad12_cross[,2])
ICC(as.matrix(cbind(socbeap[,1],gensad12_cross[,2])))
mean(socbeap[,1]-gensad12_cross[,2])
SD(socbeap[,1]-gensad12_cross[,2])

##population invariance testing
##lsas gender
lsasmale <- subset(comdata, D2==1, select=c(LSASr1:LSASr24))
lsasfemale <- subset(comdata, D2==2, select=c(LSASr1:LSASr24))
comb1male <- subset(comdata, D2==1, select=c(SADBNK1:SADBNK26,LSASr1:LSASr24))
comb1female <- subset(comdata, D2==2, select=c(SADBNK1:SADBNK26,LSASr1:LSASr24))

##males
##SAD item bank and LSAS
## generating LSAS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb1male,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parlsas<-mirt(comb1male,1, itemtype="graded", pars=sv2)
coef(parlsas, IRTpars=TRUE, simplify=TRUE)
coeflsas<-coef(parlsas, as.data.frame=TRUE)
eap<-fscores(parlsas, method="EAP")
sum(eap[1:208,])/208

##calculate lsas cross-walk tables
svlsas <-mirt(lsasmale,1, itemtype="graded", pars="values")
sv2lsas<-svlsas
sv2lsas$value<-c(coeflsas[131:228,])

sv2lsas$est<- FALSE

lsaswalk<-mirt(lsasmale,1, itemtype="graded", pars=sv2lsas)
lsaseap<-fscores(lsaswalk, method="EAP", full.scores.SE=TRUE)
sum(lsaseap[1:208,1])/208
lsaseapsum<-fscores(lsaswalk, method="EAPsum", full.scores = FALSE)

##females
##SAD item bank and LSAS
## generating LSAS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb1female,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parlsas<-mirt(comb1female,1, itemtype="graded", pars=sv2)
coef(parlsas, IRTpars=TRUE, simplify=TRUE)
coeflsas<-coef(parlsas, as.data.frame=TRUE)
eap<-fscores(parlsas, method="EAP")
sum(eap[1:813,])/813

##calculate lsas cross-walk tables
svlsas <-mirt(lsasfemale,1, itemtype="graded", pars="values")
sv2lsas<-svlsas
sv2lsas$value<-c(coeflsas[131:228,])

sv2lsas$est<- FALSE

lsaswalk<-mirt(lsasfemale,1, itemtype="graded", pars=sv2lsas)
lsaseap<-fscores(lsaswalk, method="EAP", full.scores.SE=TRUE)
SD(lsaseap[,1])
sum(lsaseap[1:813,1])/813
lsaseapsum<-fscores(lsaswalk, method="EAPsum", full.scores = FALSE)

##spin gender
spinmale <- subset(comdata, D2==1, select=c(SPIN1:SPIN17))
spinfemale <- subset(comdata, D2==2, select=c(SPIN1:SPIN17))
comb2male <- subset(comdata, D2==1, select=c(SADBNK1:SADBNK26,SPIN1:SPIN17))
comb2female <- subset(comdata, D2==2, select=c(SADBNK1:SADBNK26,SPIN1:SPIN17))

##males
##SAD item bank and SPIN
## generating SPIN irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb2male,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parspin<-mirt(comb2male,1, itemtype="graded", pars=sv2)
coef(parspin, IRTpars=TRUE, simplify=TRUE)
coefspin<-coef(parspin, as.data.frame=TRUE)
eap<-fscores(parspin, method="EAP")
sum(eap[1:208,])/208

##calculate spin cross-walk tables
svspin <-mirt(spinmale,1, itemtype="graded", pars="values")
sv2spin<-svspin
sv2spin$value<-c(coefspin[131:217,])

sv2spin$est<- FALSE

spinwalk<-mirt(spinmale,1, itemtype="graded", pars=sv2spin)
spineap<-fscores(spinwalk, method="EAP", full.scores.SE=TRUE)
sum(spineap[1:208,1])/208
spineapsum<-fscores(spinwalk, method="EAPsum", full.scores = FALSE)

##female
##SAD item bank and SPIN
## generating SPIN irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb2female,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parspin<-mirt(comb2female,1, itemtype="graded", pars=sv2)
coef(parspin, IRTpars=TRUE, simplify=TRUE)
coefspin<-coef(parspin, as.data.frame=TRUE)
eap<-fscores(parspin, method="EAP")
sum(eap[1:208,])/208

##calculate spin cross-walk tables
svspin <-mirt(spinfemale,1, itemtype="graded", pars="values")
sv2spin<-svspin
sv2spin$value<-c(coefspin[131:217,])

sv2spin$est<- FALSE

spinwalk<-mirt(spinfemale,1, itemtype="graded", pars=sv2spin)
spineap<-fscores(spinwalk, method="EAP", full.scores.SE=TRUE)
sum(spineap[1:208,1])/208
spineapsum<-fscores(spinwalk, method="EAPsum", full.scores = FALSE)

#sias gender 

siasmale <- subset(comdata, D2==1, select=c(SIAS1:SIAS20))
siasfemale <- subset(comdata, D2==2, select=c(SIAS1:SIAS20))
comb3male <- subset(comdata, D2==1, select=c(SADBNK1:SADBNK26,SIAS1:SIAS20))
comb3female <- subset(comdata, D2==2, select=c(SADBNK1:SADBNK26,SIAS1:SIAS20))

#male
##SAD item bank and SIAS
## generating SIAS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb3male,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parsias<-mirt(comb3male,1, itemtype="graded", pars=sv2)
coef(parsias, IRTpars=TRUE, simplify=TRUE)
coefsias<-coef(parsias, as.data.frame=TRUE)
eap<-fscores(parsias, method="EAP")
sum(eap[1:208,])/208

##calculate sias cross-walk tables
svsias <-mirt(siasmale,1, itemtype="graded", pars="values")
sv2sias<-svsias
sv2sias$value<-c(coefsias[131:232,])

sv2sias$est<- FALSE

siaswalk<-mirt(siasmale,1, itemtype="graded", pars=sv2sias)
coef(siaswalk, simplify=TRUE)

siaseap<-fscores(siaswalk, method="EAP", full.scores.SE=TRUE)
siaseapsum<-fscores(siaswalk, method="EAPsum", full.scores = FALSE)

#female
##SAD item bank and SIAS
## generating SIAS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb3female,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parsias<-mirt(comb3female,1, itemtype="graded", pars=sv2)
coef(parsias, IRTpars=TRUE, simplify=TRUE)
coefsias<-coef(parsias, as.data.frame=TRUE)
eap<-fscores(parsias, method="EAP")
sum(eap[1:813,])/813

##calculate sias cross-walk tables
svsias <-mirt(siasfemale,1, itemtype="graded", pars="values")
sv2sias<-svsias
sv2sias$value<-c(coefsias[131:232,])

sv2sias$est<- FALSE

siaswalk<-mirt(siasfemale,1, itemtype="graded", pars=sv2sias)
coef(siaswalk, simplify=TRUE)

siaseap<-fscores(siaswalk, method="EAP", full.scores.SE=TRUE)
siaseapsum<-fscores(siaswalk, method="EAPsum", full.scores = FALSE)

##sps gender
spsmale <- subset(comdata, D2==1, select=c(SPS1:SPS20))
spsfemale <- subset(comdata, D2==2, select=c(SPS1:SPS20))
comb4male <- subset(comdata, D2==1, select=c(SADBNK1:SADBNK26,SPS1:SPS20))
comb4female <- subset(comdata, D2==2, select=c(SADBNK1:SADBNK26,SPS1:SPS20))

#male
##SAD item bank and SPS
## generating SPS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb4male,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parsps<-mirt(comb4male,1, itemtype="graded", pars=sv2)
coef(parsps, IRTpars=TRUE, simplify=TRUE)
coefsps<-coef(parsps, as.data.frame=TRUE)
eap<-fscores(parsps, method="EAP")

##calculate sps cross-walk tables
svsps <-mirt(spsmale,1, itemtype="graded", pars="values")
sv2sps<-svsps
sv2sps$value<-c(coefsps[131:232,])

sv2sps$est<- FALSE

spswalk<-mirt(spsmale,1, itemtype="graded", pars=sv2sps)
coef(spswalk, simplify=TRUE)

spseap<-fscores(spswalk, method="EAP", full.scores.SE=TRUE)
spseapsum<-fscores(spswalk, method="EAPsum", full.scores = FALSE)

#female
##SAD item bank and SPS
## generating SPS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb4female,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parsps<-mirt(comb4female,1, itemtype="graded", pars=sv2)
coef(parsps, IRTpars=TRUE, simplify=TRUE)
coefsps<-coef(parsps, as.data.frame=TRUE)
eap<-fscores(parsps, method="EAP")

##calculate sps cross-walk tables
svsps <-mirt(spsfemale,1, itemtype="graded", pars="values")
sv2sps<-svsps
sv2sps$value<-c(coefsps[131:232,])

sv2sps$est<- FALSE

spswalk<-mirt(spsfemale,1, itemtype="graded", pars=sv2sps)
coef(spswalk, simplify=TRUE)

spseap<-fscores(spswalk, method="EAP", full.scores.SE=TRUE)
spseapsum<-fscores(spswalk, method="EAPsum", full.scores = FALSE)


#dsmsad and gender
dsmsadmale <- subset(comdata, D2==1, select=c(DSMSAD1:DSMSAD10))
dsmsadfemale <- subset(comdata, D2==2, select=c(DSMSAD1:DSMSAD10))
comb5male <- subset(comdata, D2==1, select=c(SADBNK1:SADBNK26,DSMSAD1:DSMSAD10))
comb5female <- subset(comdata, D2==2, select=c(SADBNK1:SADBNK26,DSMSAD1:DSMSAD10))

#male
##SAD item bank and dsmsad
## generating dsmsad irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb5male,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

pardsm<-mirt(comb5male,1, itemtype="graded", pars=sv2)
coef(pardsm, IRTpars=TRUE, simplify=TRUE)
coefdsm<-coef(pardsm, as.data.frame=TRUE)

##calculate dsmsad cross-walk tables
svdsm <-mirt(dsmsadmale,1, itemtype="graded", pars="values")
sv2dsm<-svdsm
sv2dsm$value<-c(coefdsm[131:182,])

sv2dsm$est<- FALSE

dsmwalk<-mirt(dsmsadmale,1, itemtype="graded", pars=sv2dsm)
coef(dsmwalk, simplify=TRUE)

dsmeap<-fscores(dsmwalk, method="EAP", full.scores.SE=TRUE)
dsmeapsum<-fscores(dsmwalk, method="EAPsum", full.scores = FALSE)

#female
##SAD item bank and dsmsad
## generating dsmsad irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb5female,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

pardsm<-mirt(comb5female,1, itemtype="graded", pars=sv2)
coef(pardsm, IRTpars=TRUE, simplify=TRUE)
coefdsm<-coef(pardsm, as.data.frame=TRUE)

##calculate dsmsad cross-walk tables
svdsm <-mirt(dsmsadfemale,1, itemtype="graded", pars="values")
sv2dsm<-svdsm
sv2dsm$value<-c(coefdsm[131:182,])

sv2dsm$est<- FALSE

dsmwalk<-mirt(dsmsadfemale,1, itemtype="graded", pars=sv2dsm)
coef(dsmwalk, simplify=TRUE)

dsmeap<-fscores(dsmwalk, method="EAP", full.scores.SE=TRUE)
dsmeapsum<-fscores(dsmwalk, method="EAPsum", full.scores = FALSE)



## age invariance 
##population invariance testing
##lsas age
lsasyoung <- subset(comdata, D1>=2 & D1<=4, select=c(LSASr1:LSASr24))
lsasold <- subset(comdata, D1>=5 & D1<=7, select=c(LSASr1:LSASr24))
comb1young <- subset(comdata, D1>=2 & D1<=4, select=c(SADBNK1:SADBNK26,LSASr1:LSASr24))
comb1old <- subset(comdata, D1>=5 & D1<=7, select=c(SADBNK1:SADBNK26,LSASr1:LSASr24))

##young
##SAD item bank and LSAS
## generating LSAS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb1young,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parlsas<-mirt(comb1young,1, itemtype="graded", pars=sv2)
coef(parlsas, IRTpars=TRUE, simplify=TRUE)
coeflsas<-coef(parlsas, as.data.frame=TRUE)


##calculate lsas cross-walk tables
svlsas <-mirt(lsasyoung,1, itemtype="graded", pars="values")
sv2lsas<-svlsas
sv2lsas$value<-c(coeflsas[131:228,])

sv2lsas$est<- FALSE

lsaswalk<-mirt(lsasyoung,1, itemtype="graded", pars=sv2lsas)
lsaseapsum<-fscores(lsaswalk, method="EAPsum", full.scores = FALSE)

##old
##SAD item bank and LSAS
## generating LSAS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb1old,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parlsas<-mirt(comb1old,1, itemtype="graded", pars=sv2)
coef(parlsas, IRTpars=TRUE, simplify=TRUE)
coeflsas<-coef(parlsas, as.data.frame=TRUE)


##calculate lsas cross-walk tables
svlsas <-mirt(lsasold,1, itemtype="graded", pars="values")
sv2lsas<-svlsas
sv2lsas$value<-c(coeflsas[131:228,])

sv2lsas$est<- FALSE

lsaswalk<-mirt(lsasold,1, itemtype="graded", pars=sv2lsas)
lsaseap<-fscores(lsaswalk, method="EAP", full.scores.SE=TRUE)
lsaseapsum<-fscores(lsaswalk, method="EAPsum", full.scores = FALSE)

##spin age
spinyoung <- subset(comdata, D1>=2 & D1<=4, select=c(SPIN1:SPIN17))
spinold <- subset(comdata, D1>=5 & D1<=7, select=c(SPIN1:SPIN17))
comb2young <- subset(comdata, D1>=2 & D1<=4, select=c(SADBNK1:SADBNK26,SPIN1:SPIN17))
comb2old <- subset(comdata, D1>=5 & D1<=7, select=c(SADBNK1:SADBNK26,SPIN1:SPIN17))

##young
##SAD item bank and SPIN
## generating SPIN irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb2young,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parspin<-mirt(comb2young,1, itemtype="graded", pars=sv2)
coef(parspin, IRTpars=TRUE, simplify=TRUE)
coefspin<-coef(parspin, as.data.frame=TRUE)

##calculate spin cross-walk tables
svspin <-mirt(spinyoung,1, itemtype="graded", pars="values")
sv2spin<-svspin
sv2spin$value<-c(coefspin[131:217,])

sv2spin$est<- FALSE

spinwalk<-mirt(spinyoung,1, itemtype="graded", pars=sv2spin)
spineapsum<-fscores(spinwalk, method="EAPsum", full.scores = FALSE)

##old
##SAD item bank and SPIN
## generating SPIN irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb2old,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parspin<-mirt(comb2old,1, itemtype="graded", pars=sv2)
coef(parspin, IRTpars=TRUE, simplify=TRUE)
coefspin<-coef(parspin, as.data.frame=TRUE)

##calculate spin cross-walk tables
svspin <-mirt(spinold,1, itemtype="graded", pars="values")
sv2spin<-svspin
sv2spin$value<-c(coefspin[131:217,])

sv2spin$est<- FALSE

spinwalk<-mirt(spinold,1, itemtype="graded", pars=sv2spin)
spineapsum<-fscores(spinwalk, method="EAPsum", full.scores = FALSE)

#sias age 

siasyoung <- subset(comdata, D1>=2 & D1<=4, select=c(SIAS1:SIAS20))
siasold <- subset(comdata, D1>=5 & D1<=7, select=c(SIAS1:SIAS20))
comb3young <- subset(comdata, D1>=2 & D1<=4, select=c(SADBNK1:SADBNK26,SIAS1:SIAS20))
comb3old <- subset(comdata, D1>=5 & D1<=7, select=c(SADBNK1:SADBNK26,SIAS1:SIAS20))

#young
##SAD item bank and SIAS
## generating SIAS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb3young,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parsias<-mirt(comb3young,1, itemtype="graded", pars=sv2)
coef(parsias, IRTpars=TRUE, simplify=TRUE)
coefsias<-coef(parsias, as.data.frame=TRUE)

##calculate sias cross-walk tables
svsias <-mirt(siasyoung,1, itemtype="graded", pars="values")
sv2sias<-svsias
sv2sias$value<-c(coefsias[131:232,])

sv2sias$est<- FALSE

siaswalk<-mirt(siasyoung,1, itemtype="graded", pars=sv2sias)
coef(siaswalk, simplify=TRUE)

siaseap<-fscores(siaswalk, method="EAP", full.scores.SE=TRUE)
siaseapsum<-fscores(siaswalk, method="EAPsum", full.scores = FALSE)

#old
##SAD item bank and SIAS
## generating SIAS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb3old,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parsias<-mirt(comb3old,1, itemtype="graded", pars=sv2)
coef(parsias, IRTpars=TRUE, simplify=TRUE)
coefsias<-coef(parsias, as.data.frame=TRUE)


##calculate sias cross-walk tables
svsias <-mirt(siasold,1, itemtype="graded", pars="values")
sv2sias<-svsias
sv2sias$value<-c(coefsias[131:232,])

sv2sias$est<- FALSE

siaswalk<-mirt(siasold,1, itemtype="graded", pars=sv2sias)
coef(siaswalk, simplify=TRUE)

siaseap<-fscores(siaswalk, method="EAP", full.scores.SE=TRUE)
siaseapsum<-fscores(siaswalk, method="EAPsum", full.scores = FALSE)

##sps age
spsyoung <- subset(comdata, D1>=2 & D1<=4, select=c(SPS1:SPS20))
spsold <- subset(comdata, D1>=5 & D1<=7, select=c(SPS1:SPS20))
comb4young <- subset(comdata, D1>=2 & D1<=4, select=c(SADBNK1:SADBNK26,SPS1:SPS20))
comb4old <- subset(comdata, D1>=5 & D1<=7, select=c(SADBNK1:SADBNK26,SPS1:SPS20))

#young
##SAD item bank and SPS
## generating SPS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb4young,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parsps<-mirt(comb4young,1, itemtype="graded", pars=sv2)
coef(parsps, IRTpars=TRUE, simplify=TRUE)
coefsps<-coef(parsps, as.data.frame=TRUE)


##calculate sps cross-walk tables
svsps <-mirt(spsyoung,1, itemtype="graded", pars="values")
sv2sps<-svsps
sv2sps$value<-c(coefsps[131:232,])

sv2sps$est<- FALSE

spswalk<-mirt(spsyoung,1, itemtype="graded", pars=sv2sps)
coef(spswalk, simplify=TRUE)

spseap<-fscores(spswalk, method="EAP", full.scores.SE=TRUE)
spseapsum<-fscores(spswalk, method="EAPsum", full.scores = FALSE)

#old
##SAD item bank and SPS
## generating SPS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb4old,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parsps<-mirt(comb4old,1, itemtype="graded", pars=sv2)
coef(parsps, IRTpars=TRUE, simplify=TRUE)
coefsps<-coef(parsps, as.data.frame=TRUE)
eap<-fscores(parsps, method="EAP")

##calculate sps cross-walk tables
svsps <-mirt(spsold,1, itemtype="graded", pars="values")
sv2sps<-svsps
sv2sps$value<-c(coefsps[131:232,])

sv2sps$est<- FALSE

spswalk<-mirt(spsold,1, itemtype="graded", pars=sv2sps)
coef(spswalk, simplify=TRUE)

spseap<-fscores(spswalk, method="EAP", full.scores.SE=TRUE)
spseapsum<-fscores(spswalk, method="EAPsum", full.scores = FALSE)


#dsmsad and age
dsmsadyoung <- subset(comdata, D1>=2 & D1<=4, select=c(DSMSAD1:DSMSAD10))
dsmsadold<- subset(comdata, D1>=5 & D1<=7, select=c(DSMSAD1:DSMSAD10))
comb5young <- subset(comdata, D1>=2 & D1<=4, select=c(SADBNK1:SADBNK26,DSMSAD1:DSMSAD10))
comb5old <- subset(comdata, D1>=5 & D1<=7, select=c(SADBNK1:SADBNK26,DSMSAD1:DSMSAD10))

#young
##SAD item bank and dsmsad
## generating dsmsad irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb5young,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

pardsm<-mirt(comb5young,1, itemtype="graded", pars=sv2)
coef(pardsm, IRTpars=TRUE, simplify=TRUE)
coefdsm<-coef(pardsm, as.data.frame=TRUE)

##calculate dsmsad cross-walk tables
svdsm <-mirt(dsmsadyoung,1, itemtype="graded", pars="values")
sv2dsm<-svdsm
sv2dsm$value<-c(coefdsm[131:182,])

sv2dsm$est<- FALSE

dsmwalk<-mirt(dsmsadyoung,1, itemtype="graded", pars=sv2dsm)
coef(dsmwalk, simplify=TRUE)

dsmeap<-fscores(dsmwalk, method="EAP", full.scores.SE=TRUE)
dsmeapsum<-fscores(dsmwalk, method="EAPsum", full.scores = FALSE)

#old
##SAD item bank and dsmsad
## generating dsmsad irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb5old,1, itemtype="graded", pars="values")
sv2<-sv
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

pardsm<-mirt(comb5old,1, itemtype="graded", pars=sv2)
coef(pardsm, IRTpars=TRUE, simplify=TRUE)
coefdsm<-coef(pardsm, as.data.frame=TRUE)

##calculate dsmsad cross-walk tables
svdsm <-mirt(dsmsadold,1, itemtype="graded", pars="values")
sv2dsm<-svdsm
sv2dsm$value<-c(coefdsm[131:182,])

sv2dsm$est<- FALSE

dsmwalk<-mirt(dsmsadold,1, itemtype="graded", pars=sv2dsm)
coef(dsmwalk, simplify=TRUE)

dsmeap<-fscores(dsmwalk, method="EAP", full.scores.SE=TRUE)
dsmeapsum<-fscores(dsmwalk, method="EAPsum", full.scores = FALSE)




