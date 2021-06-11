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
sophsdata <- read.csv("H:/CRE/NHMRC Early Career Fellowship/Common metrics/papers/social anxiety equating/r-data3.csv")

socbank <- subset(sophsdata, select=c(SADBNK1:SADBNK26))
sophs <- subset(sophsdata, select=c(SPS1_1, SPS2_2, SPS2_3, SPS2_4, SPS2_5))

comb1 <- subset(sophsdata, select=c(SADBNK1:SADBNK26,SPS1_1, SPS2_2, SPS2_3, SPS2_4, SPS2_5))


##IRT item fit and residual matrices to check for local dependence
comb1uni<-mirt(comb1,1, itemtype="graded")
coef(comb1uni, simplify=TRUE)
ifit1<-itemfit(comb1uni)
residmat1<-M2(comb1uni,residmat=TRUE)


##fixed-anchor IRT calibration
sadpars <- read.csv("H:/CRE/NHMRC Early Career Fellowship/Common metrics/papers/social anxiety equating/sadpars.csv")
rownames(sadpars)<-sadpars[,1]
sadpars[,1]<- NULL
sadpars<-as.matrix(sadpars)

##SAD item bank and SOPHS
## generating SOPHS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb1,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parsophs<-mirt(comb1,1, itemtype="graded", pars=sv2)
coef(parsophs, IRTpars=TRUE, simplify=TRUE)
coefsophs<-coef(parsophs, as.data.frame=TRUE)
eap<-fscores(parsophs, method="EAP")
sum(eap[1:3175,])/3175

##calculate sophs cross-walk tables
svsophs <-mirt(sophs,1, itemtype="graded", pars="values")
sv2sophs<-svsophs
sv2sophs$value<-c(coefsophs[131:157,])

sv2sophs$est<- FALSE

sophswalk<-mirt(sophs,1, itemtype="graded", pars=sv2sophs)
coef(sophswalk)

infosophsfixed<-plot(sophswalk, type="info")
infofixedsophs<-as.data.table(cbind(infosophsfixed$panel.args[[1]][[1]],infosophsfixed$panel.args[[1]][[2]]))

tccsophsfixed<-plot(sophswalk, type="score")
scorefixedsophs<-as.data.table(cbind(tccsophsfixed$panel.args[[1]][[1]],tccsophsfixed$panel.args[[1]][[2]]))

sophseap<-fscores(sophswalk, method="EAP", full.scores.SE=TRUE)
sum(sophseap[1:3175,1])/3175
sophseapsum<-fscores(sophswalk, method="EAPsum", full.scores = FALSE)

##calculate socbank cross-walk tables
svsocb <-mirt(socbank,1, itemtype="graded", pars="values")
sv2socb<-svsocb
sv2socb$value[1:length(sadvals)] = sadvals

sv2socb$est<- FALSE

socbwalk<-mirt(socbank,1, itemtype="graded", pars=sv2socb)
coef(socbwalk, simplify=TRUE)

socbeap<-fscores(socbwalk, method="EAP", full.scores.SE=TRUE)
sum(socbeap[1:3175,1])/3175
socbeapsum<-fscores(socbwalk, method="EAPsum", full.scores = FALSE)



##plink IRT separate calibration equating
##SAD item bank and sophs
##generating sadbank coefficients
sv <-mirt(comb1,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE
mirtsophs<-mirt(comb1,1, itemtype="graded", pars=sv2)
mirtceof<-coef(mirtsophs, IRTpars=TRUE, simplify=TRUE)
sadparsirt<-mirtceof$items[1:26,]

## generating separeately calibrated items for sadbank and sophs
comb1uni<-mirt(comb1,1, itemtype="graded")
sadparsophs<-coef(comb1uni, simplify=TRUE, IRTpars=TRUE)
sadparsophs<-as.matrix(sadparsophs$items)
##specific number of categories for both item parameters to be linked
catsad<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
catsophs<-c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
##specific the model to be used for both sets of item paramters
pmsad<-as.poly.mod(26,c("grm"),list(c(1:26)))
pmsophs<-as.poly.mod(31,c("grm"),list(c(1:31)))
##indicate the number of common items between the two sets of item parameters
commonsophs<-matrix(c(1:26, 1:26),26,2)
##combine all info into a single object
parssophs<-as.irt.pars(x=list(sadparsirt,sadparsophs), common=list(commonsophs), cat=list(catsad, catsophs), poly.mod=list(pmsad, pmsophs),grp.names=c("sadbank","sophs"))

##stocking lord
plinksophs<-plink(parssophs, rescale="SL")
summary(plinksophs)
reparsophsSL<-link.pars(plinksophs)
reparsophsSL<-reparsophsSL$sophs[27:31,1:5]
d1=-1*(reparsophsSL[,1]*reparsophsSL[,2])
d2=-1*(reparsophsSL[,1]*reparsophsSL[,3])
d3=-1*(reparsophsSL[,1]*reparsophsSL[,4])
d4=-1*(reparsophsSL[,1]*reparsophsSL[,5])
reparsophsSL2<-cbind(reparsophsSL[,1],d1,d2,d3,d4)
colnames(reparsophsSL2)<-c("a1","d1","d2","d3","d4")
object_reparsophsSL<-generate.mirt_object(reparsophsSL2, "graded")
tccsophsSF<-plot(object_reparsophsSL, type="score")
tccsophsSF
scoreSFsophs<-as.data.table(cbind(tccsophsSF$panel.args[[1]][[1]],tccsophsSF$panel.args[[1]][[2]]))
#mean/mean
plinksophs<-plink(parssophs, rescale="MM")
summary(plinksophs)
reparsophsMM<-link.pars(plinksophs)
reparsophsMM<-reparsophsMM$sophs[27:31,1:5]
d1=-1*(reparsophsMM[,1]*reparsophsMM[,2])
d2=-1*(reparsophsMM[,1]*reparsophsMM[,3])
d3=-1*(reparsophsMM[,1]*reparsophsMM[,4])
d4=-1*(reparsophsMM[,1]*reparsophsMM[,5])
reparsophsMM2<-cbind(reparsophsMM[,1],d1,d2,d3,d4)
colnames(reparsophsMM2)<-c("a1","d1","d2","d3","d4")
object_reparsophsMM<-generate.mirt_object(reparsophsMM2, "graded")
tccsophsMM<-plot(object_reparsophsMM, type="score")
scoreMMsophs<-as.data.table(cbind(tccsophsMM$panel.args[[1]][[1]],tccsophsMM$panel.args[[1]][[2]]))
#mean/sigma
plinksophs<-plink(parssophs, rescale="MS")
summary(plinksophs)
reparsophsMS<-link.pars(plinksophs)
reparsophsMS<-reparsophsMS$sophs[27:31,1:5]
d1=-1*(reparsophsMS[,1]*reparsophsMS[,2])
d2=-1*(reparsophsMS[,1]*reparsophsMS[,3])
d3=-1*(reparsophsMS[,1]*reparsophsMS[,4])
d4=-1*(reparsophsMS[,1]*reparsophsMS[,5])
reparsphsMS2<-cbind(reparsophsMS[,1],d1,d2,d3,d4)
colnames(reparsphsMS2)<-c("a1","d1","d2","d3","d4")
object_reparsophsMS<-generate.mirt_object(reparsphsMS2, "graded")
tccsophsMS<-plot(object_reparsophsMS, type="score")
scoreMSsophs<-as.data.table(cbind(tccsophsMS$panel.args[[1]][[1]],tccsophsMS$panel.args[[1]][[2]]))
#haebara
plinksophs<-plink(parssophs, rescale="HB")
summary(plinksophs)
reparsophsHB<-link.pars(plinksophs)
reparsophsHB<-reparsophsHB$sophs[27:31,1:5]
d1=-1*(reparsophsHB[,1]*reparsophsHB[,2])
d2=-1*(reparsophsHB[,1]*reparsophsHB[,3])
d3=-1*(reparsophsHB[,1]*reparsophsHB[,4])
d4=-1*(reparsophsHB[,1]*reparsophsHB[,5])
reparsophsHB2<-cbind(reparsophsHB[,1],d1,d2,d3,d4)
colnames(reparsophsHB2)<-c("a1","d1","d2","d3","d4")
object_reparsophsHB<-generate.mirt_object(reparsophsHB2, "graded")
tccsophsHB<-plot(object_reparsophsHB, type="score")
scoreHBsophs<-as.data.table(cbind(tccsophsHB$panel.args[[1]][[1]],tccsophsHB$panel.args[[1]][[2]]))


mean(sophseap[ ,1])
SD(sophseap[ ,1])

cor(socbeap[,1], sophseap[,1])

mean(socbeap[,1]-sophseap[,1])

SD(socbeap[,1]-sophseap[,1])

ICC(as.matrix(cbind(socbeap[,1],sophseap[,1])))


sophstot<-subset(sophsdata, select=c(SOPHStot))
rawsophs<-sophseapsum$Sum.Scores
sophstheta<-sophseapsum$Theta
recodemat<-data.frame(rawsophs,sophstheta)

idx<-match(sophstot$SOPHStot, recodemat$rawsophs)
sophs_cross<-as.matrix(cbind(sophstot, recodemat[idx,-1]))

mean(sophs_cross[,2])
SD(sophs_cross[,2])

cor(socbeap[,1],sophs_cross[,2])
ICC(as.matrix(cbind(socbeap[,1],sophs_cross[,2])))
mean(socbeap[,1]-sophs_cross[,2])
SD(socbeap[,1]-sophs_cross[,2])

bland.altman.plot(socbeap[,1], sophseap[,1], main="SAD-bank and SOPHS", xlab="Mean", ylab="Differences")

comdata <- read.csv("D:/Users/z3226077/OneDrive - UNSW/social anxiety equating/r-data4.csv")


##population invariance testing
##sophs gender
sophsmale <- subset(comdata, D2==1, select=c(SPS1_1, SPS2_2, SPS2_3, SPS2_4, SPS2_5))
sophsfemale <- subset(comdata, D2==2, select=c(SPS1_1, SPS2_2, SPS2_3, SPS2_4, SPS2_5))
comb1male <- subset(comdata, D2==1, select=c(SADBNK1:SADBNK26,SPS1_1, SPS2_2, SPS2_3, SPS2_4, SPS2_5))
comb1female <- subset(comdata, D2==2, select=c(SADBNK1:SADBNK26,SPS1_1, SPS2_2, SPS2_3, SPS2_4, SPS2_5))

##males
##SAD item bank and LSAS
## generating LSAS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb1male,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parsophs<-mirt(comb1male,1, itemtype="graded", pars=sv2)
coef(parsophs, IRTpars=TRUE, simplify=TRUE)
coefsophs<-coef(parsophs, as.data.frame=TRUE)
eap<-fscores(parsophs, method="EAP")
sum(eap[1:648,])/648

##calculate lsas cross-walk tables
svsophs<-mirt(sophsmale,1, itemtype="graded", pars="values")
sv2sophs<-svsophs
sv2sophs$value<-c(coefsophs[131:157,])

sv2sophs$est<- FALSE

sophswalk<-mirt(sophsmale,1, itemtype="graded", pars=sv2sophs)
sophseap<-fscores(sophswalk, method="EAP", full.scores.SE=TRUE)
sum(sophseap[1:648,1])/648
sophseapsum<-fscores(sophswalk, method="EAPsum", full.scores = FALSE)

##females
##SAD item bank and SOPHS
## generating SOPHS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb1female,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parsophs<-mirt(comb1female,1, itemtype="graded", pars=sv2)
coef(parsophs, IRTpars=TRUE, simplify=TRUE)
coefsophs<-coef(parsophs, as.data.frame=TRUE)
eap<-fscores(parsophs, method="EAP")
sum(eap[1:2527,])/2527

##calculate lsas cross-walk tables
svsophs<-mirt(sophsfemale,1, itemtype="graded", pars="values")
sv2sophs<-svsophs
sv2sophs$value<-c(coefsophs[131:157,])

sv2sophs$est<- FALSE

sophswalk<-mirt(sophsfemale,1, itemtype="graded", pars=sv2sophs)
sophseap<-fscores(sophswalk, method="EAP", full.scores.SE=TRUE)
sum(sophseap[1:2527,1])/2527
sophseapsum<-fscores(sophswalk, method="EAPsum", full.scores = FALSE)



## age invariance 
##population invariance testing
##sophs age
sophsyoung <- subset(comdata, D1>=1 & D1<=3, select=c(SPS1_1, SPS2_2, SPS2_3, SPS2_4, SPS2_5))
sophsold <- subset(comdata, D1>=4 & D1<=6, select=c(SPS1_1, SPS2_2, SPS2_3, SPS2_4, SPS2_5))
comb1young <- subset(comdata, D1>=1 & D1<=3, select=c(SADBNK1:SADBNK26,SPS1_1, SPS2_2, SPS2_3, SPS2_4, SPS2_5))
comb1old <- subset(comdata, D1>=4 & D1<=6, select=c(SADBNK1:SADBNK26,SPS1_1, SPS2_2, SPS2_3, SPS2_4, SPS2_5))

##young
##SAD item bank and LSAS
## generating LSAS irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb1young,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parsophs<-mirt(comb1young,1, itemtype="graded", pars=sv2)
coef(parsophs, IRTpars=TRUE, simplify=TRUE)
coefsophs<-coef(parsophs, as.data.frame=TRUE)
eap<-fscores(parsophs, method="EAP")
sum(eap[1:1154,])/1154

##calculate lsas cross-walk tables
svsophs<-mirt(sophsyoung,1, itemtype="graded", pars="values")
sv2sophs<-svsophs
sv2sophs$value<-c(coefsophs[131:157,])

sv2sophs$est<- FALSE

sophswalk<-mirt(sophsyoung,1, itemtype="graded", pars=sv2sophs)
sophseap<-fscores(sophswalk, method="EAP", full.scores.SE=TRUE)
sum(sophseap[1:1154,1])/1154
sophseapsum<-fscores(sophswalk, method="EAPsum", full.scores = FALSE)

##old
##SAD item bank and sophs
## generating sophs irt parameters using fixed SAD item bank parameters from literature
sv <-mirt(comb1old,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE

parsophs<-mirt(comb1old,1, itemtype="graded", pars=sv2)
coef(parsophs, IRTpars=TRUE, simplify=TRUE)
coefsophs<-coef(parsophs, as.data.frame=TRUE)
eap<-fscores(parsophs, method="EAP")
sum(eap[1:2021,])/2021

##calculate sophs cross-walk tables
svsophs<-mirt(sophsold,1, itemtype="graded", pars="values")
sv2sophs<-svsophs
sv2sophs$value<-c(coefsophs[131:157,])

sv2sophs$est<- FALSE

sophswalk<-mirt(sophsold,1, itemtype="graded", pars=sv2sophs)
sophseap<-fscores(sophswalk, method="EAP", full.scores.SE=TRUE)
sum(sophseap[1:2021,1])/2021
sophseapsum<-fscores(sophswalk, method="EAPsum", full.scores = FALSE)

