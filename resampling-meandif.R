#importing the response data
comdata <- read.csv("H:/CRE/NHMRC Early Career Fellowship/Common metrics/papers/social anxiety equating/r-data2.csv")

#subsetting data
comb1 <- subset(comdata, select=c(SADBNK1:SADBNK26,LSASr1:LSASr24))
comb2 <- subset(comdata, select=c(SADBNK1:SADBNK26,SPIN1:SPIN17))
comb3 <- subset(comdata, select=c(SADBNK1:SADBNK26,SIAS1:SIAS20))
comb4 <- subset(comdata, select=c(SADBNK1:SADBNK26,SPS1:SPS20))
comb5 <- subset(comdata, select=c(SADBNK1:SADBNK26,DSMSAD1:DSMSAD10))
comb6<- subset(comdata, select=c(SADBNK1:SADBNK26,SIAS2, SIAS4, SIAS6, SIAS8, SIAS10, SIAS13))
comb7<- subset(comdata, select=c(SADBNK1:SADBNK26,SPS4, SPS7, SPS8, SPS15, SPS16, SPS17))
comb8<- subset(comdata, select=c(SADBNK1:SADBNK26,SIAS2, SIAS4, SIAS6, SIAS8, SIAS10, SIAS13,SPS4, SPS7, SPS8, SPS15, SPS16, SPS17))
comb9<- subset(comdata, select=c(SADBNK1:SADBNK26, SPIN6,SPIN9,SPIN15))
lsas <- subset(comdata, select=c(LSASr1:LSASr24))
spin <- subset(comdata, select=c(SPIN1:SPIN17))
sias <- subset(comdata, select=c(SIAS1:SIAS20))
sias6<- subset(comdata, select=c(SIAS2, SIAS4, SIAS6, SIAS8, SIAS10, SIAS13))
sps <- subset(comdata, select=c(SPS1:SPS20))
sps6<- subset(comdata, select=c(SPS4, SPS7, SPS8, SPS15, SPS16, SPS17))
siassps12<-subset(comdata, select=c(SIAS2, SIAS4, SIAS6, SIAS8, SIAS10, SIAS13,SPS4, SPS7, SPS8, SPS15, SPS16, SPS17))
dsmsad <- subset(comdata, select=c(DSMSAD1:DSMSAD10))
minispin<-subset(comdata, select=c(SPIN6,SPIN9,SPIN15))

###SADBANK PARAMETERS AND MIRT MODEL
#reading in the sad bank parameters estimated from previous community-based sample
sadpars <- read.csv("H:/CRE/NHMRC Early Career Fellowship/Common metrics/papers/social anxiety equating/sadpars.csv")
rownames(sadpars)<-sadpars[,1]
sadpars[,1]<- NULL
sadpars<-as.matrix(sadpars)

#generating the mirt model using the sadbank parameters
sbankmirtmod <- generate.mirt_object(sadpars, itemtype="graded")

###LSAS RESAMPLING
## generating LSAS irt parameters using fixed SAD item bank parameters 
sv <-mirt(comb1,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE
parlsas<-mirt(comb1,1, itemtype="graded", pars=sv2)
coef(parlsas, IRTpars=TRUE, simplify=TRUE)
coeflsas<-coef(parlsas, as.data.frame=TRUE)

svlsas <-mirt(comb1[,27:50],1, itemtype="graded", pars="values")
sv2lsas<-svlsas
sv2lsas$value<-c(coeflsas[131:228,])

sv2lsas$est<- FALSE

lsaswalk<-mirt(lsas,1, itemtype="graded", pars=sv2lsas)

lsaseapsum<-fscores(lsaswalk, method="EAPsum", full.scores = FALSE)


#using those estimated lsas parameters from fixed calibration for form mirt model
lsaspars<- data.frame(a1= c(sv2lsas$value[sv2lsas$name == 'a1']),
                      d1= c(sv2lsas$value[sv2lsas$name == 'd1']),
                      d2= c(sv2lsas$value[sv2lsas$name == 'd2']),
                      d3= c(sv2lsas$value[sv2lsas$name == 'd3']))
lsasmirtmod <- generate.mirt_object(lsaspars, itemtype="graded")



##LSAS resampling for mean difference and ICC
MAT=matrix(nrow=100000,ncol=2)
L=1
B <-500
for (b in 1:B) {
  i=sample (1:1052,200, replace=FALSE)
  comb12=comb1[i,]

lsaseap<-fscores(lsasmirtmod, method="EAP", response.pattern=comb12[,27:50], full.scores.SE=TRUE)

socbeap<-fscores(sbankmirtmod, method="EAP", response.pattern=comb12[,1:26], full.scores.SE=TRUE)

mdif<-mean(socbeap[,27]-lsaseap[,25])

lsastot<-as.data.frame(rowSums(comb12[,27:50]))
rawlsas<-lsaseapsum$Sum.Scores
lsastheta<-lsaseapsum$Theta
recodemat<-data.frame(rawlsas,lsastheta)

idx<-match(lsastot[,1], recodemat$rawlsas)
lsas_cross<-as.matrix(cbind(lsastot, recodemat[idx,-1]))

mcross<-mean(socbeap[,27]-lsas_cross[,2])


MAT [L,1]=as.matrix(mdif)
MAT [L,2]=as.matrix(mcross)

L=L+1
}

MAT2=as.data.frame(MAT[complete.cases(MAT),])

mean(MAT2$V1)
min(MAT2$V1)
max(MAT2$V1)
SD(MAT2$V1)

mean(MAT2$V2)
min(MAT2$V2)
max(MAT2$V2)
SD(MAT2$V2)




###SPIN RESAMPLING
## generating SPIN irt parameters using fixed SAD item bank parameters 
sv <-mirt(comb2,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE
parspin<-mirt(comb2,1, itemtype="graded", pars=sv2)
coef(parspin, IRTpars=TRUE, simplify=TRUE)
coefspin<-coef(parspin, as.data.frame=TRUE)

svspin <-mirt(spin,1, itemtype="graded", pars="values")
sv2spin<-svspin
sv2spin$value<-c(coefspin[131:217,])

sv2spin$est<- FALSE

spinwalk<-mirt(spin,1, itemtype="graded", pars=sv2spin)
spineapsum<-fscores(spinwalk, method="EAPsum", full.scores = FALSE)


#using those estimated lsas parameters from fixed calibration for form mirt model
spinpars<- data.frame(a1= c(sv2spin$value[sv2spin$name == 'a1']),
                      d1= c(sv2spin$value[sv2spin$name == 'd1']),
                      d2= c(sv2spin$value[sv2spin$name == 'd2']),
                      d3= c(sv2spin$value[sv2spin$name == 'd3']),
                      d4= c(sv2spin$value[sv2spin$name == 'd4']))
spinmirtmod <- generate.mirt_object(spinpars, itemtype="graded")


##SPIN resampling for mean difference and ICC
MAT=matrix(nrow=100000,ncol=2)
L=1
B <-500
for (b in 1:B) {
  i=sample (1:1052,200, replace=FALSE)
  comb22=comb2[i,]
  
  spineap<-fscores(spinmirtmod, method="EAP", response.pattern=comb22[,27:43], full.scores.SE=TRUE)
  
  socbeap<-fscores(sbankmirtmod, method="EAP", response.pattern=comb22[,1:26], full.scores.SE=TRUE)
  
  mdif<-mean(socbeap[,27]-spineap[,18])
  
  spintot<-as.data.frame(rowSums(comb22[,27:43]))
  rawspin<-spineapsum$Sum.Scores
  spintheta<-spineapsum$Theta
  recodemat<-data.frame(rawspin,spintheta)
  
  idx<-match(spintot[,1], recodemat$rawspin)
  spin_cross<-as.matrix(cbind(spintot, recodemat[idx,-1]))
  
  mcross<-mean(socbeap[,27]-spin_cross[,2])
  
  
  MAT [L,1]=as.matrix(mdif)
  MAT [L,2]=as.matrix(mcross)
  
  L=L+1
}

MAT2=as.data.frame(MAT[complete.cases(MAT),])

mean(MAT2$V1)
min(MAT2$V1)
max(MAT2$V1)
SD(MAT2$V1)

mean(MAT2$V2)
min(MAT2$V2)
max(MAT2$V2)
SD(MAT2$V2)




#### SIAS RESAMPLING
## generating SIAS irt parameters using fixed SAD item bank parameters 
sv <-mirt(comb3,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE
parsias<-mirt(comb3,1, itemtype="graded", pars=sv2)
coef(parsias, IRTpars=TRUE, simplify=TRUE)
coefsias<-coef(parsias, as.data.frame=TRUE)

svsias <-mirt(sias,1, itemtype="graded", pars="values")
sv2sias<-svsias
sv2sias$value<-c(coefsias[131:232,])

sv2sias$est<- FALSE

siaswalk<-mirt(sias,1, itemtype="graded", pars=sv2sias)
siaseapsum<-fscores(siaswalk, method="EAPsum", full.scores = FALSE)


#using those estimated lsas parameters from fixed calibration for form mirt model
siaspars<- data.frame(a1= c(sv2sias$value[sv2sias$name == 'a1']),
                      d1= c(sv2sias$value[sv2sias$name == 'd1']),
                      d2= c(sv2sias$value[sv2sias$name == 'd2']),
                      d3= c(sv2sias$value[sv2sias$name == 'd3']),
                      d4= c(sv2sias$value[sv2sias$name == 'd4']))
siasmirtmod <- generate.mirt_object(siaspars, itemtype="graded")


##SPIN resampling for mean difference and ICC
MAT=matrix(nrow=100000,ncol=2)
L=1
B <-500
for (b in 1:B) {
  i=sample (1:1052,200, replace=FALSE)
  comb32=comb3[i,]
  
  siaseap<-fscores(siasmirtmod, method="EAP", response.pattern=comb32[,27:46], full.scores.SE=TRUE)
  
  socbeap<-fscores(sbankmirtmod, method="EAP", response.pattern=comb32[,1:26], full.scores.SE=TRUE)
  
  mdif<-mean(socbeap[,27]-siaseap[,21])
  
  siastot<-as.data.frame(rowSums(comb32[,27:46]))
  rawsias<-siaseapsum$Sum.Scores
  siastheta<-siaseapsum$Theta
  recodemat<-data.frame(rawsias,siastheta)
  
  idx<-match(siastot[,1], recodemat$rawsias)
  sias_cross<-as.matrix(cbind(siastot, recodemat[idx,-1]))
  
  mcross<-mean(socbeap[,27]-sias_cross[,2])
  
  
  MAT [L,1]=as.matrix(mdif)
  MAT [L,2]=as.matrix(mcross)
  
  L=L+1
}

MAT2=as.data.frame(MAT[complete.cases(MAT),])

mean(MAT2$V1)
min(MAT2$V1)
max(MAT2$V1)
SD(MAT2$V1)

mean(MAT2$V2)
min(MAT2$V2)
max(MAT2$V2)
SD(MAT2$V2)




#### SPS RESAMPLING
## generating SPS irt parameters using fixed SAD item bank parameters 
sv <-mirt(comb4,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE
parsps<-mirt(comb4,1, itemtype="graded", pars=sv2)
coef(parsps, IRTpars=TRUE, simplify=TRUE)
coefsps<-coef(parsps, as.data.frame=TRUE)

svsps <-mirt(sps,1, itemtype="graded", pars="values")
sv2sps<-svsps
sv2sps$value<-c(coefsps[131:232,])

sv2sps$est<- FALSE

spswalk<-mirt(sps,1, itemtype="graded", pars=sv2sps)
spseapsum<-fscores(spswalk, method="EAPsum", full.scores = FALSE)


#using those estimated lsas parameters from fixed calibration for form mirt model
spspars<- data.frame(a1= c(sv2sps$value[sv2sps$name == 'a1']),
                      d1= c(sv2sps$value[sv2sps$name == 'd1']),
                      d2= c(sv2sps$value[sv2sps$name == 'd2']),
                      d3= c(sv2sps$value[sv2sps$name == 'd3']),
                      d4= c(sv2sps$value[sv2sps$name == 'd4']))
spsmirtmod <- generate.mirt_object(spspars, itemtype="graded")


##SPS resampling for mean difference and ICC
MAT=matrix(nrow=100000,ncol=2)
L=1
B <-500
for (b in 1:B) {
  i=sample (1:1052,50, replace=FALSE)
  comb42=comb4[i,]
  
  spseap<-fscores(spsmirtmod, method="EAP", response.pattern=comb42[,27:46], full.scores.SE=TRUE)
  
  socbeap<-fscores(sbankmirtmod, method="EAP", response.pattern=comb42[,1:26], full.scores.SE=TRUE)
  
  mdif<-mean(socbeap[,27]-spseap[,21])
  
  spstot<-as.data.frame(rowSums(comb42[,27:46]))
  rawsps<-spseapsum$Sum.Scores
  spstheta<-spseapsum$Theta
  recodemat<-data.frame(rawsps,spstheta)
  
  idx<-match(spstot[,1], recodemat$rawsps)
  sps_cross<-as.matrix(cbind(spstot, recodemat[idx,-1]))
  
  mcross<-mean(socbeap[,27]-sps_cross[,2])
  
  
  MAT [L,1]=as.matrix(mdif)
  MAT [L,2]=as.matrix(mcross)
  
  L=L+1
}

MAT2=as.data.frame(MAT[complete.cases(MAT),])

mean(MAT2$V1)
min(MAT2$V1)
max(MAT2$V1)
SD(MAT2$V1)

mean(MAT2$V2)
min(MAT2$V2)
max(MAT2$V2)
SD(MAT2$V2)



#### DSMSAD RESAMPLING
## generating DSMSAD irt parameters using fixed SAD item bank parameters 
sv <-mirt(comb5,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE
pardsm<-mirt(comb5,1, itemtype="graded", pars=sv2)
coef(pardsm, IRTpars=TRUE, simplify=TRUE)
coefdsm<-coef(pardsm, as.data.frame=TRUE)

svdsm <-mirt(dsmsad,1, itemtype="graded", pars="values")
sv2dsm<-svdsm
sv2dsm$value<-c(coefdsm[131:182,])

sv2dsm$est<- FALSE

dsmwalk<-mirt(dsmsad,1, itemtype="graded", pars=sv2dsm)
dsmeapsum<-fscores(dsmwalk, method="EAPsum", full.scores = FALSE)


#using those estimated lsas parameters from fixed calibration for form mirt model
dsmpars<- data.frame(a1= c(sv2dsm$value[sv2dsm$name == 'a1']),
                     d1= c(sv2dsm$value[sv2dsm$name == 'd1']),
                     d2= c(sv2dsm$value[sv2dsm$name == 'd2']),
                     d3= c(sv2dsm$value[sv2dsm$name == 'd3']),
                     d4= c(sv2dsm$value[sv2dsm$name == 'd4']))
dsmmirtmod <- generate.mirt_object(dsmpars, itemtype="graded")


##dsmsad resampling for mean difference and ICC
MAT=matrix(nrow=100000,ncol=2)
L=1
B <-500
for (b in 1:B) {
  i=sample (1:1052,200, replace=FALSE)
  comb52=comb5[i,]
  
  dsmeap<-fscores(dsmmirtmod, method="EAP", response.pattern=comb52[,27:36], full.scores.SE=TRUE)
  
  socbeap<-fscores(sbankmirtmod, method="EAP", response.pattern=comb52[,1:26], full.scores.SE=TRUE)
  
  mdif<-mean(socbeap[,27]-dsmeap[,11])
  
  dsmtot<-as.data.frame(rowSums(comb52[,27:36]))
  rawdsm<-dsmeapsum$Sum.Scores
  dsmtheta<-dsmeapsum$Theta
  recodemat<-data.frame(rawdsm,dsmtheta)
  
  idx<-match(dsmtot[,1], recodemat$rawdsm)
  dsm_cross<-as.matrix(cbind(dsmtot, recodemat[idx,-1]))
  
  mcross<-mean(socbeap[,27]-dsm_cross[,2])
  
  
  MAT [L,1]=as.matrix(mdif)
  MAT [L,2]=as.matrix(mcross)
  
  L=L+1
}

MAT2=as.data.frame(MAT[complete.cases(MAT),])

mean(MAT2$V1)
min(MAT2$V1)
max(MAT2$V1)
SD(MAT2$V1)

mean(MAT2$V2)
min(MAT2$V2)
max(MAT2$V2)
SD(MAT2$V2)




###SIAS6 RESAMPLING
## generating SIAS6 irt parameters using fixed SAD item bank parameters 
sv <-mirt(comb6,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE
parsias6<-mirt(comb6,1, itemtype="graded", pars=sv2)
coef(parsias6, IRTpars=TRUE, simplify=TRUE)
coefsias6<-coef(parsias6, as.data.frame=TRUE)

svsias6 <-mirt(sias6,1, itemtype="graded", pars="values")
sv2sias6<-svsias6
sv2sias6$value<-c(coefsias6[131:162,])

sv2sias6$est<- FALSE

sias6walk<-mirt(sias6,1, itemtype="graded", pars=sv2sias6)
sias6eapsum<-fscores(sias6walk, method="EAPsum", full.scores = FALSE)


#using those estimated lsas parameters from fixed calibration for form mirt model
sias6pars<- data.frame(a1= c(sv2sias6$value[sv2sias6$name == 'a1']),
                     d1= c(sv2sias6$value[sv2sias6$name == 'd1']),
                     d2= c(sv2sias6$value[sv2sias6$name == 'd2']),
                     d3= c(sv2sias6$value[sv2sias6$name == 'd3']),
                     d4= c(sv2sias6$value[sv2sias6$name == 'd4']))
sias6mirtmod <- generate.mirt_object(sias6pars, itemtype="graded")


##sias6 resampling for mean difference and ICC
MAT=matrix(nrow=100000,ncol=2)
L=1
B <-500
for (b in 1:B) {
  i=sample (1:1052,100, replace=FALSE)
  comb62=comb6[i,]
  
  sias6eap<-fscores(sias6mirtmod, method="EAP", response.pattern=comb62[,27:32], full.scores.SE=TRUE)
  
  socbeap<-fscores(sbankmirtmod, method="EAP", response.pattern=comb62[,1:26], full.scores.SE=TRUE)
  
  mdif<-mean(socbeap[,27]-sias6eap[,7])
  
  sias6tot<-as.data.frame(rowSums(comb62[,27:32]))
  rawsias6<-sias6eapsum$Sum.Scores
  sias6theta<-sias6eapsum$Theta
  recodemat<-data.frame(rawsias6,sias6theta)
  
  idx<-match(sias6tot[,1], recodemat$rawsias6)
  sias6_cross<-as.matrix(cbind(sias6tot, recodemat[idx,-1]))
  
  mcross<-mean(socbeap[,27]-sias6_cross[,2])
  
  
  MAT [L,1]=as.matrix(mdif)
  MAT [L,2]=as.matrix(mcross)
  
  L=L+1
}

MAT2=as.data.frame(MAT[complete.cases(MAT),])

mean(MAT2$V1)
min(MAT2$V1)
max(MAT2$V1)
SD(MAT2$V1)

mean(MAT2$V2)
min(MAT2$V2)
max(MAT2$V2)
SD(MAT2$V2)



###SPS6 RESAMPLING
## generating SPS6 irt parameters using fixed SAD item bank parameters 
sv <-mirt(comb7,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE
parsps6<-mirt(comb7,1, itemtype="graded", pars=sv2)
coef(parsps6, IRTpars=TRUE, simplify=TRUE)
coefsps6<-coef(parsps6, as.data.frame=TRUE)

svsps6 <-mirt(sps6,1, itemtype="graded", pars="values")
sv2sps6<-svsps6
sv2sps6$value<-c(coefsps6[131:162,])

sv2sps6$est<- FALSE

sps6walk<-mirt(sps6,1, itemtype="graded", pars=sv2sps6)
sps6eapsum<-fscores(sps6walk, method="EAPsum", full.scores = FALSE)


#using those estimated lsas parameters from fixed calibration for form mirt model
sps6pars<- data.frame(a1= c(sv2sps6$value[sv2sps6$name == 'a1']),
                       d1= c(sv2sps6$value[sv2sps6$name == 'd1']),
                       d2= c(sv2sps6$value[sv2sps6$name == 'd2']),
                       d3= c(sv2sps6$value[sv2sps6$name == 'd3']),
                       d4= c(sv2sps6$value[sv2sps6$name == 'd4']))
sps6mirtmod <- generate.mirt_object(sps6pars, itemtype="graded")


##sps6 resampling for mean difference and ICC
MAT=matrix(nrow=100000,ncol=2)
L=1
B <-500
for (b in 1:B) {
  i=sample (1:1052,200, replace=FALSE)
  comb72=comb7[i,]
  
  sps6eap<-fscores(sps6mirtmod, method="EAP", response.pattern=comb72[,27:32], full.scores.SE=TRUE)
  
  socbeap<-fscores(sbankmirtmod, method="EAP", response.pattern=comb72[,1:26], full.scores.SE=TRUE)
  
  mdif<-mean(socbeap[,27]-sps6eap[,7])
  
  sps6tot<-as.data.frame(rowSums(comb72[,27:32]))
  rawsps6<-sps6eapsum$Sum.Scores
  sps6theta<-sps6eapsum$Theta
  recodemat<-data.frame(rawsps6,sps6theta)
  
  idx<-match(sps6tot[,1], recodemat$rawsps6)
  sps6_cross<-as.matrix(cbind(sps6tot, recodemat[idx,-1]))
  
  mcross<-mean(socbeap[,27]-sps6_cross[,2])
  
  
  MAT [L,1]=as.matrix(mdif)
  MAT [L,2]=as.matrix(mcross)
  
  L=L+1
}

MAT2=as.data.frame(MAT[complete.cases(MAT),])

mean(MAT2$V1)
min(MAT2$V1)
max(MAT2$V1)
SD(MAT2$V1)

mean(MAT2$V2)
min(MAT2$V2)
max(MAT2$V2)
SD(MAT2$V2)



###SPS12 RESAMPLING
## generating SPS12 irt parameters using fixed SAD item bank parameters 
sv <-mirt(comb8,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE
parsps12<-mirt(comb8,1, itemtype="graded", pars=sv2)
coef(parsps12, IRTpars=TRUE, simplify=TRUE)
coefsps12<-coef(parsps12, as.data.frame=TRUE)

svsps12 <-mirt(siassps12,1, itemtype="graded", pars="values")
sv2sps12<-svsps12
sv2sps12$value<-c(coefsps12[131:192,])

sv2sps12$est<- FALSE

sps12walk<-mirt(siassps12,1, itemtype="graded", pars=sv2sps12)
sps12eapsum<-fscores(sps12walk, method="EAPsum", full.scores = FALSE)


#using those estimated lsas parameters from fixed calibration for form mirt model
sps12pars<- data.frame(a1= c(sv2sps12$value[sv2sps12$name == 'a1']),
                      d1= c(sv2sps12$value[sv2sps12$name == 'd1']),
                      d2= c(sv2sps12$value[sv2sps12$name == 'd2']),
                      d3= c(sv2sps12$value[sv2sps12$name == 'd3']),
                      d4= c(sv2sps12$value[sv2sps12$name == 'd4']))
sps12mirtmod <- generate.mirt_object(sps12pars, itemtype="graded")


##sps12 resampling for mean difference and ICC
MAT=matrix(nrow=100000,ncol=2)
L=1
B <-500
for (b in 1:B) {
  i=sample (1:1052,200, replace=FALSE)
  comb82=comb8[i,]
  
  sps12eap<-fscores(sps12mirtmod, method="EAP", response.pattern=comb82[,27:38], full.scores.SE=TRUE)
  
  socbeap<-fscores(sbankmirtmod, method="EAP", response.pattern=comb82[,1:26], full.scores.SE=TRUE)
  
  mdif<-mean(socbeap[,27]-sps12eap[,13])
  
  sps12tot<-as.data.frame(rowSums(comb82[,27:38]))
  rawsps12<-sps12eapsum$Sum.Scores
  sps12theta<-sps12eapsum$Theta
  recodemat<-data.frame(rawsps12,sps12theta)
  
  idx<-match(sps12tot[,1], recodemat$rawsps12)
  sps12_cross<-as.matrix(cbind(sps12tot, recodemat[idx,-1]))
  
  mcross<-mean(socbeap[,27]-sps12_cross[,2])
  
  
  MAT [L,1]=as.matrix(mdif)
  MAT [L,2]=as.matrix(mcross)
  
  L=L+1
}

MAT2=as.data.frame(MAT[complete.cases(MAT),])

mean(MAT2$V1)
min(MAT2$V1)
max(MAT2$V1)
SD(MAT2$V1)

mean(MAT2$V2)
min(MAT2$V2)
max(MAT2$V2)
SD(MAT2$V2)



###MINISPIN RESAMPLING
## generating MINISPIN irt parameters using fixed SAD item bank parameters 
sv <-mirt(comb9,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE
parsmpin<-mirt(comb9,1, itemtype="graded", pars=sv2)
coef(parsmpin, IRTpars=TRUE, simplify=TRUE)
coefmspin<-coef(parsmpin, as.data.frame=TRUE)

svmspin <-mirt(minispin,1, itemtype="graded", pars="values")
sv2mspin<-svmspin
sv2mspin$value<-c(coefmspin[131:147,])

sv2mspin$est<- FALSE

mspinwalk<-mirt(minispin,1, itemtype="graded", pars=sv2mspin)
mspineapsum<-fscores(mspinwalk, method="EAPsum", full.scores = FALSE)


#using those estimated lsas parameters from fixed calibration for form mirt model
mspinpars<- data.frame(a1= c(sv2mspin$value[sv2mspin$name == 'a1']),
                       d1= c(sv2mspin$value[sv2mspin$name == 'd1']),
                       d2= c(sv2mspin$value[sv2mspin$name == 'd2']),
                       d3= c(sv2mspin$value[sv2mspin$name == 'd3']),
                       d4= c(sv2mspin$value[sv2mspin$name == 'd4']))
mspinmirtmod <- generate.mirt_object(mspinpars, itemtype="graded")


##sps12 resampling for mean difference and ICC
MAT=matrix(nrow=100000,ncol=2)
L=1
B <-500
for (b in 1:B) {
  i=sample (1:1052,200, replace=FALSE)
  comb92=comb9[i,]
  
  mspineap<-fscores(mspinmirtmod, method="EAP", response.pattern=comb92[,27:29], full.scores.SE=TRUE)
  
  socbeap<-fscores(sbankmirtmod, method="EAP", response.pattern=comb92[,1:26], full.scores.SE=TRUE)
  
  mdif<-mean(socbeap[,27]-mspineap[,4])
  
  mspintot<-as.data.frame(rowSums(comb92[,27:29]))
  rawmspin<-mspineapsum$Sum.Scores
  mspintheta<-mspineapsum$Theta
  recodemat<-data.frame(rawmspin,mspintheta)
  
  idx<-match(mspintot[,1], recodemat$rawmspin)
  mspin_cross<-as.matrix(cbind(mspintot, recodemat[idx,-1]))
  
  mcross<-mean(socbeap[,27]-mspin_cross[,2])
  
  
  MAT [L,1]=as.matrix(mdif)
  MAT [L,2]=as.matrix(mcross)
  
  L=L+1
}

MAT2=as.data.frame(MAT[complete.cases(MAT),])

mean(MAT2$V1)
min(MAT2$V1)
max(MAT2$V1)
SD(MAT2$V1)

mean(MAT2$V2)
min(MAT2$V2)
max(MAT2$V2)
SD(MAT2$V2)


## load dataset
sophsdata <- read.csv("H:/CRE/NHMRC Early Career Fellowship/Common metrics/papers/social anxiety equating/r-data3.csv")

socbank <- subset(sophsdata, select=c(SADBNK1:SADBNK26))
sophs <- subset(sophsdata, select=c(SPS1_1, SPS2_2, SPS2_3, SPS2_4, SPS2_5))

combs1 <- subset(sophsdata, select=c(SADBNK1:SADBNK26,SPS1_1, SPS2_2, SPS2_3, SPS2_4, SPS2_5))



###SOPHS RESAMPLING
## generating sophs irt parameters using fixed SAD item bank parameters 
sv <-mirt(combs1,1, itemtype="graded", pars="values")
sv2<-sv
sadvals = as.numeric(t(sadpars))
sv2$value[1:length(sadvals)] = sadvals
sv2$est[1:length(sadvals)] = FALSE
parssophs<-mirt(combs1,1, itemtype="graded", pars=sv2)
coef(parssophs, IRTpars=TRUE, simplify=TRUE)
coefsophs<-coef(parssophs, as.data.frame=TRUE)

svsophs <-mirt(sophs,1, itemtype="graded", pars="values")
sv2sophs<-svsophs
sv2sophs$value<-c(coefsophs[131:157,])

sv2sophs$est<- FALSE

sophswalk<-mirt(sophs,1, itemtype="graded", pars=sv2sophs)
sophseapsum<-fscores(sophswalk, method="EAPsum", full.scores = FALSE)


#using those estimated lsas parameters from fixed calibration for form mirt model
sophspars<- data.frame(a1= c(sv2sophs$value[sv2sophs$name == 'a1']),
                       d1= c(sv2sophs$value[sv2sophs$name == 'd1']),
                       d2= c(sv2sophs$value[sv2sophs$name == 'd2']),
                       d3= c(sv2sophs$value[sv2sophs$name == 'd3']),
                       d4= c(sv2sophs$value[sv2sophs$name == 'd4']))
sophsmirtmod <- generate.mirt_object(sophspars, itemtype="graded")


##sps12 resampling for mean difference and ICC
MAT=matrix(nrow=100000,ncol=2)
L=1
B <-500
for (b in 1:B) {
  i=sample (1:1052,200, replace=FALSE)
  combs12=combs1[i,]
  
  sophseap<-fscores(sophsmirtmod, method="EAP", response.pattern=combs12[,27:31], full.scores.SE=TRUE)
  
  socbeap<-fscores(sbankmirtmod, method="EAP", response.pattern=combs12[,1:26], full.scores.SE=TRUE)
  
  mdif<-mean(socbeap[,27]-sophseap[,6])
  
  sophstot<-as.data.frame(rowSums(combs12[,27:31]))
  rawsophs<-sophseapsum$Sum.Scores
  sophstheta<-sophseapsum$Theta
  recodemat<-data.frame(rawsophs,sophstheta)
  
  idx<-match(sophstot[,1], recodemat$rawsophs)
  sophs_cross<-as.matrix(cbind(sophstot, recodemat[idx,-1]))
  
  mcross<-mean(socbeap[,27]-sophs_cross[,2])
  
  
  MAT [L,1]=as.matrix(mdif)
  MAT [L,2]=as.matrix(mcross)
  
  L=L+1
}

MAT2=as.data.frame(MAT[complete.cases(MAT),])

mean(MAT2$V1)
min(MAT2$V1)
max(MAT2$V1)
SD(MAT2$V1)

mean(MAT2$V2)
min(MAT2$V2)
max(MAT2$V2)
SD(MAT2$V2)

