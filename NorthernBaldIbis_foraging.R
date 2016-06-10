##########################################################################################
########   ANALYSIS OF NORTHERN BALD IBIS AT FEEDING PLACES IN BIRECIK, TURKEY   #########
##########################################################################################
# written by Steffen Oppel, 20 November 2013
# updated 29 Jan 2014 to adjust to new lme4 version, removed "control=list(maxIter=6000)" from lmer function call
# recommended citation:
## Influence of feeding ecology on breeding success of the critically endangered Northern Bald Ibis Geronticus eremita in southern Turkey
## CAN YENIYURT, STEFFEN OPPEL, SÜREYYA ÐSFENDIYAROGLU, GÜLÇIN ÖZKINACI, ITRI LEVENT ERKOL, CHRISTOPHER G. R. BOWDEN
## Bird Conservation International



#####################################################################################################################################################
#############    LOAD RAW DATA FROM GITHUB      ##############################################################################
#####################################################################################################################################################
setInternet2(TRUE)
con <- url("https://github.com/steffenoppel/ibis/blob/master/NBI_input.RData?raw=true")
load(con)
close(con)

rm(list= ls()[!(ls() %in% c('NBI','NBIbreed','NBIcrop','NBI_hab','behav','att_all'))])

### DATA 'NBI': data.frame that lists the breeding output for each colour-ringed individual in each year, plus the proportion of feeding events attended and the number of sightings at various foraging areas
head(NBI)

### DATA 'NBIbreed': subset of NBI that contains only data for breeding individuals
head(NBIbreed)

### DATA 'NBIcrop': data.frame that combines 'NBI' with observations in various crop types per individual and year
head(NBIcrop)

### DATA 'NBI_hab': data.frame that lists the foraging observations for each colour-ringed individual in each year and habitat, with 'n_pecks' indicating how many pecks were necessary to acquire prey 
head(NBI_hab)

### DATA 'behav': data.frame that lists the number of birds displaying a certain behaviour in a flock at a given time and site
head(behav)

### DATA 'att_all': data.frame that lists the attendance of each colour-ringed individual at each unique feeding event in each year
head(att_all)



### LOAD NECESSARY LIBRARIES
library(lme4)
library(MuMIn)
library(AICcmodavg)
library(effects)
library(reshape)
library(Hmisc)





#######################################################################################
#
#  IS BREEDING PROPENSITY AND PERFORMANCE HIGHER FOR BIRDS THAT FEED MORE FREQUENTLY AT FEEDING STATION?
#
#######################################################################################


##### ARE NON-BREEDERS MORE FREQUENTLY SEEN AT SUPPLEMENTARY FEEDING EVENTS THAN BREEDERS? #####
bp<-glmer(breed~propPres+(1|ColourRing), family=binomial, data=NBI)
summary(bp)
aggregate(propPres~breed, NBI, FUN=mean)
aggregate(propPres~breed, NBI, FUN=min)
aggregate(propPres~breed, NBI, FUN=max)



### correlate attendance probability with nesting success
cor.test(NBIbreed$breed_succ,NBIbreed$propPres)

### test whether attendance proportion at feeding events influences number of fledglings, number of hatched eggs, and number of eggs laid
fl<-glmer(fledged~propPres+(1|Year/Nest_ID), family=poisson, data=NBIbreed)
summary(fl)

ha<-glmer(hatched~propPres+(1|Year/Nest_ID), family=poisson, data=NBIbreed)
summary(ha)

bp<-glmer(n_eggs~propPres+(1|Year/Nest_ID), family=poisson, data=NBIbreed)
summary(bp)






#######################################################################################
#
#  IS PRODUCTIVITY RELATED TO FREQUENCY OF OBSERVATION IN PARTICULAR FIELDS?
#
#######################################################################################

### PLAUSIBLE MODELS DESCRIBED IN MANUSCRIPT ###
m1<-glmer(fledged~Fidanlik+(1|Year/Nest_ID), data=NBI,family=poisson)
m2<-glmer(fledged~Mezra+(1|Year/Nest_ID), data=NBI,family=poisson)
m4<-glmer(fledged~Fidanlik+Mezra+Saray+(1|Year/Nest_ID), data=NBI,family=poisson)
m5<-glmer(fledged~Fidanlik+Mezra+(1|Year/Nest_ID), data=NBI,family=poisson)
m6<-glmer(fledged~propPres+(1|Year/Nest_ID), data=NBI,family=poisson)
m7<-glmer(fledged~Fidanlik+Saray+(1|Year/Nest_ID), data=NBI,family=poisson)
m8<-glmer(fledged~1+(1|Year/Nest_ID), data=NBI,family=poisson)

AIC_TABLE<-aictab(cand.set=list(m1,m2,m4,m5,m6,m7,m8),modnames=c('m1','m2','m4','m5','m6','m7','m8'),sort = TRUE, c.hat = 1, second.ord = TRUE, nobs = NULL)
AIC_TABLE



### assess whether tree nursery usage is correlated with presence at supplementary feeding events
cor.test(NBI$Fidanlik,NBI$propPres)
plot(NBI$Fidanlik,NBI$propPres)




################## SIMILAR ANALYSIS RELATING BREEDING PROPENSITY TO FIELD USE ###########


m1<-glmer(breed~Fidanlik+(1|ColourRing), data=NBI,family=binomial)
m2<-glmer(breed~Mezra+(1|ColourRing), data=NBI,family=binomial)
m4<-glmer(breed~Fidanlik+Mezra+Saray+(1|ColourRing), data=NBI,family=binomial)
m5<-glmer(breed~Fidanlik+Mezra+(1|ColourRing), data=NBI,family=binomial)
m6<-glmer(breed~propPres+(1|ColourRing), data=NBI,family=binomial)
m7<-glmer(breed~Fidanlik+Saray+(1|ColourRing), data=NBI,family=binomial)
m8<-glmer(breed~1+(1|ColourRing), data=NBI,family=binomial)

AIC_TABLE<-aictab(cand.set=list(m1,m2,m4,m5,m6,m7,m8),modnames=c('m1','m2','m4','m5','m6','m7','m8'),sort = TRUE, c.hat = 1, second.ord = TRUE, nobs = NULL)
AIC_TABLE





#######################################################################################
#
#  IS PRODUCTIVITY RELATED TO FREQUENCY OF OBSERVATION IN PARTICULAR CROPS?
#
#######################################################################################

NBIcropbreed<-NBIcrop[NBIcrop$breed==1,]

m1<-glmer(fledged~Eggplant+(1|Nest_ID), data=NBIcrop,family=poisson)
m2<-glmer(fledged~Fallow+(1|Nest_ID), data=NBIcrop,family=poisson)
m3<-glmer(fledged~Letture+(1|Nest_ID), data=NBIcrop,family=poisson)
m4<-glmer(fledged~Manure+(1|Nest_ID), data=NBIcrop,family=poisson)
m5<-glmer(fledged~Mint+(1|Nest_ID), data=NBIcrop,family=poisson)
m6<-glmer(fledged~Pasture+(1|Nest_ID), data=NBIcrop,family=poisson)
m7<-glmer(fledged~Seedling+(1|Nest_ID), data=NBIcrop,family=poisson)
m8<-glmer(fledged~1+(1|Nest_ID), data=NBIcrop,family=poisson)

AIC_TABLE<-aictab(cand.set=list(m1,m2,m3,m4,m5,m6,m7,m8),modnames=c('egg','fal','let','shit','mint','past','seed','null'),sort = TRUE, c.hat = 1, second.ord = TRUE, nobs = NULL)
AIC_TABLE




##### MODELS FOR BREEDING PROPENSITY DO NOT ALL CONVERGE
#m1<-glmer(breed~Eggplant+(1|Nest_ID), data=NBIcrop,family=binomial)
m2<-glmer(breed~Fallow+(1|Nest_ID), data=NBIcrop,family=binomial)
#m3<-glmer(breed~Letture+(1|Nest_ID), data=NBIcrop,family=binomial)
m4<-glmer(breed~Manure+(1|Nest_ID), data=NBIcrop,family=binomial)
m5<-glmer(breed~Mint+(1|Nest_ID), data=NBIcrop,family=binomial)
m6<-glmer(breed~Pasture+(1|Nest_ID), data=NBIcrop,family=binomial)
m7<-glmer(breed~Seedling+(1|Nest_ID), data=NBIcrop,family=binomial)
m8<-glmer(breed~1+(1|Nest_ID), data=NBIcrop,family=binomial)
#m8<-glmer(breed~Letture+Manure+Mint+(1|Nest_ID), data=breeders_crop,family=poisson)

AIC_TABLE<-aictab(cand.set=list(m2,m4,m5,m6,m7,m8),modnames=c('fal','shit','mint','past','seed','null'),sort = TRUE, c.hat = 1, second.ord = TRUE, nobs = NULL)
AIC_TABLE





#######################################################################################
#
#  DOES FORAGING SUCCESS DIFFER BETWEEN THE DIFFERENT HABITATS?
#
#######################################################################################

### first test habitat differences using the entire data set (2013 and 2014) ###
habmod<-glmer(n_pecks~-1+habitat+(1|ColourRing), data=NBI_hab,family=poisson)
summary(habmod)


### plot output effects
out<-summary(effect('habitat',habmod,se=T))
par(mar=c(5,5,0.5,2.5))
errbar(1:4,out$effect, out$lower, out$upper, xlab="habitat",ylab="n pecks needed for prey capture", axes=F, cex=1.5, cex.lab=1.8, xlim=c(0,5), ylim=c(0,15))
axis(1, at=c(0,1,2,3,4,5), labels=c("",attributes(out$effect)$dimnames$habitat,""), cex.axis=1.5, cex=1.5, cex.lab=1.5)
axis(2, at=seq(0,15,3), labels=T, las=1, cex=1.5, cex.lab=1.5, cex.axis=1.5)





### USE ONLY DATA FROM 2014 FOR DETAILED CROP ANALYSIS ####


NBI_hab<-NBI_hab[NBI_hab$Year==2014,]

m0<-glmer(n_pecks~1+(1|ColourRing), data=NBI_hab,family=poisson)
m2<-glmer(n_pecks~habitat+(1|ColourRing), data=NBI_hab,family=poisson)
m3<-glmer(n_pecks~LANDUSE+(1|ColourRing), data=NBI_hab,family=poisson)
m4<-glmer(n_pecks~CROP+(1|ColourRing), data=NBI_hab,family=poisson)

AIC_TABLE<-aictab(cand.set=list(m2,m3,m4,m0),modnames=c('hab','landuse','crop','null'),sort = TRUE, c.hat = 1, second.ord = TRUE, nobs = NULL)
AIC_TABLE

### plot output effects
out<-summary(effect('CROP',m4,se=T))
par(mar=c(5,5,0,2.5))
errbar(1:6,out$effect, out$lower, out$upper, xlab="crop type",ylab="n pecks needed for prey capture", axes=F, cex=1.5, cex.lab=1.8)
axis(1, at=c(0,1,2,3,4,5,6,7), labels=c("",attributes(out$effect)$dimnames$CROP,""), cex.axis=1.5, cex=1.5, cex.lab=1.5)
axis(2, at=seq(0,13,1), labels=T, las=1, cex=1.5, cex.lab=1.5, cex.axis=1.5)





#######################################################################################
#
#  DOES BEHAVIOUR DIFFER BETWEEN HABITATS AFTER ACCOUNTING FOR TIME?
#
#######################################################################################

### FORAGING BEHAVIOUR ####

m0<-glm(cbind(Number,Total)~month, data = behav[behav$Behaviour=='Foraging',], family=binomial)
m1<-glm(cbind(Number,Total)~CROP, data = behav[behav$Behaviour=='Foraging',], family=binomial)
m2<-glm(cbind(Number,Total)~LANDUSE, data = behav[behav$Behaviour=='Foraging',], family=binomial)
m3<-glm(cbind(Number,Total)~habitat, data = behav[behav$Behaviour=='Foraging',], family=binomial)
m4<-glm(cbind(Number,Total)~hour*month, data = behav[behav$Behaviour=='Foraging',], family=binomial)
m6<-glm(cbind(Number,Total)~month*CROP, data = behav[behav$Behaviour=='Foraging',], family=binomial)
m7<-glm(cbind(Number,Total)~month*LANDUSE, data = behav[behav$Behaviour=='Foraging',], family=binomial)
m8<-glm(cbind(Number,Total)~month*habitat, data = behav[behav$Behaviour=='Foraging',], family=binomial)
m10<-glm(cbind(Number,Total)~1, data = behav[behav$Behaviour=='Foraging',], family=binomial)


AIC_TABLE<-aictab(cand.set=list(m0,m1,m2,m3,m4,m6,m7,m8,m10),modnames=c('month','crop type','landuse','habitat','time full','month*crop type','month*landuse','month*habitat','null'),sort = TRUE, c.hat = 1, second.ord = TRUE, nobs = NULL)
AIC_TABLE


