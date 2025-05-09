
## script to process NCD RisC data and fit distribution parameters

## raw data ncdris data is from: https://ncdrisc.org/data-downloads-adiposity.html
## stored in rawdata/ are the following files:
## NCD_RisC_Lancet_2024_BMI_child_adolescent_country.csv
## NCD_RisC_Lancet_2024_BMI_female_age_specific_country.csv
## NCD_RisC_Lancet_2024_BMI_male_age_specific_country.csv

## libraries
library(here)
library(data.table)
library(readxl)

## === WHO refs for children
## want to know: BMIs for z scores -2,-1,+1,+2 in ref pop for sex and age
## from: https://www.who.int/tools/growth-reference-data-for-5to19-years/indicators/bmi-for-age
fn <- here("data/bmi-boys-z-who-2007-exp.xlsx")
WB <- read_excel(fn)
fn <- here("data/bmi-girls-z-who-2007-exp.xlsx")
WG <- read_excel(fn)
WB <- as.data.table(WB)
WG <- as.data.table(WG)
WG[,Sex := "Girls"]
WB[,Sex := "Boys"]
BG <- rbind(WB, WG) #combined version


## === read in raw data
## --- adults
fn <- here("rawdata/NCD_RisC_Lancet_2024_BMI_female_age_specific_country.csv") #female
F <- fread(fn) #NOTE not in repo as large
fn <- here("rawdata/NCD_RisC_Lancet_2024_BMI_male_age_specific_country.csv") #males
M <- fread(fn) #NOTE not in repo as large
B <- rbind(M, F) #both

## --- adolescents
fn <- here("rawdata/NCD_RisC_Lancet_2024_BMI_child_adolescent_country.csv")
A <- fread(fn) #NOTE not in repo as large

## compare variables
setdiff(names(A), names(M))
setdiff(names(M), names(A))

## compare ages
A[,unique(`Age group`)]
B[,unique(`Age group`)]

## compare years
A[,unique(Year)]
B[,unique(Year)] #both 1990 - 2022


A[ISO=="AFG" & Year==2022]

## available for all
A[,table(Year,is.na(`Prevalence of BMI < -2SD (thinness) lower 95% uncertainty interval`))]

## restrict data
DR <- A[Year==2022] #so the distribution info is available before 2017
DR <- DR[`Age group`>=5]
(keep <- names(DR)[c(1,2,3,5,6,9,18,21,24)])
DR <- DR[,..keep]

## new names
nnmz <- c("Year","Sex","age","iso3","lt.n2","gt.p2","n2.n1","n1.p2","p1.p2")
cbind(names(DR),nnmz) #CHECK
names(DR) <- nnmz     #rename
tgt <- nnmz[-c(1:4)]  #targets


## === fitting to gamma distributions
## --- adolescents

## fit to gamma
gammaStats <- function(k,theta,
                       b2112 #BMI refs for +2,+1,-1,-2 SD
                       ){
  sts <- c(pgamma(b2112[4],shape=k,scale=theta),   #low tails:  <= -2Z
           1-pgamma(b2112[1],shape=k,scale=theta), #high tails: >= +2Z
    pgamma(b2112[3],shape=k,scale=theta)-pgamma(b2112[4],shape=k,scale=theta),    #intervals: -2Z to -1Z
    pgamma(b2112[2],shape=k,scale=theta)-pgamma(b2112[3],shape=k,scale=theta),    #intervals: -1Z to +1Z
    pgamma(b2112[1],shape=k,scale=theta)-pgamma(b2112[2],shape=k,scale=theta)     #intervals: +1Z to +2Z
    )
  names(sts) <- tgt
  sts
}

## test
ba5 <- unlist(WB[Month==5*12+6,.(SD2,SD1,SD1neg,SD2neg)])
gammaStats(2,9,ba5)

## loss function
gamErr <- function(x,tgt,ref){
  x <- exp(x)
  res <- gammaStats(x[1],x[2],ref)
  mean((res/(tgt+1e-10)-1)^2) #~SSE relative, equal weight
}

## test
TGT <- unlist(DR[1,..tgt])
gammaStats(1,1,ba5)
gamErr(c(0,0),TGT,ba5)
## optimize
(out <- optim(par=c(0,0),fn=function(x) gamErr(x,TGT,ba5)))
exp(out$par)
gammaStats(exp(out$par[1]),exp(out$par[2]),ba5)
TGT #NOTE not great on item 2
## check
curve(dgamma(x,shape=exp(out$par[1]),scale=exp(out$par[2])),from=10,to=30,n=1e3)

## TODO explore other distributions (though losing the convenient analytical result)

##  ---apply this fitting across the data:
## common age to merge with WHO REFs
DR[,Month:=12*age] #low-point (otherwise miss 19-20)
DR[Month==60,Month:=61]    #otherwise missing from BG
DR[,range(Month)]; BG[,range(Month)]
akey <- unique(DR[,.(age,Month)])

## merge:
DR <- merge(DR,BG,by=c("Sex","Month"),all.x = TRUE)
unique(DR[,.(age,Month)]) #CHECK

## Loop above work flow
DRA <- DR[,{
  print(iso3);print(Sex);print(Month);
  TGT <- c(lt.n2,gt.p2,n2.n1,n1.p2,p1.p2)       #NCD RisC targets
  saref <- c(SD2, SD1, SD1neg, SD2neg);             #BMIs in WHO reference data
  out <- optim(par=c(0,0),fn=function(x) gamErr(x,TGT,saref));
  if(abs(out$convergence)>0) print("*** has not converged! ***")
  list(k=exp(out$par[1]),theta=exp(out$par[2]),cvgc=out$convergence)
},by=.(iso3,Month,Sex,age)]

DRA[,table(cvgc)] #not many

fn <- here("data/DRA.Rdata")
save(DRA,file=fn)

## TODO fix problems for 22

## checks
DRA[,mn:=k*theta]
DRA[,summary(mn)]
## DRA[,qplot(mn)]

## ## convergence problem for:
## ## Angola Girls 61
## ## SLE Girls 156
## ## Congo Girls
## ## RUS girls 204
## ## THA G 204
## ## Congo G 216
## ## Libera girls 228


## DRA[mn>25] #PNG boys 120

## ## inspect problems
## bad <- paste(c('Angola','Sierra Leone','Russian Federation','Thailand','Congo','Libera'),
##              rep('Girls',7),
##              c(61,156,204,204,216,228),
##              sep='_')
## bad <- c(bad,'PNG_Boys_120')

## DRA[,id:=paste(Country,Sex,Month,sep='_')]
## DR[,id:=paste(Country,Sex,Month,sep='_')]

## DRA[id %in% bad]
## DR[id %in% bad]

## ## replace these units with age/sex averages
## DRAM <- DRA[,.(km=median(k),thetam=median(theta)),by=.(Sex,Month)]
## DRAM[,km*thetam,by=.(Sex,Month)] #check
## DRA <- merge(DRA,DRAM,by=c('Sex','Month'))
## DRA[id %in% bad,c('k','theta'):=.(km,thetam)] #replace

## ## check
## DRA[,mn:=k*theta]
## DRA[mn>25] #PNG boys 120
## DRA[,summary(mn)] #OK

## ## drop extra vars
## DRA[,c('mn','Month','km','thetam','id'):=NULL]

## ## save
## save(DRA,file=here('rawdata/DRA.Rdata'))

## ======= adult data

## available for all
B[,table(Year,is.na(`Prevalence of BMI >=40 kg/m² (morbid obesity)`))]

## restrict
BR <- B## [Year==2022] #so the distribution info is available before 2017
keep <- names(BR)[c(1,2,4,5,
                    6,9,18,21,24,27,30,33)]
BR <- BR[,..keep]

## new names
nnmz <- c("Year","Sex","iso3","age","l185","g30","l20.g185","l25.g20",
          "l30.g25","l35.g30","l40.g35","g40")
cbind(names(BR),nnmz) #CHECK
names(BR) <- nnmz     #rename
tgt <- nnmz[-c(1:4)]  #targets


## fit to gamma
AgammaStats <- function(k,theta
                       ){
  sts <- c(pgamma(18.5,shape=k,scale=theta),   #low tails:  <= 18.5
    1-pgamma(30,shape=k,scale=theta), #high tails: >= 30
    pgamma(20,shape=k,scale=theta)-pgamma(18.5,shape=k,scale=theta),    #intervals
    pgamma(25,shape=k,scale=theta)-pgamma(20,shape=k,scale=theta),    #intervals
    pgamma(30,shape=k,scale=theta)-pgamma(25,shape=k,scale=theta),     #intervals
    pgamma(35,shape=k,scale=theta)-pgamma(30,shape=k,scale=theta),     #intervals
    pgamma(40,shape=k,scale=theta)-pgamma(35,shape=k,scale=theta),     #intervals
    1-pgamma(40,shape=k,scale=theta) #high tails: >= 40
    )
  names(sts) <- tgt
  sts
}

## test
AgammaStats(2,9)


## loss function
AgamErr <- function(x,tgt){
  x <- exp(x)
  res <- AgammaStats(x[1],x[2])
  mean((res/(tgt+1e-10)-1)^2) #~SSE relative, equal weight
}

## test
TGT <- unlist(BR[1,..tgt])
AgamErr(c(0,0),TGT)
## optimize
(out <- optim(par=c(0,0),fn=function(x) AgamErr(x,TGT)))
exp(out$par)
AgammaStats(exp(out$par[1]),exp(out$par[2]))
TGT
## check
curve(dgamma(x,shape=exp(out$par[1]),scale=exp(out$par[2])),from=10,to=40,n=1e3)


##  ---apply this fitting across the data:

## Loop above work flow
DRB <- BR[,{
  cat(iso3,", ",Sex,", ",age,", ",Year,"\n");
  TGT <- c(l185,         g30,    l20.g185,     l25.g20,     l30.g25,
           l35.g30,     l40.g35,         g40 )
  out <- optim(par=c(0,0),fn=function(x) AgamErr(x,TGT));
  if(abs(out$convergence)>0) print('*** has not converged! ***')
  list(k=exp(out$par[1]),theta=exp(out$par[2]),cvgc=out$convergence)
},by=.(iso3,Year,Sex,age)]

DRB[,table(cvgc)] #not many

fn <- here("data/DRB.Rdata")
save(DRB,file=fn)

## TODO fix 4

DRB[,unique(age)]
DRB[,unique(iso3)]
