## this uses distributional input data to calculate RRs and PAFs

## libraries
library(here)
library(data.table)
library(ggplot2)

## load relevant data
load(here("data/DRB.Rdata"))         # adult BMI distributions
load(here("data/bmirefpop.Rdata"))   # reference BMI distribution
TB <- fread(here("data/TB_burden_age_sex_2024-10-30.csv"))

## Saunders et al linear parameters
## risk per one unit increase in BMI was 14.8% (95%CI: 13.3-16.3)
t <- log(1-0.148) #risk function parameter

## risk ratio calculator
RRfun <- function(k,theta,t) (1-t*bmirefpop$theta)^bmirefpop$k/(1-t*theta)^k

## example dists
curve(dgamma(x,shape=24,scale=0.8),from=10,to=45,n=1e3,col=2,main="Example",xlab="BMI",ylab="")
curve(dgamma(x,shape=bmirefpop$k,scale=bmirefpop$theta),n=1e3,col=1,add=TRUE)

## check
K <- 1e5
bmi1 <- rgamma(K,shape=24,scale=0.7)
bmi0 <- rgamma(K,shape=bmirefpop$k,scale=bmirefpop$theta)
mean(exp(t * (bmi1-30))) / mean(exp(t * (bmi0-30))) #E_1[exp(t*X)] / E_0[exp(t*X)] with a shift for numerics
RRfun(24,0.7,t)                                     #OK

## apply to data
DRB <- DRB[age!="18-19",.(iso3,Year,Sex,age,k,theta)] #restrict
DRB[,RR:=RRfun(k,theta,t)]

## merge against TB estimates

## restrict:
TB <- TB[sex != "a" &
         !age_group %in%  c("all", "0-14", "0-4", "15-24", "15plus", "18plus", "5-14") &
         risk_factor=="all",.(iso3,sex,age=age_group,tb=best)]

## age conversion
akey <- data.table(tbage = c( NA, "25-34", "25-34",  "35-44", "35-44",  "45-54", "45-54",  "55-64", "55-64",
                             "65plus", "65plus", "65plus", "65plus", "65plus"),
                   bmage = c("20-24",  "25-29",  "30-34",  "35-39",  "40-44",  "45-49",  "50-54",  "55-59","60-64",
                             "65-69",  "70-74",  "75-79",  "80-84",  "85plus"))
akey #check


## weight TB evenly over duplicates
TBL <- merge(TB, akey, by.x = "age", by.y = "tbage", allow.cartesian = TRUE)
TBL[,K:=.N,by=.(iso3,sex,age)]
TBL[iso3=="AFG"] #OK
TBL[,tb:=tb/K]   #even weighting...TODO revisit
TBL[,K:=NULL]
TBL[,Sex:=ifelse(sex=="m","Men","Women")]

## merge
DRB <- merge(DRB[age!="20-24"], #NOTE only consider 25+ for now TODO
             TBL[,.(iso3,Sex,age=bmage,tb)],
             by=c("iso3","Sex","age"))

## === aggregation

## over time
RRbyT <- DRB[,.(RR=weighted.mean(RR,tb)),by=Year]
RRbyT[,diff(range(RR))/diff(range(Year))] #0.05 per year ~ 2% per year?

## by age and sex
RRbyAS <- DRB[,.(RR=weighted.mean(RR,tb)),by=.(Sex,age)]


## === plots

## over time
ggplot(RRbyT,aes(Year,RR)) +
  geom_line()+
  theme_linedraw()+
  expand_limits(y=c(0,NA))+
  ylab("Global weighted risk ratio")

ggsave(here("output/RR_year.png"),w=5,h=4)

## by age and sex
ggplot(RRbyAS,aes(age,RR,fill=Sex)) +
  geom_bar(stat="identity",position="dodge")+
  theme_linedraw()+
  expand_limits(y=c(0,NA))+
  ylab("Global weighted risk ratio")+
  theme(legend.position = "top")

ggsave(here("output/RR_age_sex.png"),w=7,h=5)


