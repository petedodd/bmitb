## this uses distributional input data to calculate RRs and PAFs

## libraries
library(here)
library(data.table)
library(ggplot2)
library(ggrepel)

## load relevant data
load(here("data/whokey.Rdata"))      # WHO region to iso3 key
load(here("data/DRB.Rdata"))         # adult BMI distributions
load(here("data/bmirefpop2.Rdata"))   # reference BMI distribution
bmirefpop <- bmirefpop2               #try new version
TB <- fread(here("data/TB_burden_age_sex_2024-10-30.csv"))
TBT <- fread(here("rawdata/TB_burden_countries_2024-10-30.csv"))

## Saunders et al linear parameters
## risk per one unit increase in BMI was 14.8% (95%CI: 13.3-16.3)
t <- -log(1.148) #risk function parameter
exp(-t)          #risk increase with 1 unit decreast


## risk ratio calculator
RRfun <- function(k,theta,t) (1-t*bmirefpop$theta)^bmirefpop$k/(1-t*theta)^k


## check
K <- 1e5
bmi1 <- rgamma(K,shape=24,scale=0.7)
bmi0 <- rgamma(K,shape=bmirefpop$k,scale=bmirefpop$theta)
mean(exp(t * (bmi1-30))) / mean(exp(t * (bmi0-30))) #E_1[exp(t*X)] / E_0[exp(t*X)] with a shift for numerics
RRfun(24,0.7,t)                                     #OK

## bilinear version
## see: https://search.r-project.org/CRAN/refmans/expint/html/gammainc.html
## Γ(a,x)=Γ(a)(1−P(a,x))
## https://search.r-project.org/R/refmans/stats/html/GammaDist.html
t1 <- t2 <- t
RRfunBL0 <- function(k,theta,t1,t2){
  x0 <- 25; a1 <- -t1; a2 <- -t2;
  ## stuff x Γ(k,x0*(a2+1/theta))/Γ(k) = stuff x (1-P(k,x0*(a2+1/theta))) = stuff x P(k,x0*(a2+1/theta),lower=FALSE)
  lans2 <- a2*x0 + pgamma(k,x0*(a2+1/theta),log=TRUE,lower=TRUE) - k * log(1+a2*theta)
  lans1 <- a1*x0 + pgamma(k,x0*(a1+1/theta),log=TRUE,lower=FALSE) - k * log(1+a1*theta)
  exp(lans1) + exp(lans2)
}
## NOTE works in tests

RRfunBL <- function(k,theta,t1,t2) RRfunBL0(k,theta,t1,t2) / RRfunBL0(bmirefpop$k,bmirefpop$theta,t1,t2)


## check
BL <- function(x,t1,t2){
  ans <- (x-25)
  less <- ans<0
  ans[less] <- t1*ans[less]
  ans[!less] <- t2*ans[!less]
  ans
}
xx <- seq(from=10,40,by=0.1)
plot(xx,BL(xx,t1,t2),type="l")
plot(xx,BL(xx,t1,t2-0.1),type="l")

mean(exp(BL(bmi1,t1,t1))) / mean(exp(BL(bmi0,t1,t1)))
mean(exp(t1 * (bmi1-30))) / mean(exp(t1 * (bmi0-30))) #E_1[exp(t*X)] / E_0[exp(t*X)] with a shift for numerics
RRfunBL0(24,0.7,t1,t1) / RRfunBL0(bmirefpop$k,bmirefpop$theta,t1,t1)
RRfunBL(24,0.7,t1,t1)

## 18.0% (95%CI: 16.4-19.6) for BMI<25.0kg/m2 and 6.9% (95%CI: 4.6-9.2) for BMI>=25.0kg/m2 in
t1 <- -log(1.18) #risk function parameter
exp(-t1)          #risk increase with 1 unit decreast
t2 <- -log(1.069) #risk function parameter
exp(-t2)          #risk increase with 1 unit decreast

mean(exp(BL(bmi1,t1,t2))) / mean(exp(BL(bmi0,t1,t2)))
RRfunBL(24,0.7,t1,t2) ## OK



## example dists: exaggerated
png(here("output/eg_dist.png"))

curve(dgamma(x,shape=24,scale=0.8),from=10,to=45,n=1e3,col=2,main="Example",xlab="BMI",ylab="")
curve(dgamma(x,shape=bmirefpop$k,scale=bmirefpop$theta),n=1e3,col=1,add=TRUE)

dev.off()

## distribution of where TB comes from
wts <- exp(t*(bmi0-30))
bmi0tb <- bmi0[sample(1:K,size=K,replace=TRUE,prob = wts)]
bmiz <- data.table(population=c(rep("whole population",K),rep("TB",K)),
                   BMI=c(bmi0,bmi0tb))

ggplot(bmiz,aes(x=BMI,y=after_stat(density),fill=population))+geom_histogram(alpha=0.5) +
  theme_classic() + theme(legend.position = "top",legend.title = element_blank()) + ggpubr::grids()

ggsave(file=here("output/eg_tb_vs_pop.png"),w=6,h=5)

## correct magnitude of change?
mean(bmi1)
mean(bmi0)
exp(t*(mean(bmi1)-mean(bmi0)))
RRfun(24,0.7,t) #about as good as one might expect

## apply to data
DRB <- DRB[age!="18-19",.(iso3,Year,Sex,age,k,theta)] #restrict
DRB[,RR:=RRfun(k,theta,t)]

## merge against TB estimates

## restrict:
TB <- TB[sex != "a" &
         !age_group %in%  c("all", "0-14", "0-4", "15-24", "15plus", "18plus", "5-14") &
         risk_factor=="all",.(iso3,sex,age=age_group,tb=best)]

## age conversion
akey <- data.table(tbage = c( NA, "25-34", "25-34",  "35-44", "35-44",  "45-54",
                             "45-54",  "55-64", "55-64",
                             "65plus", "65plus", "65plus", "65plus", "65plus"),
                   bmage = c("20-24",  "25-29",  "30-34",  "35-39",  "40-44",
                             "45-49",  "50-54",  "55-59","60-64",
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

DRB <- merge(DRB,whokey,by="iso3") #WHO regions

## === sense check
DRB[,mnbmi:=k*theta] #for each country year
mnbmi0 <- bmirefpop[,k*theta]
DRB[,range(mnbmi)] #TODO some absolute crazies
DRB[,range(RR)]    #TODO ditto

ggplot(DRB[RR>0.5 & RR<10],aes(mnbmi,RR,col=g_whoregion))+geom_point(shape=1)+
  facet_wrap(~g_whoregion)+theme_linedraw()+theme(legend.position="none")+
  xlab("Mean BMI in each group") + ylab("Associated relative risk")

ggsave(here("output/RR_check1.png"),w=15,h=10)


## === aggregation

## over time
RRbyT <- DRB[,.(RR=weighted.mean(RR,tb)),by=Year]
RRbyT[,diff(range(RR))/diff(range(Year))] #0.05 per year ~ 2% per year?

## over time & by country
RRbyTC <- DRB[,.(RR=weighted.mean(RR,tb)),by=.(Year,iso3)]
RRbyTC <- merge(RRbyTC,TBT[,.(iso3,Year=year,e_inc_100k,g_whoregion)],by=c("Year","iso3"))
RRbyTC[,hi:=mean(e_inc_100k)>100,by=iso3]
RRbyTC[,mn:=mean(e_inc_100k),by=iso3]
RRbyTC[Year==2022,lbl:=iso3]

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

## over time & by country

RRbyTC[iso3=="BDI"]

## for countries with mean per capita over ctp
ctp <- 15
ggplot(RRbyTC[mn>ctp],aes(RR,e_inc_100k,group=iso3,col=g_whoregion,label=lbl)) +
  geom_line()+
  geom_point(data=RRbyTC[mn>ctp][!is.na(lbl)])+
  geom_text_repel(show.legend=FALSE)+
  theme_linedraw()+
  scale_y_log10() + scale_x_log10()+
  ylab("Estimated TB incidence per 100ky  (log scale)")+
  facet_wrap(~g_whoregion)+
  theme(legend.position="none")+
  xlab("Risk-ratio from BMI distribution (log scale)")


ggsave(here("output/RR_year_country.png"),w=15,h=10)

## by age and sex
ggplot(RRbyAS,aes(age,RR,fill=Sex)) +
  geom_bar(stat="identity",position="dodge")+
  theme_linedraw()+
  expand_limits(y=c(0,NA))+
  ylab("Global weighted risk ratio")+
  theme(legend.position = "top")

ggsave(here("output/RR_age_sex.png"),w=7,h=5)


