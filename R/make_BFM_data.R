library(tidyverse)
setwd("~/Dropbox/PSID_CDS/code/R/")
# - sample: first marriages from 1980 to 1997

M <- read.csv("../../data-main/marriage/Marriage.csv") %>%
  # keep first marriages that have not ended in widowhood
  filter(MH9 == 1, MH12 != 3)

# create a marriage record for mothers/wives
Mw <- M %>%
  filter(MH4 == 2) %>% #<- women only
  mutate(MID = MH2 * 1000 + MH3, FID = MH7 * 1000 + MH8) %>%
  rename(YMAR = MH11, YDIV = MH14, YSEP = MH16) %>%
  filter(YMAR >= 1975, YMAR <= 1997) %>%
  select(MID, FID, YMAR, YDIV, YSEP)

# merge back in symmetrically with husbands/fathers to 
# ensure that we keep the marriage records that are the first for both parties
M <- M %>%
  mutate(FID = MH2*1000 + MH3,MID = MH7*1000 + MH8) %>%
  select(MID,FID) %>%
  merge(Mw)
  

# read in the childbirth data
C <- read.csv("../../data-main/childbirth/Childbirth.csv") %>%
  mutate(KID = CAH10*1000 + CAH11) %>%
  group_by(CAH3,CAH4) %>% 
  filter(sum(CAH2==2 & KID>0)==0) %>% #<- drop parents who have adoptions
  ungroup() %>%
  rename(YBIRTH = CAH15,BNUM = CAH9) %>%
  filter(CAH2==1) #<- keep only childbirth records

# create birth records for mothers
Cm <- C %>%
  filter(CAH5==2,CAH7<9000) %>% #<- women with non-missing birth info
  mutate(MID=CAH3*1000 + CAH4) %>%
  rename(YBIRTH_M=CAH7) %>%
  select(MID,KID,YBIRTH,BNUM,YBIRTH_M)

# create birth records for fathers
Cf <- C %>% 
  filter(CAH5==1,CAH7<9000) %>% #<- men with non-missing birth info
  mutate(FID=CAH3*1000 + CAH4) %>%
  rename(YBIRTH_F=CAH7) %>%
  select(FID,KID,BNUM,YBIRTH_F)

# create the sample of marriages by merging the marriage file
#  with the childbirth records
sample <- Cm %>%
  #<= drop all but the first births of the mother
  filter(BNUM == 1 | BNUM == 99) %>%
  merge(M) %>% #<- merge in the husbands of the first marriage
  # merge in the father of the same kid (will implicitly drop if not
  #  also the first birth for the father or not eventually married)
  merge(Cf) %>%
  # drop if first birth happens before marriage, 
  #  also drop if age when married is greater than 40 or less than 16
  filter(YMAR <= YBIRTH, YMAR - YBIRTH_M <= 40, YMAR - YBIRTH_M >= 16) %>%
  select(-c("KID", "BNUM")) #<- drop child-specific info

# ----- create a panel of hours, earnings, and wages
IndFile <- read.csv("../../data-main/identifiers/identifiers-panel.csv") %>%
  mutate(X=NULL)
# read in CPI
cpi <- read.csv("../../CPI-U.csv") %>%
  mutate(CPIU = CPIU / CPIU[54]) #<- normalize to year 2000

L <- read.csv("../../data-main/labor-market/earnings-panel.csv") %>%
  mutate(X=NULL) %>%
  inner_join(IndFile) %>%
  mutate(hrs = case_when(sn==1 ~ hours_head,sn==2 ~ hours_spouse),earn = case_when(sn==1 ~ earn_head,sn==2 ~ earn_spouse)) %>%
  mutate(wage = case_when(hrs>100 & earn>0 ~ earn / hrs)) %>%
  select(intnum68,pernum,year,hrs,earn,wage) %>%
  mutate(year = year-1) %>% #<- convert to calendar year for these variables
  rename(YEAR = year) %>%
  merge(cpi) %>%
  mutate(wage = wage/CPIU,earn = earn/CPIU) %>%
  select(-CPIU)

L2 <- read.csv("../../data-main/labor-market/T2clean.csv") %>%
  select(MID,FID,earn,wage,year,hrs) %>%
  rename(YEAR = year) %>%
  merge(cpi) %>%
  mutate(wage = wage/CPIU,earn = earn/CPIU) %>%
  select(-CPIU)

L2f <- L2 %>%
    select(-MID)
L2m <- L2 %>%
  select(-FID)


# get earnings for women and men separately
# fathers/husbands
Lf <- L %>%
  mutate(FID = intnum68*1000 + pernum) %>%
  select(FID,YEAR,hrs,earn,wage) %>%
  rbind(L2f) %>%
  inner_join(sample) %>%
  rename(F_hrs  = hrs,F_earn = earn,F_wage = wage)

# mothers/wives
Lm <- L %>%
  mutate(MID = intnum68*1000 + pernum) %>%
  select(MID,YEAR,hrs,earn,wage) %>%
  rbind(L2m) %>%
  inner_join(sample) %>%
  rename(M_hrs  = hrs,M_earn = earn,M_wage = wage) %>%
  arrange(MID,YEAR) %>%
  group_by(MID) %>%
  mutate(exp = cumsum(M_hrs>1000),exp1 = 0.5*cumsum(coalesce(M_hrs,0)>1000) + 0.5*cumsum(coalesce(M_hrs,0)>0),
         exp2 = cumsum(coalesce(M_hrs,0)>1000)) #<- is this good or should we sum over all available measures? I think so?
  

# ------ create a panel of state identifiers
state <- read.csv("../../data-main/state/StateVertical.csv") %>% mutate(X=NULL) #<- load state for each interview
legal_regime <- read.csv("~/Dropbox/RA_Chris/BrownFlinn/Data/StateCodesDivorce.csv") %>% #<- load data on divorce laws
  select(state_str,psid,UD_year,JC_year) %>%
  rename(state = psid)


mar_state <- IndFile %>% #<- link couples to the state of marriage using the closest interview year after getting married
  mutate(MID=intnum68*1000 + pernum) %>%
  inner_join(state) %>%
  inner_join(sample) %>%
  filter(year>=YMAR,state>0,state<=52) %>%
  arrange(MID,year) %>%
  group_by(MID) %>%
  filter(row_number()==1) %>%
  inner_join(legal_regime) %>%
  select(MID,UD_year,JC_year)

# ---- pull in education
educ <- read.csv("../../data-main/education/education.csv") %>% mutate(X=NULL)
Eh <- educ %>%
  rename(FID = ID) %>%
  merge(sample) %>%
  mutate(edH = as.integer(educ>=16)) %>%
  select(FID,edH) %>%
  filter(!is.na(edH))
Ew <- educ %>%
  rename(MID = ID) %>%
  merge(sample) %>%
  mutate(edW = as.integer(educ>=16)) %>%
  select(MID,edW) %>%
  filter(!is.na(edW))

# put together the sample of parents with horizontal data
sample <- sample %>%
  inner_join(mar_state) %>%
  inner_join(Eh) %>%
  inner_join(Ew)


# ------ create a panel for each marriage:
yrs = data.frame(x=1,YEAR = seq(1975,2010))

panel <- sample %>%
  mutate(x=1) %>%
  inner_join(yrs,relationship = "many-to-many") %>%
  filter(YEAR>=YMAR) %>%
  select(-x) %>%
  left_join(Lf) %>%
  left_join(Lm) %>%
  mutate(F_AGE = YEAR - YBIRTH_F,M_AGE = YEAR - YBIRTH_M,SEP = YEAR>=YSEP,DIV = YEAR>=YDIV,
         TSD = YEAR-YDIV,TSS = YEAR-YSEP,L = 1 + as.integer(YEAR>=UD_year)) %>%
  filter(M_AGE<60) %>%
  arrange(MID,YEAR)

# final calculations for the horizontal file
sample <- panel %>%
  group_by(MID,FID) %>%
  summarize(tlength = n(),exp0 = coalesce(exp2[1],0)) %>%
  merge(sample) %>%
  arrange(MID)

  

#write.csv(sample,"~/Dropbox/Research Projects/BrownFlinnMullins/Data/MarriageFile.csv")
write_csv(sample,"~/BFM/data/MarriageFile.csv")
#write.csv(panel,"~/Dropbox/Research Projects/BrownFlinnMullins/Data/MarriagePanel.csv")
write_csv(panel,"~/BFM/data/MarriagePanel.csv")


# Now let's create a panel of CDS children that belong the sample of marriages


# ---- Now let's make the panel of test scores and time use
# read in time diary data
# final kid panel:
# test score, tau, phi, div, sep, tsd, tss, AGE

TD <- read.csv("../../data-cds/time-diary/ActiveTimePanel.csv") %>%
  select(KID,year,tau_m,tau_f) %>%
  rename(YEAR=year)
# make an index of cds kids from the time diary
cds_index <- TD %>%
  select(KID) %>%
  unique()
S <- read.csv("../../data-cds/assessments/AssessmentPanel.csv") %>% mutate(X=NULL) %>%
  rename(YEAR=year)

# old code which was incorrectly getting age
# cds_index <- sample %>%
#   select(MID,FID) %>%
#   merge(Cm) %>%
#   select(-BNUM) %>%
#   merge(Cf) %>%
#   merge(cds_index) %>%
#   select(KID,MID,FID)
# 
# kid_panel <- cds_index %>%
#   merge(panel) %>%
#   merge(S) %>%
#   left_join(TD) %>%
#   mutate(AGE=YEAR-YBIRTH,PHI_M = tau_m / (112 - M_hrs/52), PHI_F = tau_f / (112 - F_hrs/52)) %>%
#   mutate(dgroup = case_when(YDIV>9000 ~ 1,YDIV<=YEAR ~ 2,YDIV>YEAR ~ 3),
#          sgroup = case_when(YDIV>9000 ~ 1,YSEP<=YEAR ~ 2,YSEP>YEAR ~ 3)) %>%
#   arrange(MID,KID,YEAR)

cds_index <- sample %>%
  select(MID,FID) %>%
  merge(Cm) %>%
  select(-BNUM) %>%
  merge(Cf) %>%
  merge(cds_index) %>%
  select(KID,MID,FID,YBIRTH)

kid_panel <- panel %>%
  select(-YBIRTH) %>%
  merge(cds_index) %>%
  merge(S) %>%
  left_join(TD) %>%
  mutate(AGE=YEAR-YBIRTH,PHI_M = tau_m / (112 - M_hrs/52), PHI_F = tau_f / (112 - F_hrs/52)) %>%
  mutate(dgroup = case_when(YDIV>9000 ~ 1,YDIV<=YEAR ~ 2,YDIV>YEAR ~ 3),
         sgroup = case_when(YDIV>9000 ~ 1,YSEP<=YEAR ~ 2,YSEP>YEAR ~ 3)) %>%
  arrange(MID,KID,YEAR)

write_csv(kid_panel,"~/BFM/data/KidPanelv2.csv")

break
TD %>%
  select(KID) %>%
  unique() %>%
  merge(Cm) %>%
  group_by(YBIRTH) %>%
  summarize(n())

# these two regressions show us some interesting stuff
kid_panel %>%
  mutate(div_exp = pmin(pmax(0,TSD),10)) %>%
  lm(AP_raw ~ DIV + div_exp + poly(AGE,2),data=.) %>%
  summary()

# adding the "will separate" variable shows that some of this is due to selection
kid_panel %>%
  mutate(WS = YSEP<2010,div_exp = pmin(pmax(0,TSS),10)) %>%
  lm(AP_std ~ WS + SEP + div_exp + poly(AGE,2),data=.) %>%
  summary()

# a contrast here
kid_panel %>%
  mutate(WD=YSEP<2010,TS1 = TSS>-5 & TSS<=0,TS2 = TSS>0 & TSS<=5,TS3 = TSS>5) %>%
  lm(AP_raw ~ WD + TS1 + TS2 + TS3 + poly(AGE,2),data=.) %>%
  summary()

# this is interesting, shows positive impacts in a fixed effects regression
library(fixest)
# similar result to the cross-section regression above
kid_panel %>%
  group_by(AGE) %>%
  mutate(AP_std2 = (AP_raw - mean(AP_raw,na.rm = TRUE))/sd(AP_raw,na.rm=TRUE)) %>%
  mutate(div_exp = pmin(pmax(0,TSS),10)) %>%
  feols(AP_std2 ~ SEP + div_exp | AGE + MID, data=.) %>%
  summary()

kid_panel %>%
  group_by(AGE) %>%
  mutate(AP_std2 = (AP_raw - mean(AP_raw,na.rm = TRUE))/sd(AP_raw,na.rm=TRUE)) %>%
  mutate(div_exp = pmin(pmax(0,TSS),15),div_exp2 = pmax(0,TSS)) %>%
  feols(AP_std2 ~ SEP + div_exp2 | AGE + MID, data=.) %>%
  summary()

S2 <- cds_index %>%
  merge(panel) %>%
  merge(S) %>%
  mutate(tsd = YEAR-YDIV,tss = YEAR-YSEP)
  
S2 %>%
  mutate(AGE = YEAR-YBIRTH) %>%
  filter(tsd>=-8,tsd<=8,!is.na(LW_std)) %>%
  group_by(AGE) %>%
  mutate(LW_raw = LW_raw - mean(LW_raw)) %>%
  group_by(tsd) %>%
  summarize(k =mean(LW_raw),se = sd(LW_raw)/n()) %>%
  ggplot(aes(x=tsd,y=k,ymin=k-1.96*se,ymax=k+1.96*se)) + geom_point() + geom_errorbar()

S2 %>%
  mutate(AGE = YEAR-YBIRTH) %>%
  filter(tss>=-5,tss<=5,!is.na(LW_std)) %>%
  mutate(LW_raw = AP_raw - mean(AP_raw)) %>%
  group_by(tss) %>%
  summarize(k =mean(LW_raw),se = sd(LW_raw)/n()) %>%
  ggplot(aes(x=tss,y=k,ymin=k-1.96*se,ymax=k+1.96*se)) + geom_point() + geom_errorbar()


S2 %>%
  filter(YEAR<=2002) %>%
  mutate(AGE = YEAR-YBIRTH) %>%
  arrange(KID,YEAR) %>%
  group_by(KID) %>%
  filter(n()==2) %>%
  summarize(tss=tss[1],tsd=tsd[1],dLW = LW_raw[2]-LW_raw[1]) %>%
  filter(tss>=-5,tss<=5,!is.na(dLW))
  

# this so far is the clearest looking graph we can seem to get.
#   NOTE: the pattern looks a bit differnet for LW compared to AP
kid_panel %>%
  filter(!is.na(AP_raw)) %>%
  mutate(status = case_when(YEAR>=YSEP ~ "D",YSEP<=2012 ~ "WD",YSEP>2012 ~ "ND")) %>%
  mutate(AGE = round(AGE/4)) %>%
  group_by(AGE,status) %>%
  summarize(k = mean(AP_raw),se = sd(AP_raw)/sqrt(n())) %>%
  ggplot(aes(x=AGE,y=k,color=status,ymin=k-1.96*se,ymax=k+1.96*se)) + geom_point() + geom_line()

# try running a regression with final scores and years divorced?

# ---- Husbands Income
mod <- panel %>%
  mutate(AGE_F = YEAR-YBIRTH_F) %>%
  filter(AGE_F>18) %>%
  filter(!is.na(F_wage),F_wage>0,edH==1) %>%
  lm(log(F_wage) ~ poly(AGE_F,2),data=.) #%>%
  summary()

panel %>%
  mutate(AGE_F = YEAR-YBIRTH_F) %>%
  filter(AGE_F>18) %>%
  filter(!is.na(F_wage),F_wage>0,edH==1) %>%
  group_by(AGE_F) %>%
  summarize(m = mean(log(F_wage))) %>%
  ggplot(aes(x=AGE_F,y=m)) + geom_point() + geom_line()

# ---- regressions without control, seem to suggest that the second measure of experience is preferred
panel %>%
  mutate(AGE_M = YEAR-YBIRTH_M) %>%
  filter(AGE_M>18) %>%
  filter(!is.na(M_wage),M_wage>0,edW==1) %>%
  lm(log(M_wage) ~ poly(AGE_M,2) + poly(exp1,2),data=.) %>%
  summary()

panel %>%
  mutate(AGE_M = YEAR-YBIRTH_M) %>%
  filter(!is.na(M_wage),M_wage>0,edW==1) %>%
  lm(log(M_wage) ~ poly(AGE_M,2) + poly(exp2,2),data=.) %>%
  summary()

panel %>%
  mutate(AGE_M = YEAR-YBIRTH_M) %>%
  filter(!is.na(M_wage),M_wage>0,edW==0) %>%
  lm(log(M_wage) ~ poly(AGE_M,2) + poly(exp1,2),data=.) %>%
  summary()

panel %>%
  mutate(AGE_M = YEAR-YBIRTH_M) %>%
  filter(!is.na(M_wage),M_wage>0,edW==0) %>%
  lm(log(M_wage) ~ poly(AGE_M,2) + poly(exp2,2),data=.) %>%
  summary()

panel %>%
  mutate(AGE_M = YEAR-YBIRTH_M) %>%
  filter(!is.na(M_wage),M_wage>0,edW==0) %>%
  lm(log(M_wage) ~ poly(AGE_M,2) + exp2,data=.) %>%
  summary()

panel %>%
  mutate(AGE_M = YEAR-YBIRTH_M) %>%
  filter(!is.na(M_wage),M_wage>0,edW==0) %>%
  lm(log(M_wage) ~ AGE_M + exp2,data=.) %>%
  summary()


# ---- now let's estimate a propensity
d <- panel %>%
  mutate(AGE_M = YEAR-YBIRTH_M) %>%
  filter(!is.na(M_hrs),!is.na(F_earn),!is.na(exp2),!is.na(edW))

pmod <- d %>%
  lm(M_hrs>0 ~ edW*poly(AGE_M)*poly(exp2)*poly(F_earn,2),data=.)

d$pwork <- predict(pmod)

d %>%
  filter(edW==1) %>%
  lm(log(M_wage) ~ poly(exp2,2) + poly(pwork,2),data=.) %>%
  summary()

d %>%
  filter(edW==1) %>%
  lm(log(M_wage) ~ poly(exp2,2),data=.) %>%
  summary()

# this picture suggests selection might be a concern for college but not for non-college groups
#   - it seems to be an issue only for the intercept term, however.
# Chris says don't include it. I think we should.
# We need to drop observations with missing ages from the dataset.

d %>%
  filter(M_wage>0) %>%
  mutate(pgroup = pwork>0.85) %>%
  group_by(exp2,pgroup,edW) %>%
  summarize(logw = mean(log(M_wage))) %>%
  ggplot(aes(x=exp2,y=logw,color=pgroup)) + geom_point() + geom_line() + facet_grid(. ~ edW)


# this picture shows us that the experience variable looks reasonable
# BUT: we need to think about the measurement error
d %>%
  group_by(edW,AGE_M) %>%
  filter(AGE_M>0) %>%
  summarize(exp = mean(exp2,na.rm = TRUE)) %>%
  ggplot(aes(x=AGE_M,y=exp,color=as.factor(edW))) + geom_point()
