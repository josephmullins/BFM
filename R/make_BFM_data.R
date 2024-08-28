library(tidyverse)
#setwd("~/Dropbox/PSID_CDS/code/R/")
# - sample: first marriages from 1980 to 1997

M <- read_csv("data/data-psid/Marriage.csv") %>%
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
C <- read.csv("data/data-psid/Childbirth.csv") %>%
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
IndFile <- read.csv("data/data-psid/identifiers-panel.csv") %>%
  mutate(X=NULL)
# read in CPI
cpi <- read.csv("data/CPI-U.csv") %>%
  mutate(CPIU = CPIU / CPIU[54]) #<- normalize to year 2000

L <- read.csv("data/data-psid/earnings-panel.csv") %>%
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

L2 <- read.csv("data/data-psid/T2clean.csv") %>%
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
state <- read.csv("data/data-psid/StateVertical.csv") %>% mutate(X=NULL) #<- load state for each interview
legal_regime <- read.csv("data/StateCodesDivorce.csv") %>% #<- load data on divorce laws
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
educ <- read.csv("data/data-psid/education.csv") %>% mutate(X=NULL)
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


write_csv(sample,"data/MarriageFile.csv")
write_csv(panel,"data/MarriagePanel.csv")

# Now we create a panel of CDS children that belong to this sample of marriages

# ---- Start with test scores and time use
# read in time diary data
# final kid panel:
# test score, tau, phi, div, sep, tsd, tss, AGE

TD <- read.csv("data/data-cds/ActiveTimePanel.csv") %>%
  select(KID,year,tau_m,tau_f) %>%
  rename(YEAR=year)
# make an index of cds kids from the time diary
cds_index <- TD %>%
  select(KID) %>%
  unique()
S <- read.csv("data/data-cds/AssessmentPanel.csv") %>% mutate(X=NULL) %>%
  rename(YEAR=year)


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

write_csv(kid_panel,"data/KidPanelv2.csv")
