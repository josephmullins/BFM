
library(tidyverse)

# make factor shares:
d <- read_csv("../output/factor_shares.csv")

g <- d %>%
  ggplot(aes(x=Age,y=value,ymin=value-1.96*se,ymax=value+1.96*se)) + geom_point() +
  geom_line() + geom_errorbar(width=0.5) + facet_grid(. ~ Input) + theme_minimal() + ylab("Factor Share")
ggsave("../output/figures/factor_shares.png",g,width=5,height=3)

# tfp
d <- read_csv("../output/rel_tfp.csv")

g <- d %>%
  ggplot(aes(x=MarriageQuality,y=TFP,ymin=TFP-1*se,ymax=TFP+1*se)) + geom_errorbar(width=0.2) + 
  geom_line() + geom_point() + geom_hline(yintercept=0,linetype="dashed") + 
  theme_minimal() + ylab("TFP - Married relative to Divorced") + xlab("Marriage Quality")
ggsave("..output/figures/relative_TFP.png",g,width=4,height=3)

# fit of test scores
d <- read_csv("../output/modelfit_testcores.csv")
g <- d %>%
  mutate(group = case_when(Dgroup==1 ~ "Never Divorce",Dgroup==2 ~ "Will Divorce",Dgroup==3 ~ "Divorced")) %>%
  ggplot(aes(x=Age,y=value,linetype=case,color=group)) + geom_line() + 
  theme_minimal() + theme(legend.position="bottom",legend.title = element_blank()) + 
  ylab("Average AP Test Score")
ggsave("../output/figures/modelfit_testscores.png",g)

d <- read_csv("../output/modelfit_testcores_relative.csv")
g <- d %>%
  mutate(group = case_when(Dgroup==1 ~ "Never Divorce",Dgroup==2 ~ "Will Divorce",Dgroup==3 ~ "Divorced")) %>%
  ggplot(aes(x=Age,y=value,linetype=case,color=group)) + geom_line() + 
  theme_minimal() + theme(legend.position="bottom",legend.title = element_blank()) + 
  ylab("Average AP score relative to Never Divorced group")
ggsave("../output/figures/modelfit_testscores_relative.png",g)

d <- read_csv("../output/modelfit_sd.csv")
g <- d %>%
  ggplot(aes(x=Age,y=sd,linetype=case)) + geom_line() + 
  theme_minimal() + theme(legend.position="bottom",legend.title = element_blank()) + 
  ylab("Std. Dev of AP test score")
ggsave("../output/figures/modelfit_sd.png",g)
