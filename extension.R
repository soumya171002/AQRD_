library(tidyverse)
library(fixest)
library(broom)
library(panelView)
library(modelsummary)
library(janitor)
library(rbounds)
library(MatchIt)

data <- read.csv("data/Final.csv") |> 
  clean_names() 
fit <- feols(nftp_in_tons ~ treatment+gdp| state + year, data)
tidy(fit)
panelview(
          ntfp ~ treatment,
          data = data,
          index = c("state","year"))
ggsave("PV.png",width = 5, height = 2)


mod <- modelsummary(fit,
              gof_map = c("nobs", "FE: state","FE: year","r.squared","rmse"),
             coef_rename = c("PESA states", "State GDP"),
             stars = TRUE,
             output = "gt")

mod

data |> 
  ggplot(aes(x=year,y = nftp_in_tons,color = state)) + 
  geom_line()+
  geom_vline(xintercept = 2017)+
  labs(
    x = "Non-timber forest produce in tons",
    y = "year",
    title = "Effect of PESA on Non-timber forest produce"
  )
ggsave("DIDI.png", width = 10, height = 4)
#Propensity score matching 
data1 <- read.csv("data/PSM.csv") |> 
  clean_names()

m.out <- matchit(treatment~ gap_in_literacy + st_percent+population,data1, method = NULL, distance = "glm")
summary(m.out)
m.out1 <- matchit(treatment~ gap_in_literacy + st_percent+population,data1, method = "nearest", distance = "glm", replace = TRUE)
summary(m.out1, un = FALSE)
plot(m.out1, type = "jitter", interactive = FALSE)
plot(m.out1, type = "density", interactive = TRUE)
plot(m.out1)
plot(summary(m.out1), xlim = c(0,0.7))
plot(m.out1, interactive = FALSE)
m.out2 <- matchit(treatment~ gap_in_literacy + st_percent+population,data1, method = "nearest", distance = "mahalanobis") 
summary(m.out2)
plot(m.out2, interactive = FALSE)
plot(summary(m.out2), xlim = c(0,0.7))
m.out3 <- matchit(treatment~ gap_in_literacy + st_percent+population,data1, method = "full", distance = "mahalanobis")
plot(summary(m.out3), xlim = c(0,0.7))
m.data <- match.data(m.out1)
head(m.data)

#compare original without matching
fit0 <- lm(forest_cover_loss ~ treatment+gap_in_literacy + st_percent+population,data1)
fit <- lm(forest_cover_loss~ treatment,m.data, weights = weights)
summary(fit)
m.data1 <- match.data(m.out2)
fit1 <- lm(forest_cover_loss~ treatment,m.data1, weights = weights)

summary(fit1)
m.data2 <-match.data(m.out3)
fit2 <- lm(forest_cover_loss~ treatment,m.data2, weights = weights)
modelsummary(list(fit0,fit,fit1,fit2),
             gof_map = c("nobs","rmse","r.squared"),
             coef_rename = c("Intercept","PESA state","literacy gap","ST percent","population"),
             stars = TRUE, output = "gt")


vcf <- read.csv("data/vcf_data_complete.csv")
vcf <- vcf%>%
  filter(year >= 1990)

vcf1 <- vcf %>%
  mutate(
    t = year - 1995,
    ex_ante_med = quantile(cover_1990, 0.5)
  )

# Filter data and add never_treated column in one step
vcf2 <- vcf1 %>%
  filter(cover_1990 > ex_ante_med) %>%
  group_by(cellid) %>%
  mutate(
    never_treated = max(D) == 0
  )


vcffra <- vcf2 |> 
  filter(state == "Chhattisgarh"| state == "Madhya Pradesh"| state == "Orissa"|state == "Jharkhand") |> 
  mutate(Dfra = sch * (year >= 2008))
fit1 <- feols(forest_index~ Dfra| cellid + styear, vcffra)
fit2 <- feols(forest_index ~ D |cellid+ styear, vcffra)


modelsummary(list(fit1,fit2))
vcf |> 
  distinct(state)
vcfrich <- vcf2 |> 
  filter(state == "Maharashtra"|state == "Gujarat" ) |> 
  mutate(Dfra = sch * (year >= 2008))
fit3 <- feols(forest_index ~ D|cellid + styear, vcfrich)
fit4 <- feols(forest_index ~ Dfra|cellid + styear, vcfrich)
modelsummary(list(fit1, fit2, fit3,fit4))

vcfaswing <- vcf2 |> 
  filter(state == "Rajasthan" |state == "Andhra Pradesh"|state == "Himachal Pradesh") |> 
  mutate(Dfra = sch * (year >= 2008))
fit5 <- feols(forest_index ~ D |cellid + styear, vcfaswing)
fit6 <- feols(forest_index ~ Dfra|cellid + styear, vcfaswing)

mod1 <- modelsummary(list(fit1, fit2, fit3,fit4,fit5,fit6),
             gof_map = c("nobs","rmse","FE: cellid","FE: styear","r.squared"),
             coef_rename = c("Foest rights act","PESA"),
             stars = TRUE,
             output = "gt")
  

