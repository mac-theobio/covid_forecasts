library(McMasterPandemic)
library(readr)
library(dplyr)
library(tidyr)
library(shellpipes)

commandEnvironments()

url <- "https://raw.githubusercontent.com/wzmli/COVID19-Canada/master/git_push/clean.Rout.csv"
    
dd <- read_csv(url)
    
all <- (dd
# %>% filter(Province=="ON")
	%>% select(Province,Date,Hospitalization,ICU,Ventilator,deceased,newConfirmations,bestTotal)
	%>% group_by(Province)
	%>% mutate(newDeaths=c(NA,diff(deceased))
## ON hosp includes ICU, our model compartment is just acute care
		, Hospitalization=Hospitalization-ICU
		, newtotal = bestTotal
		, newTests = diff(c(NA,newtotal))
		, newConfirmations = ifelse((newTests == 0) & (newConfirmations == 0), NA, newConfirmations)
	)
	%>% mutate(newConfirmations = ifelse((Province == "BC") & (weekdays(Date)%in%c("Monday","Sunday")), NA, newConfirmations))
	%>% select(-c(deceased,newTests,newtotal))
	%>% pivot_longer(names_to="var",-c(Date,Province))
	%>% setNames(tolower(names(.)))
	%>% ungroup()
)
    
## translate variable names to internally used values
## drop unused variables
keep_vars <- c("H","ICU","death","report","bestTotal")
    
all_sub <- (all
	%>% mutate_at("var",trans_state_vars)
	%>% filter(var %in% keep_vars)
)

print(all_sub)

saveVars(all_sub)
