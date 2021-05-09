library(McMasterPandemic)
library(dplyr)
library(readr)
library(gsheet)
library(tidyr)
library(zoo)
library(ggplot2)
library(shellpipes)

commandEnvironments()

trim_date <- as.Date("2020-09-15")


## using latest date from clean.RData
if(!exists("end_date")){
end_date <- as.Date(max(all_sub$date))
}

all_inputs <- csvRead()

immunity_lag <- 14
vac_effect <- 0.6
data_scale <- 1

vacdat <- data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/1PjkemMdFSZgA-M8Esr6rbNjHiyfcXcBxPeMjselJIso/edit#gid=0"))


info <- all_inputs %>% filter(province == "ON")
print(info)
prov = info[["province"]]
params <- fix_pars(read_params(paste0(prov,".csv"))
	, target = c(R0 = 1.3 , Gbar=6)
)
params[["obs_disp"]] <- 40
params[["vacc"]] <- 1e-10
   
   # Retrieve break dates:
bd <- as.Date(unlist(strsplit(info[["break_dates"]],split = ";")))
n.bd = length(bd)

## cleaning vaccination data
clean_vac <- (vacdat
	%>% transmute(date = as.Date(data...date) + immunity_lag
		, daily_vac = diff(c(0,data...total_vaccinations))
		)
	%>% filter(date %in% bd)
	%>% mutate(daily_vac = daily_vac/info$population/1e-10
		, daily_vac = ifelse(is.na(daily_vac) | daily_vac==0,1,vac_effect*daily_vac*data_scale)
		, daily_vac = ifelse(date < as.Date("2021-01-01"), 1, daily_vac)
	)
)

print(clean_vac)
	
summary(params)
params[["N"]] <- info[["population"]]
		
lgf <- function(x){log(x/(1-x))}
opt_pars <- list(params = c(log_beta0= log(params[["beta0"]])
	# ,logit_mu = lgf(params[["mu"]])
	# , logit_phi1 = lgf(params[["phi1"]])
	)
	, rel_beta0 = c(0.8,0.9,0.7,0.6,1,1,0.8)
#	, rel_beta0 = rep(1,n.bd)
	, rel_vacc = clean_vac$daily_vac
	)
	
## Giving vaccination rates a narrow prior 
## The fix_pars option does not work yet
## The narrow prior is a temp hack

priors= list(~dnorm(rel_vacc[1], mean=clean_vac$daily_vac[1],sd=0.001)
	, ~dnorm(rel_vacc[3], mean=clean_vac$daily_vac[3],sd=0.001)
	, ~dnorm(rel_vacc[4], mean=clean_vac$daily_vac[4],sd=0.001)
	, ~dnorm(rel_vacc[5], mean=clean_vac$daily_vac[5],sd=0.001)
	, ~dnorm(rel_vacc[6], mean=clean_vac$daily_vac[6],sd=0.001)
	, ~dnorm(rel_vacc[7], mean=clean_vac$daily_vac[7],sd=0.001)
)
	

# Subset the data to the requested 
# province and variables:
province_dat <- (all_sub
	%>% group_by(var)
	%>% filter(province == info[["province"]])
	%>% filter(grepl(var, info[["vars"]]))
)

# Define time window
start_date <- trim_date
start_date_offset <- 60
date_vec <- as.Date(start_date:end_date)
date_df <- data.frame(date = rep(date_vec,length(unique(province_dat[["var"]]))), 
	var  = rep(unlist(strsplit(info[["vars"]],"/")),each=length(date_vec))
)
		
dat <- (left_join(date_df, province_dat)
	%>% mutate(value = ifelse(value == 0, NA, value))
)
	
## ==== Model calibration ====
	
fitdat <- dat
	
res <- calibrate(base_params  = params
	, priors = priors
	, debug_plot = FALSE
	, debug      = FALSE
	, data       = fitdat
	, opt_pars   = opt_pars
	, DE_cores   = info$DE_cores
	, sim_args   = list(ndt = 2)
	, time_args  = list(break_dates = bd)
	, start_date_offset = start_date_offset
	, use_DEoptim = FALSE
)
	
	
## parameters getting passed straight through 
## from the data frame to calibrate_comb()

res_list <- list(fit=res, inputs=info, trimdat = fitdat, fulldat=dat)

saveRDS(object=res_list, file=paste0("./cachestuff/",end_date,".",prov,".vac_breaks.RDS"))
	
