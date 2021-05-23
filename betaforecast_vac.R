library(McMasterPandemic)
library(tidyverse)
library(parallel)
library(shellpipes)
library(cowplot)
library(zoo)
library(shellpipes)

commandEnvironments()

end_date <- as.Date("2021-05-20")

flist <- list.files(path="cachestuff/",pattern=as.character(end_date))

mod <- readRDS(paste0("cachestuff/",flist))

cc <- coef(mod$fit,"fitted")
reopen_factor <- cc$rel_beta0[length(cc$rel_beta0)-2]/cc$rel_beta0[length(cc$rel_beta0)]

print(reopen_factor)

vac <- read_csv("https://data.ontario.ca/dataset/752ce2b7-c15a-4965-a3dc-397bf405e7cc/resource/8a89caa9-511c-4568-af89-7f2174b4378c/download/vaccine_doses.csv")

vaccdat <- (vac
  %>% transmute(Date = report_date+14
    , Symbol = "vacc"
    , Relative_value = (0.6*previous_day_doses_administered/14570000)/1e-10
  )
  %>% filter(complete.cases(.))
)

lift_frame <- data.frame(province = c("ON")
                         , close_date = c("2021-04-08")
                         , reopen_date = c("2021-06-02")
                         , voc_start = c("2020-12-19")
                         , pegprop = c(0.0016)
                         , scale_factor = c(15)
                         , vaccdate = c("2021-01-01")
                         , vaccrate = c(0.002)
)

lift_frame <- data.frame(province = c("ON")
                         , close_date = c("2021-04-08")
                         , reopen_date = c("2021-06-16")
                         , voc_start = c("2020-12-19")
                         , pegprop = c(0.0016)
                         , scale_factor = c(15)
                         , vaccdate = c("2021-01-01")
                         , vaccrate = c(0.002)
)


betaforecast <- function(x,voc=FALSE, close_factor=1,reopen_factor=1, Rmult=1,last_vac_factor=1,nsim=200,vacc=TRUE){
 	tempmod <- readRDS(paste0("./cachestuff/",x))
 	
 	end_date <- tempmod$fit$forecast_args$end_date
 	bd <- tempmod$fit$forecast_args$time_args$break_dates
 	last_break <- bd[length(bd)]
# 	flip_date <- as.Date((as.numeric(end_date)+as.numeric(last_break))/2)
 	flip_date <- as.Date("2021-04-23")
 	
 	## getting all the switch dates
  scale_factor <- (lift_frame %>% filter(province == tempmod$inputs$province))[,"scale_factor"]
 	voc_start <- (lift_frame %>% filter(province == tempmod$inputs$province))[,"voc_start"] %>% as.Date()
 	close_date <- (lift_frame %>% filter(province == tempmod$inputs$province))[,"close_date"] %>% as.Date()
 	reopen_date <- (lift_frame %>% filter(province == tempmod$inputs$province))[,"reopen_date"] %>% as.Date()
 	pegprop <- (lift_frame %>% filter(province == tempmod$inputs$province))[,"pegprop"] 
 	vaccdate <- (lift_frame %>% filter(province == tempmod$inputs$province))[,"vaccdate"] %>% as.Date() 
 	vaccrate <- (lift_frame %>% filter(province == tempmod$inputs$province))[,"vaccrate"] 
 	
 	
 	## flip date is when do we switch to beta_avg instead of beta_hat
 	
 	
	dateframe <- data.frame(Date = seq.Date(from = tempmod$fit$forecast_args$start_date
		, to=tempmod$fit$forecast_args$end_date + 200,by=1)
		, Symbol = "beta0"
		, Relative_value = 1
		, VoCprop = 0
	)
	
	vaccframe <- dateframe
	
  ## Create a dataframe such that the relative values are in place
	for(i in 1:length(tempmod$fit$forecast_args$time_args$break_dates)){
		for(j in 1:nrow(dateframe)){
			if(dateframe[j,"Date"] >= as.Date(tempmod$fit$forecast_args$time_args$break_dates[i])){
				dateframe[j,"Relative_value"] <- coef(tempmod$fit,"fitted")$rel_beta0[[i]]
			}
		}
	}
	
	for(i in 1:length(tempmod$fit$forecast_args$time_args$break_dates)){
	  for(j in 1:nrow(vaccframe)){
	    if(vaccframe[j,"Date"] >= as.Date(tempmod$fit$forecast_args$time_args$break_dates[i])){
	      vaccframe[j,"Relative_value"] <- coef(tempmod$fit,"fitted")$rel_vacc[[i]]
	    }
	  }
	}
	
	## If we are reopening, increase transmission by reopen_factor

	dateframe <- (dateframe
    %>% mutate(Relative_value = ifelse(Date >= as.Date(close_date),Relative_value*close_factor,Relative_value)
	    , Relative_value = ifelse(Date >= as.Date(reopen_date), Relative_value*reopen_factor,Relative_value))
	)
	
			
	## creating a dummy frame, probably don't need
	pp <- coef(tempmod$fit,"all")
	state <- make_state(params=pp)
	rr <- run_sim(pp
			, state
			, start_date=min(dateframe$Date)
			, end_date=max(dateframe$Date)
			, params_timevar = dateframe
			)
	
	cc <- coef(tempmod$fit,"fitted")
	all_cc <- coef(tempmod$fit,"all")
	logistic_fun <- function(r,t, p0=0.01){
	  p0/(p0+(1-p0)*exp(-r*t))
	}
	
	## calculate little r for each bt and bt*Rmult
	if(voc){
	  p0 = pegprop
	}
	if(!voc){
	  p0=0
	}
	
	dateframe0 <- (dateframe
	  %>% rowwise()
	  %>% mutate(ref_beta = Relative_value * as.numeric(cc$params["beta0"])
	    , VoC_beta = ref_beta*Rmult
	    , ref_r = get_r(update(all_cc,beta0=ref_beta))
	    , VoC_r = get_r(update(all_cc,beta0=VoC_beta))
	    , rdiff = VoC_r - ref_r
	    , VoCprop = ifelse(Date >= as.Date(voc_start),p0,VoCprop)
	    , bt_voc = 1
	   )
	)
	
	dateframe <- dateframe0 # %>% mutate(VoCprop = ifelse(VoCprop==p0,0.00004,VoCprop))
	
	## Stepping through and calculating bt_voc
	i <- 1
	while(i < nrow(dateframe)){
	  lf <- logistic_fun(r=dateframe$rdiff[i], t=1,p0=dateframe$VoCprop[i])
	  dateframe$bt_voc[i] <- ((1-lf)*dateframe$ref_beta[i] + lf*dateframe$VoC_beta[i])/cc$params["beta0"]
	  if(dateframe$VoCprop[i] != 0){
	  dateframe$VoCprop[i+1] <- lf
	  }
	  i <- i + 1
	}
	
	## Repegging to ~ mid point of the last bt_avg interval
	pegscale <- 1
	
	if(!is.null(flip_date)){
	  tempdf <- dateframe %>% dplyr::filter(Date == as.Date(as.character(flip_date)))
	  pegscale <- as.numeric(tempdf["bt_voc"]/tempdf["Relative_value"])
	}
	
	dateframe <- (dateframe
	  %>% mutate(tempvalue = bt_voc/pegscale
	    , Relative_value = ifelse(tempvalue>=Relative_value,tempvalue,Relative_value)
	    , final_bt = Relative_value * as.numeric(cc$params["beta0"])
	 )
	)
	dateframe = as.data.frame(dateframe)
	fa <- tempmod$fit$forecast_args
	fa$end_date <- max(dateframe$Date)
	last_vac <- mean(tail(vaccdat,7)[["Relative_value"]]) * last_vac_factor
	tempframe <- dateframe %>% select(Date,Symbol,Relative_value)
	if(vacc){
	  # vaccframe <- (tempframe
	  #               %>% select(Date)
	  #               %>% left_join(.,vaccdat)
	  #               %>% mutate(Symbol = "vacc"
	  #                          , Relative_value = ifelse(Date <= as.Date("2021-01-13"),1,Relative_value)
	  #                          , Relative_value = ifelse(Date > max(vaccdat$Date),last_vac,Relative_value)
	  #                          , Relative_value = Relative_value*0.6
	  #               )
	  # )
	  vaccframe <- (vaccframe 
	     %>% mutate(Symbol = "vacc"
	       , Relative_value = ifelse(Date > end_date+16,last_vac,Relative_value)
	       )
	     %>% select(-VoCprop)
	  )
	}
	
	tempframe <- bind_rows(tempframe,vaccframe)
	fa$params_timevar <- tempframe
	fa$end_date <- max(dateframe$Date)
	fa$sim_fun <- run_sim
	fa$time_args <- NULL
	fa$sim_args <- NULL
	fa$extra_pars <- NULL
	rr <- forecast_ensemble(tempmod$fit
									 , nsim=nsim
									 , forecast_args = fa
									, scale_Sigma = scale_factor 
									, seed = 1
	)
	
	
	rr2 <- (rr 
		%>% mutate(province = tempmod$inputs$province
					, new_strain_fraction = p0
					, voc_start = voc_start
					, reopen_date = reopen_date
					, reopen_factor = reopen_factor
					, close_date = close_date
					, close_factor = close_factor
					, vacc_factor = last_vac_factor
					, flip_date = flip_date
					, R0 = get_R0(coef(tempmod$fit,"all"))
					, pop = tempmod$inputs$population
					, vacc=vacc
		)
		%>% left_join(.,dateframe, by=c("date"="Date"))
	)
	return(rr2)
}

## No Lift
# sim0<- mclapply(flist,function(y){betaforecast(x=y,voc=FALSE, close_factor = 1,reopen_factor = 1, Rmult = 1.5,vacc=TRUE)},mc.cores=4)
sim1<- mclapply(flist,function(y){betaforecast(x=y,voc=TRUE, close_factor = 1,reopen_factor = reopen_factor, Rmult = 1.5,nsim=20,vacc=TRUE,last_vac_factor = 1)},mc.cores=4)
sim2<- mclapply(flist,function(y){betaforecast(x=y,voc=TRUE, close_factor = 1,reopen_factor = reopen_factor, Rmult = 1.5,nsim=20,vacc=TRUE,last_vac_factor = 1.5)},mc.cores=4)

sim3<- mclapply(flist,function(y){betaforecast(x=y,voc=TRUE, close_factor = 1,reopen_factor = reopen_factor, Rmult = 1.5,nsim=200,vacc=TRUE,last_vac_factor = 1)},mc.cores=4)
sim4<- mclapply(flist,function(y){betaforecast(x=y,voc=TRUE, close_factor = 1,reopen_factor = reopen_factor, Rmult = 1.5,nsim=200,vacc=TRUE,last_vac_factor = 1.5)},mc.cores=4)

sim5<- mclapply(flist,function(y){betaforecast(x=y,voc=TRUE, close_factor = 1,reopen_factor = reopen_factor, Rmult = 1.5,nsim=200,vacc=TRUE,last_vac_factor = 1)},mc.cores=4)
sim6<- mclapply(flist,function(y){betaforecast(x=y,voc=TRUE, close_factor = 1,reopen_factor = reopen_factor, Rmult = 1.5,nsim=200,vacc=TRUE,last_vac_factor = 1.5)},mc.cores=4)


betaforecast_dat <- bind_rows(sim1,sim2)
use_local_data_repo <- FALSE
source("clean.R")

betaforecast_dat2 <- (all_sub
	%>% transmute(date, var, province,obs = value)
	%>% left_join(betaforecast_dat,.)
	%>% filter(date >= as.Date("2020-09-15"))
	# %>% mutate(obstype = ifelse(date>= as.Date("2020-12-19"),"new_obs","fitted"))
	%>% mutate( #new_strain_fraction = factor(new_strain_fraction)
		Vaccination = ifelse(vacc_factor == 1, "Current","50% increase")
		)
	# %>% filter(var %in% c("report"))
)

# Ontario
Ontario_dat <- (betaforecast_dat2 %>% filter(province == "ON") %>% filter(var %in% c("report","S")))
write.csv(Ontario_dat,paste0("OMT/forecast/",end_date,"_VOC_vac.csv"))
write.csv(Ontario_dat,paste0("cachestuff/VOC_vac.csv"))


