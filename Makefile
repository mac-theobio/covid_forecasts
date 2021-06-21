## This is covid_forecasts
## https://github.com/mac-theobio/covid_forecasts

current: target
-include target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt"

######################################################################

## OMT repo

OMT:
	git clone https://github.com/wzmli/MacOMT_report.git $@

cachestuff:
	git clone https://github.com/mac-theobio/forecast_cache.git $@ 
## 

Sources += *.R *.csv
clean.Rout: clean.R
	$(pipeR)

breaks_calibrate_vac.Rout: breaks_calibrate_vac.R breaks.csv clean.rda
	$(pipeR)

betaforecast_vac.Rout: betaforecast_vac.R 
	$(pipeR)


### Makestuff

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff

## Want to chain and make makestuff if it doesn't exist
## Compress this ¶ to choose default makestuff route
Makefile: makestuff/Makefile
makestuff/Makefile:
	git clone $(msrepo)/makestuff
	ls makestuff/Makefile

-include makestuff/os.mk

## -include makestuff/pipeR.mk

-include makestuff/pipeR.mk
-include makestuff/git.mk
-include makestuff/visual.mk
