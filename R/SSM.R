

#########################################################
## r
## based on Schoolfield, R. M., Sharpe, P. J. & Magnuson, C. E. 
## 1981. Non-linear regression of biological temperature-dependent 
## rate models based on absolute reaction-rate theory. Journal of 
## Theoretical Biology, 88, 719-731.
##########################################################

.SSM<-function (T, parms) {

	if (all(names(parms)!="Rho25")) {
	  newx <- abs(parms[(names(parms)!="rK") & (names(parms)!="K") & (names(parms)!="Scale")])
    scale <- ifelse(is.na(parms["Scale"]), 1, parms["Scale"])
	  tableT <- data.frame(Temperature=as.numeric(names(newx)), R=newx)
	  ml <- loess(R ~ Temperature, data=tableT)
	  newtable<- data.frame(Temperature=T)
	  rT <- abs(predict(ml, newtable)*1E-5)*scale
	  rT_L <- rT
	
	} else {

	R <- 8.314472
	rho25 <- parms["Rho25"]/1E7
	dha <- parms["DHA"]*1E3
	dhh <- parms["DHH"]*1E3
  unsurT <- 1/T
  unsur298 <- 1/298

	if (is.na(parms["DHL"])) {
# Je suis en 4 paramètres
#	"T12H", "DHA",  "DHH", "Rho25"
		t12H=abs(parms["T12H"])
		rT<-(rho25*(T/298)*exp((dha/R)*(unsur298-unsurT)))/(1+exp((dhh/R)*((1/t12H)-unsurT)))
	} else {
# 28/7/2012 - T12H changé en DT	
#	"T12L", "DT", "DHA",  "DHH", "DHL", "Rho25"
		dhl <- parms["DHL"]*1E3
		t12L <- parms["T12L"]
		t12H <- t12L+abs(parms["DT"])
		rT<-(rho25*(T/298)*exp((dha/R)*(unsur298-unsurT)))/(1+exp((dhl/R)*((1/t12L)-unsurT))+exp((dhh/R)*((1/t12H)-unsurT)))
	}

	rT_L <- rT

	if (!is.na(parms["Rho25_L"])) {
	
	rho25_L <- parms["Rho25_L"]/1E7
	dha_L <- parms["DHA_L"]*1E3
	dhh_L <- parms["DHH_L"]*1E3

	if (is.na(parms["DHL_L"])) {
#	"T12H", "DHA",  "DHH", "Rho25"
		t12H_L <- abs(parms["T12H_L"])
		rT_L <- (rho25_L*(T/298)*exp((dha_L/R)*(unsur298-unsurT)))/(1+exp((dhh_L/R)*((1/t12H_L)-unsurT)))
	} else {
# 28/7/2012 - T12H changé en DT	
#	"T12L", "T12H", "DHA",  "DHH", "DHL", "Rho25"
		dhl_L <- parms["DHL_L"]*1E3
		t12L_L <- parms["T12L_L"]
		t12H_L <- t12L+abs(parms["DT_L"])
		rT_L <- (rho25_L*(T/298)*exp((dha_L/R)*(unsur298-unsurT)))/(1+exp((dhl_L/R)*((1/t12L_L)-unsurT))+exp((dhh_L/R)*((1/t12H_L)-unsurT)))
	}
	
	
	
	
	}
	}

return(list(as.numeric(rT), as.numeric(rT_L)))

}
