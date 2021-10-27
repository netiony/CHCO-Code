timeseries <- Seatbelts[, "VanKilled"]
regressors <- cbind(PetrolPrice=Seatbelts[, c("PetrolPrice")],linearTrend=seq(along=timeseries)/12)
timeseries_until1981 <- window(timeseries, end=1981+11/12)
regressors_until1981 <- window(regressors, end=1981+11/12)

seatbeltsfit <- tsglm(ts=timeseries_until1981,model=list(past_obs=c(1, 12)), link="log", distr="pois",xreg=regressors_until1981)
summary(seatbeltsfit, B=100)

seatbeltsfit_alldata <- tsglm(ts=timeseries, link="log",model=list(past_obs=c(1, 12)),xreg=regressors, distr="pois")
interv_test(seatbeltsfit_alldata, tau=170, delta=1, est_interv=TRUE)
