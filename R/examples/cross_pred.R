library(splines)
library(dlnm)
library(data.table)
devtools::load_all()

### example of application in time series analysis - see vignette("dlnmTS")

df = data.table(chicagoNMMAPS)
# create the crossbasis objects and summarize their contents
cb1.pm <- crossbasis(chicagoNMMAPS$pm10,
  lag = 15, argvar = list(fun = "lin"),
  arglag = list(fun = "poly", degree = 4)
)
cb1.temp <- crossbasis(chicagoNMMAPS$temp,
  lag = 3, argvar = list(df = 5),
  arglag = list(fun = "strata", breaks = 1)
)
summary(cb1.pm)
summary(cb1.temp)

# run the model and get the predictions for pm10
model1 <- glm(death ~ cb1.pm + cb1.temp + ns(time, 7 * 14) + dow,
  family = quasipoisson(), chicagoNMMAPS
)
pred1.pm <- crosspred(cb1.pm, model1, at = 0:20, bylag = 0.2, cumul = TRUE)
str(pred1.pm)

# plot the lag-response curves for specific and incremental cumulative effects
plot(pred1.pm, "slices",
  var = 10, col = 3, ylab = "RR", ci.arg = list(density = 15, lwd = 2),
  main = "Lag-response curve for a 10-unit increase in PM10"
)

plot(pred1.pm, "slices",
  var = 10, col = 2, cumul = TRUE, ylab = "Cumulative RR",
  main = "Lag-response curve of incremental cumulative effects"
)
