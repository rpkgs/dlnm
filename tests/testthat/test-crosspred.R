library(splines)
library(dlnm)
library(data.table)

test_that("crosspred works", {
  expect_equal(2 * 2, 4)
})

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
colnames(pred1.pm$matfit)
