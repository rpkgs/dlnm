x = chicagoNMMAPS$temp
r <- crossbasis(x, lag = 3,
  argvar = list(df = 5), # ns
  arglag = list(fun = "strata", breaks = c(0, 10))
)

# summary(cb1.pm)
str(r)

onebasis(x, "ns", df = 5) %>% str()
onebasis(x, "strata", breaks = c(0, 10)) %>% str()
