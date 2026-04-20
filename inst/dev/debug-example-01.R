devtools::load_all(".")
library(segmented)
data(plant)
a = subset(plant,group == "RKW")
plot(y = a$y,x = a$time)
spec = dsp::dsp_spec(family = "gaussian",
                     model = "bspline",
                     times = a$time)
model = dsp::dsp_fit(y = a$y,model_spec = spec)
model$model_spec
figure_n = "Figure"

aa = model$mcmc_output
lines(y = colMeans(aa$mu), x = a$time,lwd = 4)
lines(y = apply(aa$mu,2,quantile,0.95), 
    x = a$time ,lwd = 2,lty = 2)
lines(y = apply(aa$mu,2,quantile,0.05), 
    x = a$time ,lwd = 2,lty = 2)
sum = summary(model)

model$model_spec

