###Examples#################


#First example:

I1<-int_slow_periodic(function(t) 1/t,c(1),c(0,1,0),100,2)

integrate(function(t){1/t*(cos(100*t))},1,2)

Int.int_slow_periodic(I1,1,2)

#Second example

I2<- int_slow_periodic(function(t) t^2,c(1,2),c(1,0.5,0,0,1),100,2)

integrate(function(t)t^2*(1+0.5*cos(100*t)+sin(2*100*t)),1,2)

#As expected the error is around 1%
(Int.int_slow_periodic(I2,1,2)-integrate(function(t)t^2*(1+0.5*cos(100*t)+sin(2*100*t)),1,2)$value)/integrate(function(t)t^2*(1+0.5*cos(100*t)+sin(2*100*t)),1,2)$value
