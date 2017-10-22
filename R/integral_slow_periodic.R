#Please read the README file in order to understand the code and the inputs required.



vect_sinu<- function(modes,modes_coef, omega, n){

  #Creates a vector with the sines/cosines of the Fourier modes involved.

  n_modes <- length(modes)

  vect_c <- sapply(1:n_modes, function(x) paste0(modes_coef[2*x],'*','cos(',modes[x],'*',omega,'*t)',
                                                 '-', modes_coef[2*x+1],'*','cos(',modes[x],'*',omega,'*t+pi/2)'))
  vect_s <- sapply(1:n_modes, function(x) paste0(modes_coef[2*x],'*','sin(',modes[x],'*',omega,'*t)',
                                                 '-', modes_coef[2*x+1],'*','sin(',modes[x],'*',omega,'*t+pi/2)'))

  return(list(parse(text = vect_c),parse(text = vect_s)))
}


matrix_sinu <- function(modes, omega, n){

  #Creates a matrix with the coefficients that multiply the sines/cosines
  #that corresponds to the integration coefficients of the periodic part when integrating by parts.

  n_modes <- length(modes)

  matrix_cos <- matrix(1:n*n_modes, ncol = n, nrow = n_modes)
  matrix_sin <- matrix(1:n*n_modes, ncol = n, nrow = n_modes)

  for(k in 1:(n)){
    for(j in 1:n_modes){
      matrix_cos[j,k] <- (-1)^(k+1)*Re((1i)^(k))/(modes[j]*omega)^(k)
      matrix_sin[j,k] <- -(-1)^(k+1)*Re((1i)^(k+1))/(modes[j]*omega)^(k)
    }
  }
  return(list(matrix_cos,matrix_sin))
}

I_periodic_part <- function(modes,modes_coef, omega, n){
  #Builds the integrals of the periodic part.
  ms <- matrix_sinu(modes,omega, n)
  vs <- vect_sinu(modes,modes_coef, omega, n)

  return(function(t){t(ms[[1]]*eval(vs[[1]]) +ms[[2]]*eval(vs[[2]]))})
}

I_slow_part<- function(f, n){
  #Creates a vector with the form of the slow-part for each iteration of the integration by parts.

  slow_f <-body(f)

  vec_slow_f <- c(slow_f)

  if(n!=1){
    for(k in 1:(n-1)){
      print(k)
      slow_f <- D(slow_f,'t')
      vec_slow_f <- c(vec_slow_f,slow_f)
    }
  }
  return(vec_slow_f)
}




int_slow_periodic <- function(f,modes,modes_coef, omega, n){
  #Creates a class with the relevant information to solve the integrals
  #with the form F(t)P(t) as it is explained in the documentation

  periodic<- I_periodic_part(modes,modes_coef, omega, n)
  slowly <- I_slow_part(f,n)

  values<- list(f=f,modes=modes,modes_coef=modes_coef,
                omega=omega, n=n,periodic=periodic, slowly=slowly)

  attr(values,'class')<- 'int_slow_periodic'
  values
}

Int.int_slow_periodic <- function(obj,ti,tf){
  #Method to finally solve the integral between ti and tf.
  I0 <- obj$modes_coef[1]*integrate(obj$f,ti,tf)$value
  IP <- function(obj,t){
    m<- length(obj$slowly)
    v_slowly <- 1:m
    for(i in 1:m){
      v_slowly[i]<- eval(obj$slowly[[i]])
    }
    return(sum(sapply(t,obj$periodic)*v_slowly))}


  return(I0+IP(obj,tf)-IP(obj,ti))
}



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


