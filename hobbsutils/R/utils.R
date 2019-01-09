.onAttach <- function(libname, pkgname) {
  packageStartupMessage("><(((O> ><(((O> ><(((O> Welcome to the Hobbslab utility package <O)))>< <O)))>< <O)))>< ")

}

ec_2_sal = function(temp, cond){
  if (any(temp > 35)) {
    warning('Temperature is high, ensure that units are in degrees C', call. = F, immediate. = T)
  }
  ref_cond = 42914
  cond_rat = cond/ref_cond
  rt = 0.6766097 + (0.0200564*temp) + (0.0001104259*(temp^2)) + ((-6.9698*10^-7)*(temp^3)) + ((1.0031*10^-9)*temp^4)
  Rt = cond_rat/rt
  dS = ((temp-15)/(1+0.0162*(temp-15)))*(0.0005+(-0.0056)*(Rt^0.5)+(-0.0066)*Rt+(-0.0375)*(Rt^1.5)+(0.0636)*(Rt^2)+(-0.0144)*(Rt^2.5))
  sal = ifelse(is.na(cond), NA,
         ifelse (cond > 3000, 0.008+ (-0.1692)*(Rt^0.5)+25.3851*Rt+14.0941*(Rt^1.5)+(-7.0261)*(Rt^2)+2.7081*(Rt^2.5)+dS,
                 (0.008+ (-0.1692)*(Rt^0.5)+25.3851*Rt+14.0941*(Rt^1.5)+(-7.0261)*(Rt^2)+2.7081*(Rt^2.5)+dS)-(0.008/(1+(1.5*(400*Rt))+((400*Rt)^2))-(0.0005*(temp-15)/(1+0.0162*(temp-15)))/(1+((100*Rt)^0.5)+((100*Rt)^1.5)))))
  return(sal)
}


sr_2_sal = function(sr, srfw = 0.705264, srmar = 0.70918,confw = 74.6, conmar = 6819,salfw = 0.1,salmar = 31.8){
  if (any(sr < min(srfw, srmar)| sr > max(srfw, srmar))) {
    warning('Some of your measured strontium ratio values are outside the bounds of your two endmembers, make sure that srfw and srmar are set correctly',
            call. = F, immediate. = T)
  }
    sal = (((salfw*srmar*conmar) - (salfw*sr*conmar) - (salmar*srmar*conmar) + (salmar*sr*conmar))/
             ((sr*confw) - (sr*conmar) - (srfw*confw) + (srmar*conmar))) + salmar
    return(sal)

}

sal_2_sr = function(sal, srfw = 0.705264, srmar = 0.70918, confw = 74.6, conmar = 6819, salfw = 0.1, salmar = 31.8) {
  if(any(sal < min(salfw, salmar) | sal > max(salfw, salmar))) {
    warning('Some of your measured salinity values are outside the bounds of your two endmembers, make sure that salfw and salmar are set correctly',
            call. = F, immediate. = T)
  }
    sr = ((((srfw*confw*sal)-(srfw*confw*salmar))/(salfw - salmar))+(srmar*conmar)-(((srmar*conmar*sal)-(srmar*conmar*salmar))/(salfw - salmar)))/
      ((((confw*sal)-(confw*salmar))/(salfw - salmar))+(conmar)-(((conmar*sal)-(conmar*salmar))/(salfw-salmar)))
    return(sr)

}

bim = function(fl, hl, gt) {
  fl + ((gt - max(gt))*(fl - hl))/(max(gt) - min(gt))
}
