.onAttach <- function(libname, pkgname) {
  packageStartupMessage("><(((O> ><(((O> ><(((O> Welcome to the Hobbslab utility package <O)))>< <O)))>< <O)))>< ")

}


sr_2_sal = function(sr, srfw = 0.705264, srmar = 0.70918,confw = 74.6, conmar = 6819,salfw = 0.1,salmar = 31.8){
  # if(sr < min(srfw, srmar)| sr > max(srfw, srmar) {
  #   warning('Your measured strontium ratio is outside the bounds of your two endmembers, make sure that srfw and srmar are set correctly',
  #           call. = F, immediate. = T)
  #   return(NULL)
  # } else{
  sal = (((salfw*srmar*conmar) - (salfw*sr*conmar) - (salmar*srmar*conmar) + (salmar*sr*conmar))/
           ((sr*confw) - (sr*conmar) - (srfw*confw) + (srmar*conmar))) + salmar
  return(sal)
}

sal_2_sr = function(sal, srfw = 0.705264, srmar = 0.70918, confw = 74.6, conmar = 6819, salfw = 0.1, salmar = 31.8) {
  # if(sal < min(salfw, salmar) | sal > max(salfw, salmar)) {
  #   warning('Your measured salinity is outside the bounds of your two endmembers, make sure that salfw and salmar are set correctly',
  #           call. = F, immediate. = T)
  #   return(NULL)
  # } else {
  sr = ((((srfw*confw*sal)-(srfw*confw*salmar))/(salfw - salmar))+(srmar*conmar)-(((srmar*conmar*sal)-(srmar*conmar*salmar))/(salfw - salmar)))/
    ((((confw*sal)-(confw*salmar))/(salfw - salmar))+(conmar)-(((conmar*sal)-(conmar*salmar))/(salfw-salmar)))
  return(sr)
  # }
}

bim = function(fl, hl, gt) {
  fl + ((gt - max(gt))*(fl - hl))/(max(gt) - min(gt))
}
