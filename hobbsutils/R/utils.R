.onAttach <- function(libname, pkgname) {
  packageStartupMessage("><(((O> ><(((O> ><(((O> Welcome to the Hobbslab utility package <O)))>< <O)))>< <O)))>< ")

}

ec_2_sal = function(temp, cond){
  if (any(temp > 35, na.rm = T)) {
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

sal_2_spc = function(sal){
  if (any(sal > 40, na.rm = T)){
    warning("Salinity is high, ensure that data is correct")
  }
  J1 = -16.072
  J2 = 4.1495
  J3 = -0.5345
  J4 = 0.0261
  spc = (sal/35)*(53087) + sal*(sal-35) *
    (J1 + (J2 * (sal^(0.5))) + (J3*sal) + (J4*(sal^(1.5))))
  return(spc)
}


sr_2_sal = function(sr, srfw = 0.705264, srmar = 0.70918,confw = 74.6, conmar = 6819,salfw = 0.1,salmar = 31.8, sallim, fill = "NA"){
  if (any(sr < min(srfw, srmar)| sr > max(srfw, srmar), na.rm = T)) {
    warning('Some of your measured strontium ratio values are outside the bounds of your two endmembers, make sure that srfw and srmar are set correctly',
            call. = F, immediate. = T)
  }
  if (is.na(sallim)) stop('You have not set a high salinity limit (sallim argument)')
    sal = (((salfw*srmar*conmar) - (salfw*sr*conmar) - (salmar*srmar*conmar) + (salmar*sr*conmar))/
             ((sr*confw) - (sr*conmar) - (srfw*confw) + (srmar*conmar))) + salmar
    sal = ifelse(sal > sallim, sallim, sal)
    sal = ifelse(sr > srmar & fill == 'NA', NA, ifelse(sr>srmar & fill == 'sallim',sallim,sal))
    return(sal)

}

sal_2_sr = function(sal, srfw = 0.705264, srmar = 0.70918, confw = 74.6, conmar = 6819, salfw = 0.1, salmar = 31.8) {
  if(any(sal < min(salfw, salmar) | sal > max(salfw, salmar), na.rm = T)) {
    warning('Some of your measured salinity values are outside the bounds of your two endmembers, make sure that salfw and salmar are set correctly',
            call. = F, immediate. = T)
  }
    sr = ((((srfw*confw*sal)-(srfw*confw*salmar))/(salfw - salmar))+(srmar*conmar)-(((srmar*conmar*sal)-(srmar*conmar*salmar))/(salfw - salmar)))/
      ((((confw*sal)-(confw*salmar))/(salfw - salmar))+(conmar)-(((conmar*sal)-(conmar*salmar))/(salfw-salmar)))
    return(sr)

}

o2_2_sal = function(oxy_rat, source = 'ingram') {
  sal = ifelse(source == 'ingram',((oxy_rat+10.9)/0.32),
               ifelse(source == 'mclg', ((oxy_rat+10.17)/0.29),
                      'improper source selected'))
  return(sal)
}

vsmow_2_vpdb = function(x) {
  return((0.97001*x)-29.99)
}

vpdb_2_vsmow = function(x) {
  return((1.03091*x)+30.91)
}

bim = function(fl, hl, gt) {
  fl + ((gt - max(gt))*(fl - hl))/(max(gt) - min(gt))
}

membermix = function(sr, conc, sal, mix) {
  if(sum(mix) != 1) {
    warning('Your mixture does not sum to 100%',
            call. = F, immediate. = T)
  }
  srmix = sum(sr*conc*mix)/sum(conc*mix)
  srconc = sum(conc*mix)
  salmix = sum(sal*mix)
  return(list(sr = srmix, conc = srconc, sal = salmix))
}

l2l = function(from,to,measurement,lengths,species){
  if(species %in% c('LONSME', 'DELSME')) {
    if(measurement == 'SL') {
      if(from == 'ETOH' & to == 'FIELD') {
        calclength = (lengths* 1.05793) + 0.97646
      } else if(from == 'FIELD' & to == 'ETOH') {
        calclength = (lengths * 0.91736) + 97900
      }
    } else if(measurement== 'FL') {
      if(from == 'ETOH' & to == 'FIELD') {
        calclength = (lengths* 1.0225) + 0.336769
      } else if(from == 'FIELD' & to == 'ETOH') {
        calclength = (lengths * 0.972129) + 0.227886
      }
    }
  }
  return(calclength)
}
