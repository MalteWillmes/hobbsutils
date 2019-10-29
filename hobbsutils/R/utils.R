.onAttach <- function(libname, pkgname) {
  packageStartupMessage("><(((O> ><(((O> ><(((O> Welcome to the Hobbslab utility package <O)))>< <O)))>< <O)))>< ")

}

ec_2_sal = function(temp, cond){
  if (any(temp > 35, na.rm = T)) {
    warning('Temperature is high, ensure that units are in degrees C', call. = F, immediate. = T)
  }
  if(length(temp == 1)) {
    temp = rep(temp, length.out = length(cond))
  } else if ( length(temp == length(cond))) {
    temp = temp
  } else {
    warning('Temperature must be of length 1 or match number of conductivity readings')
  }
  ref_cond = 42914
  cond_rat = cond/ref_cond
  rt = 0.6766097 + (0.0200564*temp) + (0.0001104259*(temp^2)) + ((-6.9698*10^-7)*(temp^3)) + ((1.0031*10^-9)*temp^4)
  Rt = cond_rat/rt
  dS = ((temp-15)/(1+0.0162*(temp-15)))*(0.0005+(-0.0056)*(Rt^0.5)+(-0.0066)*Rt+(-0.0375)*(Rt^1.5)+(0.0636)*(Rt^2)+(-0.0144)*(Rt^2.5))
  sal = c()
  sal[is.na(cond)] = NA
  sal[cond >3000 & !is.na(cond)] = 0.008 + (-0.1692)*(Rt[cond>3000 & !is.na(cond)]^0.5) + 25.3851*Rt[cond > 3000 & !is.na(cond)]+14.0941*(Rt[cond > 3000 & !is.na(cond)]^1.5)+(-7.0261)*(Rt[cond > 3000 & !is.na(cond)]^2)+2.7081*(Rt[cond > 3000 & !is.na(cond)]^2.5)+dS[cond > 3000 & !is.na(cond)]
  sal[cond <=3000 & !is.na(cond)] = (0.008+ (-0.1692)*(Rt[cond <= 3000 & !is.na(cond)]^0.5)+25.3851*Rt[cond <= 3000 & !is.na(cond)]+14.0941*(Rt[cond <= 3000 & !is.na(cond)]^1.5)+(-7.0261)*(Rt[cond <= 3000 & !is.na(cond)]^2)+2.7081*(Rt[cond <= 3000 & !is.na(cond)]^2.5)+dS[cond <= 3000 & !is.na(cond)]) -
    (0.008/(1+(1.5*(400*Rt[cond <= 3000 & !is.na(cond)]))+((400*Rt[cond <= 3000 & !is.na(cond)])^2))-(0.0005*(temp[cond <= 3000 & !is.na(cond)]-15)/(1+0.0162*(temp[cond <= 3000 & !is.na(cond)]-15)))/(1+((100*Rt[cond <= 3000 & !is.na(cond)])^0.5)+((100*Rt[cond <= 3000 & !is.na(cond)])^1.5)))
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
    if(fill == 'NA') {
      sal[sal > sallim] = NA
    } else if (fill == 'sallim') {
      sal [sal > sallim] = sallim
    }
    sal[sr < (srfw - 0.00005) | sr > srmar + 0.00005] = NA
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
  if(source == 'ingram') {
    sal = (oxy_rat + 10.9)/0.32
  } else if(source == 'mclg') {
    sal = (oxy_rat + 10.17)/0.29
  } else{
    warning('Non-supported source')
  }
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


l2l = function(from, to = 'flf',length) {
  if(from == 'slf' & to == 'flf') {
    calclength = 1.09*(length)-0.97
  } else if (from == 'slf' & to == 'tlf') {
    calclength = 1.21*length-2.17
  } else if (from == 'tlf' & to == 'flf') {
    calclength = (length-0.36)/1.08
  } else if (from == 'tlf' & to == 'slf') {
    calclength = (length+2.17)/1.21
  } else if (from == 'flf' & to == 'slf') {
    calclength = (length+0.97)/1.09
  } else if (from == 'flf' & to == 'tlf') {
    calclength = 1.08*length+0.36
  } else if (from == 'sle' & to == 'flf') {
    intercalc = (length-0.64)/0.94
    calclength = 1.09*(intercalc)-0.97
  } else if (from == 'sle' & to == 'tlf') {
    intercalc = (length-0.64)/0.94
    calclength = 1.21*intercalc-2.17
  } else if (from == 'sle' & to == 'slf') {
    calclength = (length-0.64)/0.94
  } else if (from == 'fle' & to == 'flf') {
    calclength = (length-3.09)/0.94
  } else if (from == 'fle' & to == 'slf') {
    intercalc = (length-3.09)/0.94
    calclength = (intercalc+0.97)/1.09
  } else if (from == 'fle' & to == 'tlf') {
    intercalc = (length-3.09)/0.94
    calclength = 1.08*intercalc+0.36
  } else {
    stop('Unsupported conversion', call. = F)
  }
  if(any(length < 40, na.rm = T)) {
    warning("Some lengths outside of established conversion bounds", call. = F)
  }
  return(round(calclength,1))
}


colmatch <- function(x,y, match = F, join = F){
  if (match == T & join == F) {
    if (all(colnames(x) %in% colnames(y)) & all(colnames(y) %in% colnames(x))) {
      return('All column names match between both objects')
    } else if (any(colnames(x) %in% colnames(y))) {
      cat("\n\n",'columns ')
      cat(colnames(x)[colnames(x) %in% colnames(y)], sep = " and ")
      cat( " are found in both objects\n")
    } else if(!any(colnames(x) %in% colnames(y))) {
      cat("\n\n",'No columns ')
      cat(colnames(x)[!colnames(x) %in% colnames(y)], sep = " and ")
      cat( " are matching between both objects\n")
    }
  } else if(match == F & join == F) {
    if (all(colnames(x) %in% colnames(y)) & all(colnames(y) %in% colnames(x))) {
      return('All column names match between both objects')
    } else if (!all(colnames(x) %in% colnames(y))
               & all(colnames(y) %in% colnames(x))) {
      if(length(colnames(x)[!colnames(x) %in% colnames(y)]) > 1) {
        cat("\n\n",'columns ')
        cat(colnames(x)[!colnames(x) %in% colnames(y)], sep = " and ")
        cat( " are missing from second object\n")
      } else {
        cat("\n\n",'column ')
        cat(colnames(x)[!colnames(x) %in% colnames(y)], sep = " and ")
        cat( " is missing from second object\n")
      }
    } else if (all(colnames(x) %in% colnames(y)) & !all(colnames(y) %in% colnames(x))){
      if(length(colnames(y)[!colnames(y) %in% colnames(x)]) > 1) {
        cat("\n\n",'columns ')
        cat(colnames(y)[!colnames(y) %in% colnames(x)], sep = " and ")
        cat( " are missing from first object\n")
      } else {
        cat("\n\n",'column ')
        cat(colnames(y)[!colnames(y) %in% colnames(x)], sep = " and ")
        cat( " is missing from first object\n")
      }
    } else if (!all(colnames(x) %in% colnames(y)) & !all(colnames(y) %in% colnames(x))){
      if(length(colnames(y)[!colnames(y) %in% colnames(x)]) > 1 & length(colnames(x)[!colnames(x) %in% colnames(y)]) > 1) {
        cat("\n\n",'columns ')
        cat(colnames(y)[!colnames(y) %in% colnames(x)], sep = " and ")
        cat(' are missing from first object \n and')
        cat('\n','columns ')
        cat(colnames(x)[!colnames(x) %in% colnames(y)], sep = " and ")
        cat( " are missing from second object\n")
      } else if(length(colnames(y)[!colnames(y) %in% colnames(x)]) > 1 & length(colnames(x)[!colnames(x) %in% colnames(y)]) == 1){
        cat("\n\n",'columns ')
        cat(colnames(y)[!colnames(y) %in% colnames(x)], sep = " and ")
        cat(' are missing from first object \n and')
        cat('\n','column ')
        cat(colnames(x)[!colnames(x) %in% colnames(y)], sep = " and ")
        cat( " is missing from second object\n")
      } else if(length(colnames(y)[!colnames(y) %in% colnames(x)]) == 1 & length(colnames(x)[!colnames(x) %in% colnames(y)]) > 1){
        cat("\n\n",'column ')
        cat(colnames(y)[!colnames(y) %in% colnames(x)], sep = " and ")
        cat(' is missing from first object \n and')
        cat('\n','columns ')
        cat(colnames(x)[!colnames(x) %in% colnames(y)], sep = " and ")
        cat( " are missing from second object\n")
      } else {
        cat("\n\n",'column ')
        cat(colnames(y)[!colnames(y) %in% colnames(x)], sep = " and ")
        cat(' is missing from first object \n and')
        cat('\n','column ')
        cat(colnames(x)[!colnames(x) %in% colnames(y)], sep = " and ")
        cat( " is missing from second object\n")
      }
    }
  } else if (match == T & join == T) {
    colnames(x[colnames(x) %in% colnames(y)])
    } else {
    cat('Imporper selection of match argument, select either "TRUE/T" or "FALSE/F. If join is set to true, match must also be T)')
  }

}
