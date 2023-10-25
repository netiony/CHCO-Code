#' Read IVGTT Modelling Data
#' 
#' `modelDataWideToTall()` creates data in tall format for modelling measures
#' 
#' @description This function will input the modelling data as created from a
#'              call to function `readModelData()` and return in tall format
#'              for the secretion, potentiation and dose-response variables.
#'
#' @params df   The input dataframe from a call to function `readModelData()`
#' @params x    Character string giving the name of the variable in the dataset
#'              to be converted to tall format.  Options: 'secretion', 'dose-response'
#'              or 'potentiation'   
#' @return df.tall  dataframe with the variable in tall format
#'
#' @examples
#' modeldata = readModelData()
#' secretion.long = modelDataWideToTall(df=modeldata, x='secretion')
#' potentiation.long = modelDataWideToTall(df=modeldata, x='potentiation')
#' doseresponse.long = modelDataWideToTall(df=modeldata, x='dose-response')
#' 
#' @export
modelDataWideToTall = function(df, x='secretion'){
  # Test that the parameter x is one of the valid options
  if(!(x %in% c('secretion', 'dose-response', 'potentiation'))) stop('Input parameter x is not valid')
  
  # Transform secretion data
  if(x == 'secretion'){
    index.vars = grep('^sec', names(df))
    timepoints = gsub('sec.', '', names(df[index.vars]))
    
    df.tall = reshape(df,
                      direction = 'long',
                      varying = names(df[index.vars]), 
                      v.names = 'secretion',
                      timevar = 'timepoint',
                      times   = as.integer(timepoints))
    
    df.tall = df.tall[,!(colnames(df.tall) %in% c(grep(("^pot."), names(modeldata), value = TRUE), grep("^dos.", names(modeldata), value=TRUE))) ]
    
    attr(df.tall$secretion, 'label') = "Insulin Secretion (pmol/min/m^2)"
    attr(df.tall$secretion, 'units') = "pmol/min/m^2"
    attr(df.tall$timepoint, 'label') = "Timepoint of IVGTT (min)"
    attr(df.tall$timepoint, 'units') = "min"
    
  }
  
  # Transform potentiation data
  if(x == 'potentiation'){
    index.vars = grep('^pot', names(df))
    timepoints = gsub('pot.', '', names(df[index.vars]))
    
    df.tall = reshape(df,
                      direction = 'long',
                      varying = names(df[index.vars]), 
                      v.names = 'potentiation',
                      timevar = 'timepoint',
                      times   = as.integer(timepoints))
    
    df.tall = df.tall[,!(colnames(df.tall) %in% c(grep(("^sec."), names(modeldata), value = TRUE), grep("^dos.", names(modeldata), value=TRUE)))]
  
    
    attr(df.tall$potentiation, 'label') = "Potentiation"
    attr(df.tall$potentiation, 'units') = ""
    attr(df.tall$timepoint, 'label') = "Timepoint of IVGTT (min)"
    attr(df.tall$timepoint, 'units') = "min"
    
  }
  
  
  # Transform dose-response data
  if(x == 'dose-response'){
    index.vars = grep('^dos', names(df))
    timepoints = gsub('dos.', '', names(df[index.vars]))
    
    df.tall = reshape(df,
                      direction = 'long',
                      varying = names(df[index.vars]), 
                      v.names = 'dose.response',
                      timevar = 'gluc.concentration',
                      times   = as.numeric(timepoints))
    
    df.tall = df.tall[,!(colnames(df.tall) %in% c(grep(("^pot."), names(modeldata), value = TRUE), grep("^sec.", names(modeldata), value=TRUE)))]
    
    attr(df.tall$dose.response, 'label') = "Insulin Secretion (pmol min-1m-2)"
    attr(df.tall$dose.response, 'units') = "pmol min-1m-2"
    attr(df.tall$gluc.concentration, 'label') = "Glucose loading concentration (mmol)"
    attr(df.tall$gluc.concentration, 'units') = "mmol"
  }
  
  return(df.tall)
  
}