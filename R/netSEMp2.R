##' This function builds an netSEM model using principle 2. 
##' 
##' Principle 2 resembles the multiple regression principle in the way multiple predictors are considered simultaneously. Specifically, the first-level predictors to the system level variable, such as, Time and unit level variables, acted on the system level variable collectively by an additive model. This collective additive model can be found with a generalized stepwise variable selection (using the stepAIC() function in R, which performs variable selection on the basis of AIC or BIC) and this proceeds iteratively.
##' 
##' Data is analysed first using Principle 1 to find the best models. If needed, transformations based on the best models are applied to the predictors. Starting from the system response variable, each variable is regressed on all other variables except for the system response in an additive multiple regression model, which is reduced by a stepwise selection using stepAIC(). Then, for each selected variable, fitted regression for 6 selected functional forms (8 if strength type problem) and pick the best.
##'
##' @title network Structural Equation Modelling (netSEM) - Principle 2


##' @param x A dataframe. By default it considers all columns as exogenous variables, except the first column which stores the system endogenous variable.
##' @param exogenous, by default it consideres all columns as exogenous variables except column number 1, which is the main endogenous response.
##' @param endogenous A character string of the column name of the main endogenous OR a numeric number indexing the column of the main endogenous.
##' @param str A boolean, whether or not this is a 'strength' type problem.
##' @param criterion By default uses AIC to identify best model. Alternative choice is BIC.
##' @return A list of the following items:
##'
##' \itemize{
##' \item "data":
##' \item "res.print": A dataframe.
##' \item "res.model": A model object.
##'}

##' @seealso \link[netSEM]{netSEMp1} 
##'
##' @export
##' 
##' 
##' @examples
##' \dontrun{
##' ## Load the sample acrylic data set
##' data(acrylic)
##' 
##' ## Using AIC criterion
##' ans <- netSEMp2(acrylic, "IrradTot", "YI")
##' 
##' ## Using BIC criterion
##' ans <- netSEMp2(acrylic, "IrradTot", "YI", criterion = "BIC")
##' }

####### v0.7.0 ##################
netSEMp2 <-
  function(x,
           exogenous,
           endogenous,
           str = FALSE,
           criterion = "AIC") {
    # Default is when str = FALSE (not considering strength type problem)
    
    # Check for exogenous & endogenous input. Error if missing.
    # Rearrange dataframe based on exogenous & endogenous specification.
    # Transformed dataframe structure is: Stressor, Response,...
    if (missing(exogenous) | missing(endogenous)) {
      stop("exogenous and endogenous variables need to be specified")
    } else {
      x <- x %>%
        dplyr::relocate(all_of(endogenous), all_of(exogenous))
    }
    x <<- x
    # Run principle 1
    p1.result <- netSEMp1(x, exogenous, endogenous)
    p1.bestModels <- p1.result[["bestModels"]]
    p1.allModels <- p1.result[["allModels"]]
    
    Resp <- p1.bestModels$Resp
    unique_Resp <- unique(Resp)
    
    lm.best.model <- list()
    lm.best.text.final <- list()
    
    Res.print.local <- matrix(NA, nrow = 0, ncol = 9)
    colnames(Res.print.local) <- c("endogenous",  "Variable",  "GModel", 
                                   "GAdj-R2",       "IModel",  "IAdj-R2",
                                   "Gcoeff",       "GModelB",  "Model" )
    
    Res.print <- matrix(NA, nrow = 0, ncol = 9)
    colnames(Res.print) <- c("endogenous",  "Variable",  "GModel", 
                             "GAdj-R2",       "IModel",  "IAdj-R2",
                             "Gcoeff",       "GModelB",  "Model" )
    
    ################## START: APPLY TRANSFORMATION TO PREDICTORS #################
    
    
    for ( numb_Resp in 1: length(unique_Resp) ) {
      x <- p1.result$data
      # New dataframe to save transformed columns
      x1 <- x
      # Filter principle 1 best model dataframe when Response is unique_Resp[numb_Resp]
      p1.bestModels.new <- p1.bestModels %>% 
        dplyr::filter(Resp == unique_Resp[numb_Resp])
      unique_Var <- unique(p1.bestModels.new$Var)
      trans <- rep(NA, length = length(unique_Var))
      trans1 <- rep(NA, length = length(unique_Var))
      names(trans) <- unique_Var
      
      for ( i in 1: sum(Resp == unique_Resp[numb_Resp]) ) {
        # If "SL", x ---> x
        if ( p1.bestModels.new$best_model[i] == "Linear" ) {
          trans1[i] <- NA 
          names(trans1)[i] <- unique_Var[i]
          # Find coefficients
          for (ii in 1:nrow(p1.bestModels)) {
            if ((p1.allModels[[ii]]$exogenous) == unique_Var[i] & (p1.allModels[[ii]]$endogenous) == unique_Resp[numb_Resp]) {
              Quad_value <- coef(p1.allModels[[ii]]$Linear)
            }
          }
          if ( Quad_value[2]>0 ) {
            trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),"+",
                              format(Quad_value[2],scientific=FALSE,digits=2),"*",unique_Var[i],sep="")
          } else {
            trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),
                              format(Quad_value[2],scientific=FALSE,digits=2),"*",unique_Var[i],sep="")
          }
        }
        
        # If "SQuad", x ---> x^2
        if ( p1.bestModels.new$best_model[i] == "SQuad" ) {
          trans1[i] <- paste(unique_Var[i],"_t=", unique_Var[i],"^2",sep="")
          names(trans1)[i] <- unique_Var[i]
          x1[,ncol(x1)+1] <- (x1[,unique_Var[i]])^2
          names(x1)[ncol(x1)] <- paste(unique_Var[i],"_t",sep="")
          # Find coefficients
          for (ii in 1:nrow(p1.bestModels)) {
            if ((p1.allModels[[ii]]$exogenous) == unique_Var[i] & (p1.allModels[[ii]]$endogenous) == unique_Resp[numb_Resp]) {
              Quad_value <- coef(p1.allModels[[ii]]$SQuad)
            }
          }
          if ( Quad_value[2]>0 ) {
            trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),"+",
                              format(Quad_value[2],scientific=FALSE,digits=2),"*",unique_Var[i],"^2",sep="")
          } else {
            trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),
                              format(Quad_value[2],scientific=FALSE,digits=2),"*",unique_Var[i],"^2",sep="")
          }
        }
        
        
        # If "Exp", x ---> exp(x)
        if ( p1.bestModels.new$best_model[i] == "Exp" ) {
          trans1[i] <- paste(unique_Var[i],"_t=exp(", unique_Var[i],")",sep="")
          names(trans1)[i] <- unique_Var[i]
          x1[,ncol(x1)+1] <- exp(x1[,unique_Var[i]])
          names(x1)[ncol(x1)] <- paste(unique_Var[i],"_t",sep="")
          Quad_value <- coef(p1.allModels[[i]]$Exponential)
          if ( Quad_value[2]>0 ) {
            trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),"+",
                              format(Quad_value[2],scientific=FALSE,digits=2),"*e^",unique_Var[i],sep="")
          } else {
            trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),
                              format(Quad_value[2],scientific=FALSE,digits=2),"*e^",unique_Var[i],sep="")
          }
        }
        
        
        # If "Log", x ---> log(x)
        if ( p1.bestModels.new$best_model[i] == "Log" ) {
          trans1[i] <- paste(unique_Var[i],"_t=log(", unique_Var[i],")",sep="")
          names(trans1)[i] <- unique_Var[i]
          x1[,ncol(x1)+1] <- log(x1[,unique_Var[i]])
          names(x1)[ncol(x1)] <- paste(unique_Var[i],"_t",sep="")
          # Find coefficients
          for (ii in 1:nrow(p1.bestModels)) {
            if ((p1.allModels[[ii]]$exogenous) == unique_Var[i] & (p1.allModels[[ii]]$endogenous) == unique_Resp[numb_Resp]) {
              Quad_value <- coef(p1.allModels[[ii]]$Log)
            }
          }
          if ( Quad_value[2]>0 ) {
            trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),"+",
                              format(Quad_value[2],scientific=FALSE,digits=2),"*log_",unique_Var[i],sep="")
          } else {
            trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),
                              format(Quad_value[2],scientific=FALSE,digits=2),"*log_",unique_Var[i],sep="")
          }
        }
        
        # If "SqRt", x ---> sqrt(x)
        if ( p1.bestModels.new$best_model[i] == "SqRoot" ) {
          trans1[i] <- paste(unique_Var[i],"_t=sqrt(", unique_Var[i],")",sep="")
          names(trans1)[i] <- unique_Var[i]
          x1[,ncol(x1)+1] <- sqrt(x1[,unique_Var[i]])
          names(x1)[ncol(x1)] <- paste(unique_Var[i],"_t",sep="")
          # Find coefficients
          for (ii in 1:nrow(p1.bestModels)) {
            if ((p1.allModels[[ii]]$exogenous) == unique_Var[i] & (p1.allModels[[ii]]$endogenous) == unique_Resp[numb_Resp]) {
              Quad_value <- coef(p1.allModels[[ii]]$SqRoot)
            }
          }
          if ( Quad_value[2]>0 ) {
            trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),"+",
                              format(Quad_value[2],scientific=FALSE,digits=2),"*sqrt_",unique_Var[i],sep="")
          } else {
            trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),
                              format(Quad_value[2],scientific=FALSE,digits=2),"*sqrt_",unique_Var[i],sep="")
          }
        }
        
        # If "InvSqRt", x ---> invsqrt(x)
        if ( p1.bestModels.new$best_model[i] == "InvSqRoot" ) {
          trans1[i] <- paste(unique_Var[i],"_t=invsqrt(", unique_Var[i],")",sep="")
          names(trans1)[i] <- unique_Var[i]
          x1[,ncol(x1)+1] <- 1/sqrt(x1[,unique_Var[i]])
          names(x1)[ncol(x1)] <- paste(unique_Var[i],"_t",sep="")
          # Find coefficients
          for (ii in 1:nrow(p1.bestModels)) {
            if ((p1.allModels[[ii]]$exogenous) == unique_Var[i] & (p1.allModels[[ii]]$endogenous) == unique_Resp[numb_Resp]) {
              Quad_value <- coef(p1.allModels[[ii]]$InvSqRoot)
            }
          }
          if ( Quad_value[2]>0 ) {
            trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),"+",
                              format(Quad_value[2],scientific=FALSE,digits=2),"*invsqrt_",unique_Var[i],sep="")
          } else {
            trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),
                              format(Quad_value[2],scientific=FALSE,digits=2),"*invsqrt_",unique_Var[i],sep="")
          }
        }
        
        
        # If "Quad", x ---> (x+b/2c)^2
        if ( p1.bestModels.new$best_model[i] == "Quad" ) {
          x1[,ncol(x1)+1] <- x1[,unique_Var[i]]^2
          names(x1)[ncol(x1)] <- paste(unique_Var[i],"__2",sep="")
          
          # Find coefficients
          for (ii in 1:nrow(p1.bestModels)) {
            if ((p1.allModels[[ii]]$exogenous) == unique_Var[i] & (p1.allModels[[ii]]$endogenous) == unique_Resp[numb_Resp]) {
              Quad_value <- coef(p1.allModels[[ii]]$Quad)
            }
          }
          
          trans1[i] <- NA
          names(trans1)[i] <- unique_Var[i]
          if ( (Quad_value[2]>0)&(Quad_value[3]>0) ) {
            trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),"+",
                              format(Quad_value[2],scientific=FALSE,digits=2),'*',unique_Var[i],'+',
                              format(Quad_value[3],scientific=FALSE,digits=2),'*',unique_Var[i],"^2",sep="")
          }
          if ( (Quad_value[2]>0)&(Quad_value[3]<0) ) {
            trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),"+",
                              format(Quad_value[2],scientific=FALSE,digits=2),'*',unique_Var[i],
                              format(Quad_value[3],scientific=FALSE,digits=2),'*',unique_Var[i],"^2",sep="")
          }
          if ( (Quad_value[2]<0)&(Quad_value[3]>0) ) {
            trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),
                              format(Quad_value[2],scientific=FALSE,digits=2),'*',unique_Var[i],'+',
                              format(Quad_value[3],scientific=FALSE,digits=2),'*',unique_Var[i],"^2",sep="")
          }
          if ( (Quad_value[2]<0)&(Quad_value[3]<0) ) {
            trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),
                              format(Quad_value[2],scientific=FALSE,digits=2),'*',unique_Var[i],
                              format(Quad_value[3],scientific=FALSE,digits=2),'*',unique_Var[i],"^2",sep="")
          }
        }
        
        
        # If "CP", x ---> x + pmax(x-tau,0)
        if ( p1.bestModels.new$best_model[i] == "CP" ) {
          
          # Find tau & coefficients
          for (ii in 1:nrow(p1.bestModels)) {
          if ((p1.allModels[[ii]]$exogenous) == unique_Var[i] & (p1.allModels[[ii]]$endogenous) == unique_Resp[numb_Resp]) {
          tau <- (p1.allModels[[ii]]$CP)$psi.history[[3]]
          Quad_value <- coef(p1.allModels[[ii]]$CP)
          }
          }
          
          x1[,ncol(x1)+1] <- pmax(x1[,unique_Var[i]] - tau, 0)
          ## Replace NA with 0 ##
          x1[, ncol(x1)][is.na(x1[, ncol(x1)])] <- 0
          names(x1)[ncol(x1)] <- paste(unique_Var[i],"_c",sep="")
          
          names(trans1)[i] <- unique_Var[i]
          if(tau > 0){
            trans1[i] <- paste(unique_Var[i],"_c=(", unique_Var[i],"-",format(tau,scientific=FALSE,digits=2),")_c",sep=""  )
            if ( (Quad_value[2]>0)&(Quad_value[3]>0) ) {
              trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),"+",
                                format(Quad_value[2],scientific=FALSE,digits=2),'*',unique_Var[i],'+',
                                format(Quad_value[3],scientific=FALSE,digits=2),'*(',unique_Var[i],"-",format(tau,scientific=FALSE,digits=2),")_c",sep="")
            }
            if ( (Quad_value[2]>0)&(Quad_value[3]<0) ) {
              trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),"+",
                                format(Quad_value[2],scientific=FALSE,digits=2),'*',unique_Var[i],
                                format(Quad_value[3],scientific=FALSE,digits=2),'*(',unique_Var[i],"-",format(tau,scientific=FALSE,digits=2),")_c",sep="")
            }
            if ( (Quad_value[2]<0)&(Quad_value[3]>0) ) {
              trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),
                                format(Quad_value[2],scientific=FALSE,digits=2),'*',unique_Var[i],'+',
                                format(Quad_value[3],scientific=FALSE,digits=2),'*(',unique_Var[i],"-",format(tau,scientific=FALSE,digits=2),")_c",sep="")
            }
            if ( (Quad_value[2]<0)&(Quad_value[3]<0) ) {
              trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),
                                format(Quad_value[2],scientific=FALSE,digits=2),'*',unique_Var[i],
                                format(Quad_value[3],scientific=FALSE,digits=2),'*(',unique_Var[i],"-",format(tau,scientific=FALSE,digits=2),")_c",sep="")
            }
          }else{
            trans1[i] <- paste(unique_Var[i],"_c=(", unique_Var[i],"+",format(-tau,scientific=FALSE,digits=2),")_c",sep=""  )
            if ( (Quad_value[2]>0)&(Quad_value[3]>0) ) {
              trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),"+",
                                format(Quad_value[2],scientific=FALSE,digits=2),'*',unique_Var[i],'+',
                                format(Quad_value[3],scientific=FALSE,digits=2),'*(',unique_Var[i],"+",format(-tau,scientific=FALSE,digits=2),")_c",sep="")
            }
            if ( (Quad_value[2]>0)&(Quad_value[3]<0) ) {
              trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),"+",
                                format(Quad_value[2],scientific=FALSE,digits=2),'*',unique_Var[i],
                                format(Quad_value[3],scientific=FALSE,digits=2),'*(',unique_Var[i],"+",format(-tau,scientific=FALSE,digits=2),")_c",sep="")
            }
            if ( (Quad_value[2]<0)&(Quad_value[3]>0) ) {
              trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),
                                format(Quad_value[2],scientific=FALSE,digits=2),'*',unique_Var[i],'+',
                                format(Quad_value[3],scientific=FALSE,digits=2),'*(',unique_Var[i],"+",format(-tau,scientific=FALSE,digits=2),")_c",sep="")
            }
            if ( (Quad_value[2]<0)&(Quad_value[3]<0) ) {
              trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),
                                format(Quad_value[2],scientific=FALSE,digits=2),'*',unique_Var[i],
                                format(Quad_value[3],scientific=FALSE,digits=2),'*(',unique_Var[i],"+",format(tau,scientific=FALSE,digits=2),")_c",sep="")
            }
          }
          
        }
        
        # If "nls", x ---> exp(cx)
        if ( p1.bestModels.new$best_model[i] == "nls" ) {
          # Find coefficients
          for (ii in 1:nrow(p1.bestModels)) {
            if ((p1.allModels[[ii]]$exogenous) == unique_Var[i] & (p1.allModels[[ii]]$endogenous) == unique_Resp[numb_Resp]) {
              Quad_value <- coef(p1.allModels[[ii]]$NLS)
            }
          }
          trans1[i] <- paste(unique_Var[i],"_t=exp(", format(Quad_value[3],scientific=FALSE,digits=2),"*",unique_Var[i],")",sep="")
          names(trans1)[i] <- unique_Var[i]
          x1[,ncol(x1)+1] <- exp(Quad_value[3]*x1[,unique_Var[i]])
          names(x1)[ncol(x1)] <- paste(unique_Var[i],"_t",sep="")
          if ( Quad_value[2]>0 ) {
            trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),"+",
                              format(Quad_value[2],scientific=FALSE,digits=2),"*exp(",
                              format(Quad_value[3],scientific=FALSE,digits=2),unique_Var[i],")",sep="")
          } else {
            trans[i] <- paste(unique_Resp[numb_Resp],"=",format(Quad_value[1],scientific=FALSE,digits=2),
                              format(Quad_value[2],scientific=FALSE,digits=2),"*exp(",
                              format(Quad_value[3],scientific=FALSE,digits=2),unique_Var[i],")",sep="")
          }
        }
        
        
      }
      
      Var_char <- paste(names(x1 %>%
                                dplyr::select(-c(1, unique_Resp[numb_Resp]))), collapse = "+") # all variables except for endogenous
      Rel <- paste0(unique_Resp[numb_Resp], "~", Var_char) # Relation formula
      
      lm.full.model <- (do.call("lm", list(Rel, data = x1 %>% na.omit())))
      lm.best <- MASS::stepAIC(lm.full.model, direction = "backward", trace = FALSE) # variable selection
      lm.best.model[[numb_Resp]] <- lm.best
      R2 <- summary(lm.best)$r.squared
      adj.R2 <- summary(lm.best)$adj.r.squared
      best.vars <- names(coef(lm.best))[-1] 
      
      coeff_number <- round(lm.best$coeff,digits=2)
      names(coeff_number) <- NULL
      coeff_names <- names(lm.best$coeff)
      
      for (n in 1:length(best.vars)){
        Res.print.newrow <- matrix(NA, nrow = 1, ncol = 9)
        Res.print.newrow <- as.data.frame(Res.print.newrow)
        colnames(Res.print.newrow) <- c("endogenous",  "Variable", "GModel",
                                        "GAdj-R2",     "IModel",   "IAdj-R2",
                                        "Gcoeff",       "GModelB", "Model" )
        Res.print.newrow[1,1] <- unique_Resp[numb_Resp]
        
        coeff <- data.frame(coefficients = lm.best$coeff[-1], intercept = lm.best$coeff[1])
        
        # Include tau values in CP
        if (sum(stringr::str_detect(rownames(coeff), "_c")) > 0) {
          for (i in 1:length(trans1)) {
            if (is.na(trans1[i]) == FALSE) {
              original_var <- stringr::str_split(trans1[[i]][1], "=")[[1]][1]
              trans_var <- stringr::str_split(trans1[[i]][1], "=")[[1]][2]
              rownames(coeff) <- stringr::str_replace(rownames(coeff), original_var, trans_var)
            }
          }
        }
        
        lm.best.text <- paste(unique_Resp[numb_Resp],"=",sep="")
        for (i in 1:length(rownames(coeff))) {
          if (is.na(coeff$coefficients[i]) == T) {
            lm.best.text <- lm.best.text
          } else {
          
          if (coeff$coefficients[i] > 0) {
            lm.best.text <- paste0(lm.best.text, "+", signif(coeff$coefficients[i], digits = 3),"*",rownames(coeff)[i])
          }
          if (coeff$coefficients[i] < 0) {
            lm.best.text <- paste0(lm.best.text, signif(coeff$coefficients[i], digits = 3),"*",rownames(coeff)[i])
          } 
          }
          
        }
        
        if (is.na(coeff$intercept[i]) == T) {
          lm.best.text <- lm.best.text
        } else {
        if (coeff$intercept[i] > 0) {
        lm.best.text <- paste0(lm.best.text, "+", signif((coeff$intercept[i]), digits = 3))
        }
        if (coeff$intercept[i] < 0) {
          lm.best.text <- paste0(lm.best.text, signif((coeff$intercept[i]), digits = 3))
        }
        }
        
        Res.print.newrow[1,3] <- lm.best.text
        
        Res.print.newrow[1,4] <- format(adj.R2,scientific=FALSE,digits=2)
        Res.print.newrow[1,7] <- paste(format(lm.best$coeff,digits=2),collapse=" ")
        best.vars.rep <- best.vars
        for (iitemp in 1:length(best.vars)) {
          if (regexpr("_t",best.vars[iitemp])!=-1) {
            temp_str <- trans1[substr(best.vars[iitemp],1,regexpr("_t",best.vars[iitemp])-1)]
            best.vars.rep[iitemp] <- substr(temp_str,regexpr("=",temp_str)+1,nchar(temp_str))
          }
        }
        Res.print.newrow[1,8] <- paste(unique_Resp[numb_Resp],"=",gsub("__2","^2",paste(best.vars.rep,collapse="+")),sep="")
        best.vars.rep2 <- best.vars
        for (j in 1:length(best.vars)) {
          if ((length(grep("__2",best.vars[j]))==0)&(length(grep("_t",best.vars[j]))==0)&(length(grep("_c",best.vars[j]))==0)) {
            best.vars.rep2[j] <- best.vars[j]
          }else{
            if (length(grep("_t",best.vars[j]))!=0) {
              new_string <- substr(best.vars[j],1,regexpr("_t",best.vars[j])[1]-1)
            }
            if (length(grep("__2",best.vars[j]))!=0) {
              new_string <- substr(best.vars[j],1,regexpr("__2",best.vars[j])[1]-1)
            }
            if (length(grep("_c",best.vars[j]))!=0) {
              new_string <- substr(best.vars[j],1,regexpr("_c",best.vars[j])[1]-1)
            }
            best.vars.rep2[j] <- new_string
          }
        }
        Res.print.newrow[1,9] <- paste(unique_Resp[numb_Resp],"=f(",paste(unique(best.vars.rep2),collapse=","),")",sep="")
        if ((length(grep("__2",best.vars[n]))==0)&(length(grep("_t",best.vars[n]))==0)&(length(grep("_c",best.vars[n]))==0)) {
          Res.print.newrow[1,2]<- best.vars[n]
          #Res.print.newrow[1,3] <- lm.best.text.final[numb_Resp]
          Res.print.newrow[1,5] <- trans[best.vars[n]]

          # Individual_summary <- summary(p1.allModels[Var == best.vars[n],Resp == unique_Resp[numb_Resp],p1.bestModels[Var == best.vars[n],Resp == unique_Resp[numb_Resp]]$best_model[1]][[1]])
          # #Res.print.newrow[1,6] <- paste(round(Individual_summary$coeff[,4],digits=3),collapse="  ")
          # if (p1.bestModels[best.vars[n],Resp]!="nls") {
          #   #Res.print.newrow[1,7] <- paste(round(Individual_summary$r.squared,digits=3),
          #   #                               round(Individual_summary$adj.r.squared,digits=3), sep="  ")
          #   Res.print.newrow[1,6] <- round(Individual_summary$adj.r.squared,digits=3)
          # }

          #Res.print.newrow[1,7] <- which(order(summary(lm.best)$coeff[-1,4])==n)
          Res.print <<- rbind(Res.print, Res.print.newrow)
          Res.print.local <- rbind(Res.print.local,Res.print.newrow)
        } else {
          if (length(grep("_t",best.vars[n]))!=0) {
            new_string <- substr(best.vars[n],1,regexpr("_t",best.vars[n])[1]-1)
          }
          if (length(grep("__2",best.vars[n]))!=0) {
            new_string <- substr(best.vars[n],1,regexpr("__2",best.vars[n])[1]-1)
          }
          if (length(grep("_c",best.vars[n]))!=0) {
            new_string <- substr(best.vars[n],1,regexpr("_c",best.vars[n])[1]-1)
          }
          Res.print.newrow[1,2]<- new_string
          Res.print.newrow[1,5] <- trans[new_string]
          
          #Individual_summary <- summary(p1.allModels[new_string,unique_Resp[numb_Resp],p1.bestModels[new_string,unique_Resp[numb_Resp]]][[1]])
          #Res.print.newrow[1,6] <- paste(round(Individual_summary$coeff[,4],digits=3),collapse="  ")
          # if (p1.bestModels[new_string,unique_Resp[numb_Resp]]!="nls") {
          #   #Res.print.newrow[1,7] <- paste(round(Individual_summary$r.squared,digits=3),
          #   #                               round(Individual_summary$adj.r.squared,digits=3), sep="  ")
          #   #Res.print.newrow[1,6] <- round(Individual_summary$adj.r.squared,digits=3)
          # }
          #Res.print.newrow[1,7] <- which(order(summary(lm.best)$coeff[-1,4])==n)
          insert_ind <- 0
          Res.temp <- Res.print[-1,]
          if (nrow(Res.temp)>0) {
            for ( mm in 1:nrow(Res.temp) ) {
              if ( Res.temp[mm,1]==unique_Resp[numb_Resp] ) {
                if ( Res.temp[mm,2]==new_string ) {
                  insert_ind <- 1
                }
              }
            }
          }
          if (insert_ind==0) {
            Res.print <<- rbind(Res.print, Res.print.newrow)
            Res.print.local <- rbind(Res.print.local,Res.print.newrow)
          }
        }
      }

    }
    
    Res.print.local <- Res.print.local %>% 
      dplyr::distinct() %>% 
      dplyr::select(-c("IAdj-R2"))
    
    ################## END: APPLY TRANSFORMATION TO PREDICTORS #################
    
    res <- list(data = x, res.print = Res.print.local, res.full.model = butcher::axe_env(lm.best.model, verbose = FALSE), res.model = lm.full.model)
    class(res) <- c("netSEMp2","list")
    invisible(res)
}
