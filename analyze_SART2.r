expname <- "SART"
rtcol <- 7
successcol <- 6
status_correct <- 1
trimSD <- NA
lowestRT <- NA
highestRT <- NA
followsuccess <- FALSE
conditions <- c(2)
datafilename <- paste("exp_datafiles_", expname, ".txt", sep = "")
d <- read.table(datafilename, as.is = TRUE)

outputMean <- NULL
outputMedian <- NULL
outputMin <- NULL
outputMax <- NULL
outputN <- NULL
outputAccuracy <- NULL
outputErrorCom <- NULL
outputErrorOmi <- NULL
outputSD <- NULL

for (i in 1:length(d[, 1])) {
  if (!is.na(d[i, 1]) & file.info(d[i, 1])$size > 0) {
    
    print(i)
    
    if (exists("x")) {
      rm(x)
    }
    
    if (exists("exclude_lastlines") && exclude_lastlines > 0) {
      o <- system(paste("head -n", exclude_lastlines * -1, d[i, 1]), intern = TRUE)
      x <- read.table(textConnection(o), fill = TRUE)
    }
    
    if (exists("lastlines") && lastlines > 0) {
      o <- system(paste("tail -n", lastlines, d[i, 1]), intern = TRUE)
      x <- read.table(textConnection(o), fill = TRUE)
    }
    
    if (!exists("x")) {
      x <- read.table(d[i, 1], fill = TRUE)
    }
    
    if (exists("include_me"))
      x <- x[x[, 1] == include_me, ]
    if (exists("exclude_me"))
      x <- x[x[, 1] != exclude_me, ]
    
    ntrials <- length(x[, 1])
    if (exists("blockcol")) {
      block <- x[, blockcol]
      
      inclusion <- if (exists("IncludeBlocks")) {
        block %in% IncludeBlocks
      } else {
        rep(TRUE, ntrials)
      }
      
      exclusion <- if (exists("ExcludeBlocks")) {
        block %in% ExcludeBlocks
      } else {
        rep(FALSE, ntrials)
      }
      
      blockselection <- inclusion & !exclusion
    } else {
      blockselection <- rep(TRUE, ntrials)
    }
    
    x <- droplevels(x[blockselection, ])
    
    sink("Rout.txt")
    
    outputdataMean <- numeric()
    outputdataMedian <- numeric()
    outputdataMin <- numeric()
    outputdataMax <- numeric()
    outputdataN <- numeric()
    outputdataAccuracy <- numeric()
    outputdataErrorCom <- numeric()
    outputdataErrorOmi <- numeric()
    outputdataSD <- numeric()
    
    psytkTrim <- function(x, xsd = 3, low = NA, high = NA) {
      y <- x
      if (!is.na(low))
        y <- y[y >= low]
      if (!is.na(high))
        y <- y[y <= high]
      m <- mean(y)
      if (!is.na(xsd)) {
        s <- sd(y) * xsd
        y <- y[y > (m - s) & y < (m + s)]
      }
      return(y)
    }
    
    ntrials <- length(x[, 1])
    trialnum <- 1:ntrials
    
    if (exists("successcol")) {
      CORRECT <- if (exists("status_correct")) {
        x[, successcol] == status_correct
      } else {
        x[, successcol] == 1
      }
    } else {
      CORRECT <- rep(TRUE, ntrials)
    }
    
    if (exists("blockcol")) {
      sameblock <- c(FALSE, x[2:ntrials, blockcol] == x[1:(ntrials - 1), blockcol])
    } else {
      sameblock <- rep(TRUE, ntrials)
    }
    
    allvalues <- function(a) {
      return(a)
    }
    
    if (exists("conditions")) {
      z <- by(trialnum, x[, conditions], allvalues)
      conditionnames <- apply(expand.grid(dimnames(z)), 1, paste, collapse = "_")
    } else {
      z <- list(trialnum)
      conditionnames <- "all_data"
    }
    
    for (i in 1:length(conditionnames)) {
      tc <- conditionnames[i]
      selection <- trialnum %in% z[[i]]
      
      if (followsuccess) {
        selection <- selection & sameblock & c(FALSE, CORRECT[1:(ntrials - 1)])
      }
      
      if (!is.na(lowestRT) | !is.na(highestRT) | !is.na(trimSD)) {
        tmpRtData <- psytkTrim(x[selection & CORRECT & x[, rtcol] < 900, rtcol], 3, lowestRT, highestRT)
      } else {
        tmpRtData <- x[selection & CORRECT & x[, rtcol] < 900, rtcol]
      }
      
      cat("Condition:", tc, "\n")
      cat("-------------------------------------------\n")
      cat("Total trials  :", sum(selection), "trials.\n")
      
      cat("Total correct :", sum(selection & CORRECT), "trials (errors excluded, for RT data below).\n")
      
      if (!is.na(lowestRT) | !is.na(highestRT) | !is.na(trimSD)) {
        cat("After RT trim :", length(tmpRtData), "\n")
      }
      
      cat("Error count   :", sum(!CORRECT & selection), "\n")
      cat("Mean value    :", mean(tmpRtData), "\n")
      cat("Median value  :", median(tmpRtData), "\n")
      cat("Min value     :", min(tmpRtData), "\n")
      cat("Max value     :", max(tmpRtData), "\n")
      cat("Standard Deviation:", sd(tmpRtData), "\n")
      
      total_trials <- sum(selection)
      correct_trials <- sum(selection & CORRECT)
      incorrect_trials <- sum(!CORRECT & selection)
      omission_errors <- sum(selection & !CORRECT & x[, rtcol] >= 900)
      commission_errors <- incorrect_trials - omission_errors
      
      accuracy <- (correct_trials / total_trials) * 100
      error_com <- (commission_errors / total_trials) * 100
      error_omi <- (omission_errors / total_trials) * 100
      
      cat("Accuracy      :", accuracy, "percent\n")
      cat("Error of Commission:", error_com, "percent\n")
      cat("Error of Omission  :", error_omi, "percent\n\n")
      
      outputdataMean <- c(outputdataMean, mean(tmpRtData))
      outputdataMedian <- c(outputdataMedian, median(tmpRtData))
      outputdataMin <- c(outputdataMin, min(tmpRtData))
      outputdataMax <- c(outputdataMax, max(tmpRtData))
      outputdataN <- c(outputdataN, total_trials)
      outputdataAccuracy <- c(outputdataAccuracy, accuracy)
      outputdataErrorCom <- c(outputdataErrorCom, error_com)
      outputdataErrorOmi <- c(outputdataErrorOmi, error_omi)
      outputdataSD <- c(outputdataSD, sd(tmpRtData))
    }
    
    sink()
    
    outputMean <- rbind(outputMean, outputdataMean)
    outputMedian <- rbind(outputMedian, outputdataMedian)
    outputMin <- rbind(outputMin, outputdataMin)
    outputMax <- rbind(outputMax, outputdataMax)
    outputN <- rbind(outputN, outputdataN)
    outputAccuracy <- rbind(outputAccuracy, outputdataAccuracy)
    outputErrorCom <- rbind(outputErrorCom, outputdataErrorCom)
    outputErrorOmi <- rbind(outputErrorOmi, outputdataErrorOmi)
    outputSD <- rbind(outputSD, outputdataSD)
  }
}

print(outputMean)

colnames(outputMean) <- conditionnames
colnames(outputMedian) <- conditionnames
colnames(outputMin) <- conditionnames
colnames(outputMax) <- conditionnames
colnames(outputN) <- conditionnames
colnames(outputAccuracy) <- conditionnames
colnames(outputErrorCom) <- conditionnames
colnames(outputErrorOmi) <- conditionnames
colnames(outputSD) <- conditionnames

outputMean2 <- matrix(ncol=length(conditionnames), nrow=length(d[,1]))
outputMedian2 <- matrix(ncol=length(conditionnames), nrow=length(d[,1]))
outputMin2 <- matrix(ncol=length(conditionnames), nrow=length(d[,1]))
outputMax2 <- matrix(ncol=length(conditionnames), nrow=length(d[,1]))
outputN2 <- matrix(ncol=length(conditionnames), nrow=length(d[,1]))
outputAccuracy2 <- matrix(ncol=length(conditionnames), nrow=length(d[,1]))
outputErrorCom2 <- matrix(ncol=length(conditionnames), nrow=length(d[,1]))
outputErrorOmi2 <- matrix(ncol=length(conditionnames), nrow=length(d[,1]))
outputSD2 <- matrix(ncol=length(conditionnames), nrow=length(d[,1]))

colnames(outputMean2) <- conditionnames
colnames(outputMedian2) <- conditionnames
colnames(outputMin2) <- conditionnames
colnames(outputMax2) <- conditionnames
colnames(outputN2) <- conditionnames
colnames(outputAccuracy2) <- conditionnames
colnames(outputErrorCom2) <- conditionnames
colnames(outputErrorOmi2) <- conditionnames
colnames(outputSD2) <- conditionnames

counter <- 1
for (i in 1:length(d[, 1])) {
  if (is.na(d[i, 1]) | file.info(d[i, 1])$size == 0) {
    outputMean2[i, ] <- rep(NA, length(conditionnames))
    outputMedian2[i, ] <- rep(NA, length(conditionnames))
    outputMin2[i, ] <- rep(NA, length(conditionnames))
    outputMax2[i, ] <- rep(NA, length(conditionnames))
    outputN2[i, ] <- rep(NA, length(conditionnames))
    outputAccuracy2[i, ] <- rep(NA, length(conditionnames))
    outputErrorCom2[i, ] <- rep(NA, length(conditionnames))
    outputErrorOmi2[i, ] <- rep(NA, length(conditionnames))
    outputSD2[i, ] <- rep(NA, length(conditionnames))
  } else {
    outputMean2[i, ] <- outputMean[counter, ]
    outputMedian2[i, ] <- outputMedian[counter, ]
    outputMin2[i, ] <- outputMin[counter, ]
    outputMax2[i, ] <- outputMax[counter, ]
    outputN2[i, ] <- outputN[counter, ]
    outputAccuracy2[i, ] <- outputAccuracy[counter, ]
    outputErrorCom2[i, ] <- outputErrorCom[counter, ]
    outputErrorOmi2[i, ] <- outputErrorOmi[counter, ]
    outputSD2[i, ] <- outputSD[counter, ]
    counter <- counter + 1
  }
}

# You can write these to csv if needed with:
# write.csv(outputMean2, "file_mean.csv", row.names = FALSE)
# write.csv(outputMedian2, "file_median.csv", row.names = FALSE)
# write.csv(outputMin2, "file_min.csv", row.names = FALSE)
# write.csv(outputMax2, "file_max.csv", row.names = FALSE)
# write.csv(outputN2, "file_n.csv", row.names = FALSE)
# write.csv(outputAccuracy2, "file_accuracy.csv", row.names = FALSE)
# write.csv(outputErrorCom2, "file_error_com.csv", row.names = FALSE)
# write.csv(outputErrorOmi2, "file_error_omi.csv", row.names = FALSE)
# write.csv(outputSD2, "file_sd.csv", row.names = FALSE)
