cat("\n Sparky Systems data parser.  ", "Last edited: October 11, 2024.\n\n")

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Parses SparkyJr summary data generated in HPChemStation
# and ouputs fit sample data based on least squares linear
# regression.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Critical problems:
# TC:  The parser is generating an extra column in the "long" output file.
# 
# 
# 
# 
#
# Suggested improvements:
# ALL:  Output std names when reporting residuals and expected values
# ALL:  Add sample type column and identify standards with "STD"
# ALL:  Calculate average AC CO2 conc and subtract from fit sample values
# ALL:  Look-up matching FIDTCD or SRI file based on folder name for ECDFID.txt file[requires testing from S526]
# ALL:  Progress bar
# ALL:  Run log:  Sequence table name MMDDYYC, system (JR, S3, TC), date-time, user, completed
# 
# TC:  Implement auto calibration option for TotalChrom stover incubation data.
# TC:  Eliminate NA column at the end of SumTable
# 

# What's new:

# Added 1% mix standard to CH4 >>  CH4_PID_std[PCNT1] = 10000.0


# ****************************************
# 072017: Attempted to use access network drives using //xxx.xxx.xxx.xxx/, but didn't work.
# 072117: Removed Sparky TC from GC selection.  Data processing code still in script
#         Correct error in plot scales
#         Generate plots in background, without using call to windows()
#         Switch to ggplot()
#         best.lm selected based on NRMSE first, then number of standards used in the calibration.
# 120717: Renames columns including units, replaces negative values with NA, saves results in _dist.csv
#
# 122217: Generates results in ppm and ug.  No longer requires user unit selection.
#         Clears ALL environmental variables at the end of run.  So when run in RStudio, it will delete variables from other scripts.
#
# 081621: best.lm [Line 205] selects (N - 1) unique standards to test for best fit rather than all possible combinations from min.cal to N.
#
# 070124: Based on version 22.11s.  When Sparky Jr was resurrected, FIDTCD were designated System1 and ECD, System2, which is the reverse of what 22.11s
#         expects.  22.11s also expects the sample names to be in ECDFID.txt, which they are not.  24.01s will need to get sample names from FIDTCD.txt.
#
# ****************************************


# To execute type: source("SPKYNetFit20.05s.R", print.eval = TRUE)

# Clear all variables
rm(list = ls())

library(tcltk2)
library(ggplot2)
library(rmarkdown)

# Clear RStudio console
cat("\014")

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Constants

setClass(Class="Regress",
         representation(
           data = "data.frame",
           lm = "lm",
           rsq = "numeric",
           pValue = "numeric",
           cooks.lm = "lm",
           cooks.rsq = "numeric",
           cooks.distance = "numeric",
           cooks.pValue = "numeric",
           cooksKeep = "logical",
           cooksDcutoff = "numeric",
           compound = "character",
           compUnits = "character",
           runOrder = "character",
           var = "character"
         )
)


setClass(Class="mydirlist",
         representation(
           d1="character",
           d2="character",
           dO="character",
           dP="character",
           dR="character",
           dW="character"
         )
)


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Setup input and output directories
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Constants

parserVersion <- "sparkyFit24.10.R"
markdownVersion <- "sparkyFit24.10.Rmd"
myFilters = Filters # Substitute .ps in original Filters array with .csv
myFilters[3,1] = "Excel csv files (*.csv)"
myFilters[3,2] = "*.csv"
rownames(myFilters)[3] = "csv"
runDate <- as.character(Sys.Date())
runTime <- as.character(Sys.time())

nCompounds = 11
# computer = Sys.info()["nodename"]
# min.cal = 2 # Minimum number of unique stds used in calibration regression
radiobuttondone <- tclVar(0)


stdnames = c("HE","5K", "3K", "AMB", "1%", "5%", "10%", "AC", "AIR", "640PPB")
compoundnames = c("_CH4_", "_C2H2_", "_C2H4_", "_C2H6_", "_CO2_", "_N2_", "_N2O_","_O2_")
qaqcnames = c("T", "TR", "TS", "Q", "QA", "QC", "-AMB")
colNames_ug = c("CH4 (ng C)", "C2H2 (ng C)", "C2H4 (ng C)", "C2H6 (ng C)", "CO2 (ug C)", "N2 (%)", "N2O (ng N)", "O2 (%)")
colNames_ppm = c("CH4 (ppm)", "C2H2 (ppm)", "C2H4 (ppm)", "C2H6 (ppm)", "CO2 (ppm)", "N2 (%)", "N2O (ppb)", "O2 (%)")


# fitCompTable =  matrix(nrow = nCompounds, ncol = 3)
# fitCompTable[1,] = c("cal.CO2","Area_CO2", "Height_CO2")
# fitCompTable[2,] = c("cal.N2O","Area_N2O", "Height_N2O")
# fitCompTable[3,] = c("cal.CH4","Area_CH4", "Height_CH4")
# fitCompTable[4,] = c("cal.O2","Area_O2", "Height_O2")
# fitCompTable[5,] = c("cal.N2","Area_N2", "Height_N2")
# fitCompTable[6,] = c("cal.C2H2","Area_C2H2", "Height_C2H2")
# fitCompTable[7,] = c("cal.C2H4","Area_C2H4", "Height_C2H4")
# fitCompTable[8,] = c("cal.C2H6","Area_C2H6", "Height_C2H6")
# fitCompTable[9,] = c("cal.CO2_ECD","Area_CO2_ECD", "Height_CO2_ECD")
# fitCompTable[10,] = c("cal.N2O_ECD","Area_N2O_ECD", "Height_N2O_ECD")
# fitCompTable[11,] = c("cal.CH4_PID","Area_N2O_PID", "Height_N2O_PID")
# 
# colnames(fitCompTable) = c("cal_name", "Area_name", "Height_name")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Initialize and fill standard arrays with values corresponding to stds list based
# on Sample Name.  For now hardcoding might be the easiest way to do this.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Old standard vector definitions

gc = NA
calmode = NA

# concentration to mass conversion factors

ug_C_CO2 = 12.27/5000.0
ng_C_CH4 = 12.27/5.0
ng_N_N2O = 28.63/5000.0
ug_O_O2 = 1308.79/20


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function definitions
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function chooseUnits allows user to choose between concentration and
# mass units for GHG amounts.  ppm are preferred for gas samples; ug for
# insitu incubations.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

chooseCalMode = function(){
  print("chooseCalMode")
  rbVal = as.numeric(tclvalue(rbValue))
  tkdestroy(tt)
  if (rbVal==1){
    # No longer does anything.  Used to display these cutsie messages...      
    # tkmessageBox(title = "Calibration mode:", message="Auto calibration?  Good luck!")
  }
  if (rbVal==2){
    # tkmessageBox(title = "Calibration mode:", message="Manual calibration? Sorry about that.")
  }
  tclvalue(radiobuttondone) <- 1
  calmode <<- rbVal
}# End of function chooseUnits



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function chooseSystems allows user to choose between Sparky Jr. and
# Sparky III
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

chooseSystem = function(){
  print("chooseSystem")
  rbVal = as.numeric(tclvalue(rbValue))
  tkdestroy(tt)
  if (rbVal==1)
    tkmessageBox(title = "Sparky Jr. data processing", message="ECD and FID/TCD data.")
  if (rbVal==2)
    tkmessageBox(title = "Sparky III data processing", message="ECD, FID and PDHID data.")
  tclvalue(radiobuttondone) <- 1
  gc <<- rbVal
}# End of function chooseSystem



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function cooksD.lm calculates linear fit using Cook's distance to remove
# overly influencial points.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# cooksD.lm(calData$conc, calData$height, "height", compName)
# y = calData$conc
# x = calData$area
# var = "area"
# compound = compName

cooksD.lm = function(y, x, var, compound, compUnits, runOrder) {
  
  fitData <- na.omit(data.frame(y = y, x = x))
  cooksKeep <- !is.na(fitData$x)
  names(cooksKeep) <- rownames(fitData)
  
  sampleSize <- nrow(fitData)
  fit.lm <- lm(y ~ x, data = fitData)    
  fit.rsq = summary(fit.lm)$r.squared
  fit.pValue <- summary(fit.lm)$coefficients[4]

  
  if (sampleSize > 2) {
    cooksDcutoff <- 4 / sampleSize
    cooksD <- cooks.distance(fit.lm)
    
    cooksKeep <- (cooksD < cooksDcutoff)
    
    if (sum(!cooksKeep) < 1) {
      
      # No influential data points detected.  Discard the point with the single largest Cook's Distance.
      
      #cooksKeep <- cooksD < max(cooksD)
    }
    cooksData <- data.frame(y = fitData$y[cooksKeep], x = fitData$x[cooksKeep])
    cooks.lm <- lm(y ~ x, data = cooksData)
    cooks.rsq <- summary(cooks.lm)$r.squared
    cooks.pValue <- summary(cooks.lm)$coefficients[4]
    
  } else {
    
    cooks.lm <- fit.lm
    cooks.rsq <- fit.rsq
    cooks.pValue <- fit.pValue
    
  }
  return(new("Regress",
             data = fitData,
             lm = fit.lm,
             rsq = fit.rsq,
             pValue = fit.pValue,
             cooks.lm = cooks.lm,
             cooks.rsq = cooks.rsq,
             cooks.distance = cooksD,
             cooks.pValue = cooks.pValue,
             cooksKeep = c(cooksKeep),
             cooksDcutoff = cooksDcutoff,
             compound = compound,
             compUnits = compUnits,
             runOrder = runOrder,
             var = var)
  )
  
} # End of cooksD.lm


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Function fitSumTable uses calibration curve to fit area or height data
# according to @var slot from bestCal().  Adds results along with upper
# and lower 95% CI.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 


# calibration = bestFit@cooks.lm
# rsq = bestFit@cooks.rsq
# parameter = bestFit@var
# compound = compName
# allData = SumTable


fitSumTable = function(calibration, parameter, compound, allData, unitConversion){
  print(paste("fitSumTable[", compound, "]"))
  nSumCols = ncol(allData)
  insertCol = min(which(names(allData) == paste("Amount_", compound, sep = "")), na.rm = TRUE)
  #  comp.name = which(name_cal == fitCompTable[,"cal_name"])
  if (parameter == "area") {
    x_data = as.numeric(allData[, paste("Area_", compound, sep = "")]) #All Area data
  }else{
    x_data = as.numeric(allData[, paste("Height_", compound, sep = "")]) #All Height data
  }
  fit.results = predict(calibration, data.frame(x = x_data), se.fit = TRUE, level = 0.95, interval = "confidence")
  colnames(fit.results[[1]]) <- paste(colnames(fit.results[[1]]), compound, parameter, sep = "_")
  if (!is.na(unitConversion)) {
    
    #Apply the conversion factor to fit.results[[1]] to get results.converted
    #Generate column names for results.converted
    #cbind fit.results[[1]] and results.converted
    
  }
  
  fitTable = cbind(allData[,1:insertCol], fit.results[[1]], allData[, (insertCol + 1):nSumCols])
  return(fitTable)
} # End of fitSumTable()



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Function setSPKYJRDir() sets up directory assigments for working with
# Sparky Jr data sets.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 


setSPKYJRDir = function(){
  print("setSPKYJRDir")
  
  dir1 <- "C:/"
  dir2 <- "C:/"
  
  
  dirParser <- "C:/Martin/data processing/"
  dirWork <- getwd()
  dirOut <- "C:/Labwork228/SparkyJr/"
  if (!is.na(file.info(dirOut)$isdir)){
    # Has access to Soils 228 lab computer
  } else {
    dirOut <- paste(netAddrS228, "SparkyJr/", sep = "")
    if (!is.na(file.info(dirOut)$isdir)){
      print(paste(dirOut, "available.\n\n"))
      # Has access to Soils 228 lab computer via network
      dirParser <- paste(netAddrS228, parserVersion, "/", sep = "")
      if (is.na(file.info(dirParser)$isdir)){
        dir.create(dirParser, recursive = TRUE)
      }
    } else {
      # Set dirOut to local Labwork folder
      dirOut <- "C:/LabworkLocal/SparkyJr/"
      if (is.na(file.info(dirOut)$isdir)){
        dir.create(dirOut, recursive = TRUE)
      }
      dirParser <- paste("C:/LabworkLocal/", parserVersion, "/", sep = "")
      if (is.na(file.info(dirParser)$isdir)){
        dir.create(dirParser, recursive = TRUE)
      }
    }
  }
  dirReport <- dirOut
  return(new("mydirlist", d1 = dir1, d2 = dir2, dO = dirOut, dP = dirParser, dR = dirReport, dW = dirWork))
}# End of function setSPKYJRDir()


# End of function definitions
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Select GC system for data to be processed
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 


gc = NA
tt <- tktoplevel()
tktitle(tt) <- "Select GC system"
rb1 <- tkradiobutton(tt)
rb2 <- tkradiobutton(tt)
rb3 <- tkradiobutton(tt)
rbValue <- tclVar(1)
tkconfigure(rb1,variable=rbValue,value=1)
tkconfigure(rb2,variable=rbValue,value=2)
tkconfigure(rb3,variable=rbValue,value=3)
tkgrid(tklabel(tt,text="Select GC system:"))
tkgrid(tklabel(tt,text="Sparky Jr. 2024"),rb1)
#tkgrid(tklabel(tt,text="Sparky III"),rb2)
#tkgrid(tklabel(tt,text="Sparky TC"), rb3)
OK.but <- tkbutton(tt, text="OK", command = chooseSystem)
tkgrid(OK.but)
tkfocus(tt)
tkwait.variable(radiobuttondone)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Import Sparky Jr data
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

if (gc == 1){
  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Sparky Jr. specific constants
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  systemName = "SPKYJR"
  TCDList = "CO2N2"
  
  # mydir = setSPKYJRDir()
  
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Select ECD and FID data files to be parsed.
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  #ECD data file contains sample list and must be processed first.
  cat("Imported ECD data from:\n\n")
  
  infile2 = choose.files("C:/*.txt", filters = myFilters[c("txt","All"),],
                         caption = "Choose SPARKY JR. ECD datafile")
  infile2 <- gsub("\\\\", "/", infile2)
  ECDData = read.csv(infile2, sep = " ", header = FALSE, as.is = TRUE, strip.white = TRUE, blank.lines.skip = TRUE)
  cat(infile2, "\n\n")
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Setup output directory based on sequence name.
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  getSequence <- strsplit(readLines(infile2, 4)[3], "\\\\")
  thisSequence <- getSequence[[1]][length(getSequence[[1]])]
  sequenceName <- strsplit(thisSequence, ".S")[[1]][1]
  
# Moved sumReport setup from here to after the FID data is read.  
  

  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Locate first Compound Summary Table in data file and use the column
  # separator character "|" to identify column widths.  Column widths
  # are stored in Seqfwf and used by read.fwf()to extract sequence
  # information: Run, Vial, Inj, Date and Time, Filename.
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  runECD = which(ECDData[, 1] == "Run")
  fw = runECD + 2
  l = ECDData[fw[1],1] 
  m=strsplit(l,"|")
  Seqfw = which(m[[1]][] == "|")
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Manual fix for start of File Name column
  
  Seqfw[4] = Seqfw[4] - 1
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Seqfwf = Seqfw[1]
  for (i in 1:length(Seqfw)-1){
    newSeqfw = Seqfw[i+1] - Seqfw[i]
    Seqfwf = c(Seqfwf, newSeqfw)
  }
  lastSeqfw = length(m[[1]]) - Seqfw[i+1]
  Seqfwf = c(Seqfwf, lastSeqfw)
  
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  #Now re-read ECD data with correct fwf for sequence information
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  
  Sequence = read.fwf(infile2, widths = Seqfwf, as.is = TRUE, strip.white = TRUE, blank.lines.skip = TRUE)
  runSeq = which(Sequence[, 1] == "Run")
  nSamples = runSeq[2] - runSeq[1] - 6
  nSeqCols = ncol(Sequence)
  beginSeq = runSeq[1] + 2
  endSeq = runSeq[2] - 5
  
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Determine fwf info for compound information in ECD datafile:
  # Run, Type, RetTime, Amount, Area, Height, Width, Symm.
  # Store in ECDfwf.
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  l = ECDData[fw[2],1] 
  m=strsplit(l,"|")
  ECDfw = as.numeric(which(m[[1]][] == "|"))
  ECDfwf = ECDfw[1]
  for (i in 1:length(ECDfw)-1){
    newECDfw = ECDfw[i+1]-ECDfw[i]
    ECDfwf = c(ECDfwf, newECDfw)
  }
  lastECDfw = length(m[[1]]) - ECDfw[i+1]
  ECDfwf = c(ECDfwf, lastECDfw)
  
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Re-read ECD datafile using ECDfwf to extract compound information.
  # Overwrite ECDData array with fwf compound data.
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  ECDData = read.fwf(infile2, widths = ECDfwf, as.is = TRUE, strip.white = TRUE, blank.lines.skip = TRUE)
  runECD = which(ECDData[, 1] == "Run")
  endECD = which(ECDData[, 1] == "Mean")
  nECDCols = ncol(ECDData)
  nCompECD = length(endECD)
  beginECD = runECD + 2
  endECD = endECD - 2
  compECD = runECD - 2
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  #Now get fwf for FID data
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  
  # Try to locate FIDData automatically
  infile1 = infile2
  infile1 = sub("/2", "/1", infile1)
  infile1 = sub("/ECDFID.txt", "/FIDTCD.txt", infile1)
  if (file.exists(infile1)){
    cat("Importing FIDTCD data.")
  } else {
    infile1 = choose.files("C:/*.txt", filters = Filters[c("txt","All"),],
                           caption = "Choose SPARKY JR. FID datafile")
    infile1 <- gsub("\\\\", "/", infile1)
  }
  FIDData = read.csv(infile1, sep = " ", header = FALSE, as.is = TRUE, strip.white = TRUE, blank.lines.skip = TRUE)
  cat("Imported FID data from:\n\n")
  cat(infile1, "\n\n")
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Locate first Compound Summary Table in data file and use the column
  # separator character "|" to identify column widths.  Column widths
  # are stored in Seqfwf and used by read.fwf()to extract sequence
  # information: Run, Vial, Inj, Date and Time, Filename.
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  runFID = which(FIDData[, 1] == "Run")
  fw = runFID + 2
  l = FIDData[fw[1],1] 
  m=strsplit(l,"|")
  Seqfw = which(m[[1]][] == "|")
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Manual fix for start of File Name column
  
  Seqfw[4] = Seqfw[4] - 1
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Seqfwf = Seqfw[1]
  for (i in 1:length(Seqfw)-1){
    newSeqfw = Seqfw[i+1] - Seqfw[i]
    Seqfwf = c(Seqfwf, newSeqfw)
  }
  lastSeqfw = length(m[[1]]) - Seqfw[i+1]
  Seqfwf = c(Seqfwf, lastSeqfw)
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  #Now re-read FID data with correct fwf for sequence information
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  
  Sequence = read.fwf(infile1, widths = Seqfwf, as.is = TRUE, strip.white = TRUE, blank.lines.skip = TRUE)
  runSeq = which(Sequence[, 1] == "Run")
  nSamples = runSeq[2] - runSeq[1] - 6
  nSeqCols = ncol(Sequence)
  beginSeq = runSeq[1] + 2
  endSeq = runSeq[2] - 5
  
  
  runFID = which(FIDData[, 1] == "Run")
  fw = runFID + 2
  l = FIDData[fw[2],1] 
  m=strsplit(l,"|")
  FIDfw = as.numeric(which(m[[1]][] == "|")
  )
  FIDfwf = FIDfw[1]
  for (i in 1:length(FIDfw)-1){
    newFIDfw = FIDfw[i+1]-FIDfw[i]
    FIDfwf = c(FIDfwf, newFIDfw)
  }
  lastFIDfw = length(m[[1]]) - FIDfw[i+1]
  FIDfwf = c(FIDfwf, lastFIDfw)
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  #Now re-read FID data with correct fwf for compound information
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  
  FIDData = read.fwf(infile1, widths = FIDfwf, as.is = TRUE, strip.white = TRUE, blank.lines.skip = TRUE)
  runFID = which(FIDData[, 1] == "Run")
  endFID = which(FIDData[, 1] == "Mean")
  nFIDCols = ncol(FIDData)
  nCompFID = length(endFID)
  beginFID = runFID + 2
  endFID = endFID - 2
  compFID = runFID - 2
  nSumCols = nSeqCols + nECDCols * nCompECD + nFIDCols * nCompFID + 1
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  dataHeader = matrix(nrow = 1, ncol = nSumCols)
  SumTable = matrix(nrow = nSamples, ncol = nSumCols)
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  headerRow = min(runSeq)
  dataHeader = Sequence[headerRow, ]
  SumTable = Sequence[beginSeq:endSeq, ]
  
  for (k in 1:nCompECD) {
    startCol = ncol(SumTable)
    for (x in beginECD[k+1]:endECD[k]) {
      i = as.integer(ECDData[x, 1])
      for (j in 1:nECDCols) {
        SumTable[i, (j+startCol)] = ECDData[x, j]
      }
    }
    for (j in 1:nECDCols)
      dataHeader[(j+startCol)] =  paste(ECDData[runECD[k+1],j], strsplit(ECDData[compECD[k+1], 3], " ")[[1]][1], "ECD", sep = "_")
  }
  for (k in 1:nCompFID) {
    startCol = ncol(SumTable)
    for (x in beginFID[k+1]:endFID[k]) {
      i = as.integer(FIDData[x, 1])
      for (j in 1:nFIDCols) {
        SumTable[i, (j+startCol)] = FIDData[x, j]
      }
    }
    for (j in 1:nFIDCols) {
      compName = strsplit(FIDData[compFID[k+1], 3], " ")[[1]][1]
      if (length(grep(compName, TCDList, ignore.case = TRUE)) != 0)
        dataHeader[(j+startCol)] =  paste(FIDData[runFID[k+1],j], compName, sep = "_")
      else
        dataHeader[(j+startCol)] =  paste(FIDData[runFID[k+1],j], compName, sep = "_")
    }
  }
  
  colnames(SumTable) = dataHeader
  
  in1 <- strsplit(infile1, "/")[[1]][length(strsplit(infile1, "/")[[1]])]
  in2 <- strsplit(infile2, "/")[[1]][length(strsplit(infile2, "/")[[1]])]
  
  tkmessageBox(title="SparkyJr: Parsing successful!", message=paste("Parsed files:", in2,"and", in1), icon="info", type="ok")
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  #Set up report directories and output files
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  
  defaultResultsDir <- sub("ECDFID.txt", "", infile2)
  defaultResultsDir <- gsub("/","\\\\", defaultResultsDir)
  
  outfile <- choose.files(paste(defaultResultsDir, "*.csv", sep = ""), caption = "Choose SparkyJr 2024 parser output file")
  tempSplit <- strsplit(outfile, "\\\\")
  outFilename <- tempSplit[[1]][length(tempSplit[[1]])]
  resultsDir <- strsplit(outfile, paste("\\\\", outFilename, sep=""))[[1]][1]
  outfile <- gsub("\\\\", "/", outfile)
  
  
  if(is.na(file.info(resultsDir)$isdir) == FALSE){
    #Output folder exists.  Do nothing
    tkmessageBox(title = "Report folder", message = paste("Results will be written to", outfile), icon="warning", type="ok")
  } else {
    #Output folder does not exist.  Create new folder.
    dir.create(resultsDir)
    tkmessageBox(title = "New report folder", message = paste("Results will be written to", outfile), icon="warning", type="ok")
  }
  
  # Fit summary report filename: sumReport
  
  sumReport <- paste(resultsDir, "\\sparkyjr ", sequenceName, " report ", runDate, ".html", sep = "") 
  
} # End of SparkyJr data import



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Choose whether to perform automated or manual calibration
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

calmode = NA
tt <- tktoplevel()
tktitle(tt) <- "Select Auto or Manual calibration"
rb1 <- tkradiobutton(tt)
rb2 <- tkradiobutton(tt)
rbValue <- tclVar(1)
tkconfigure(rb1,variable=rbValue,value=1)
tkconfigure(rb2,variable=rbValue,value=2)
tkgrid(tklabel(tt,text="Select calibration mode:"))
tkgrid(tklabel(tt,text="Auto"),rb1)
tkgrid(tklabel(tt,text="Only parse data"),rb2)
OK.but <- tkbutton(tt, text="OK", command = chooseCalMode)
tkgrid(OK.but)
tkfocus(tt)
tkwait.variable(radiobuttondone)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Perform auto-calibration of gc data.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if (calmode == 1){
  
  # Calculate calibration curves and fit data to best fit using Cook's distance
  # to remove one or more standards while maintaining maximum range.

  # Identify which columns in SumTable correspond to Area_ and Height_
  
  peakAreaIndex <- which(grepl("Area_", names(SumTable)))
  peakHeightIndex <- which(grepl("Height_", names(SumTable)))
  
  # compoundNames[] contains the names of the compounds that HPCHEM generated results
  # for in the HPCHEM reports ECDFID.txt and FIDTCD.txt.  They will be re-calibrated
  # if there are matching calibration data available in the reports, otherwise only the 
  # HPCHEM generated results will be exported.

  compoundNames <- unlist(strsplit(names(SumTable[peakAreaIndex]), "Area_"))
  compoundNames <- compoundNames[which("" != compoundNames)]
  
  # Calculate the calibration curve for each compound
  
  # Get standard concentrations for the compound.  If they do not exist
  # report "missing concentration table for <compoundName>" and move to next
  # compound.
  
  compoundNamesUpper <- toupper(compoundNames)
  compoundCount <- length(compoundNames)  
  
  sampleNames <- SumTable[,"Sample Name"]
  nSamples <- length(sampleNames)
  stdNames <- c("HE","5K", "3K", "AMB", "1%", "10%", "640PPB", "AC", "AIR", "15PPM")
  qaqcnames <- c("T", "TB", "TBLANK", "TR", "TS", "Q", "QA", "QC", "-AMB")
  
  calibrationNames = c("CH4", "Methane", "ETHANE", "ACETYLE", "ETHYLEN", "CO", "HYDROGE", "CO2", "N2O", "N2", "O2")
  calibrationUnits = c("ppm", "ppm", "ppm", "ppm", "ppm", "ppm", "ppm", "ppm", "ppb", "%", "%")

  CH4stds = c(NA, 5.0, 2.0, 1.0, 10000.0, NA, 2.0, NA, NA, 15)
  Methanestds = c(NA, 5.0, 2.0, 1.0, 10000.0, NA, 2.0, NA, NA, 15)
  ETHANEstds = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 15)
  ACETYLEstds = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 15)
  ETHYLENstds = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 15)
  COstds = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  HYDROGEstds = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)  
  CO2stds = c(NA, 5000.0, 3000.0, 400.0, 10000.0, 100000.0, 400.0, NA, NA, NA)
  N2Ostds = c(NA, 5000.0, 800.0, 200.0, 10000.0, NA, 640.0, NA, NA, NA)
  N2stds = c(NA, 80.0, 90.0, 95.0, 1.0, 90.0, 99.0, NA, NA, NA)
  O2stds = c(NA, 20.0, 10.0, 5.0, 2.0, NA, NA, NA, NA, NA)
  
  conversionFactor <- c(ng_C_CH4, ng_C_CH4, NA, NA, NA, NA, NA, ug_C_CO2, ng_N_N2O, NA, ug_O_O2)
  
  # stdConc is an array of concentrations for compStd with NA for all other samples.

  stdConc = NA
  stdN = length(stdNames)
  cName <- NA
  cNameSubstrCount <- 3
  maxLenCalName <- max(nchar(calibrationNames))
  calNamesUpper <- toupper(calibrationNames)
  loopCounter <- TRUE
  
  for (i in 1:compoundCount){
    stdConc <- NA
    compName <- compoundNames[i]
    cName <- substr(compoundNamesUpper[i], 1, cNameSubstrCount)
    cMatches <- cName == substr(calNamesUpper, 1, cNameSubstrCount)
    calCheck <- sum(cMatches)
    
    # Check for existing calibration standard values
    
    if (calCheck > 0) {
      if (calCheck == 1) {
        
        # Single calibration present
        
      } else {
        
        cat(paste("Matches for", cName, ":", calCheck, "\n\n"))
        
        # More than 1 matching calibrationNames[] and compoundNames[].  At the time of writing
        # this usually happened with Methane, ETHANE and ETHYLENE.  A longer substring should
        # fix the problem.
        
        k <- 1
        while ((calCheck > 1) && k < maxLenCalName) {
          cName <- substr(compoundNamesUpper[i], 1, cNameSubstrCount + k)
          cMatches <- cName == substr(calNamesUpper, 1, cNameSubstrCount + k)
          calCheck <- sum(cMatches)       
          
          # cName <- substr(compoundNamesUpper[i], 1, cNameSubstrCount + k)
          # cMatches <- cName == substr(toupper(calibrationNames), 1, cNameSubstrCount)
          
          #      calCheck <- grep(cName, calibrationNames, ignore.case = TRUE)
          k <- k + 1
        }
      }
    } else {
      
      # No calibration data available for this compound
      
      cat(paste("No calibration data available for ", compoundNames[i], ".\n\n"))
    }
    
    cat(paste("calCheck : ", calCheck, cName, compoundNames[i], "\n\n", sep = " "))
    cat(paste(i, cName, "matches:", calCheck, compName, "\n\n"))
    
    # Generate array of standard values, Height and Area for column compoundNames[i]
    calName <- calibrationNames[cMatches]
    calUnits <- calibrationUnits[cMatches]
    
    for (std in 1:stdN){
      stdIndex <- grep(stdNames[std], sampleNames, ignore.case = TRUE) 
      if (length(stdIndex) > 0)
        
        # There is at least 1 sample that contains the standard keyword stdnames[std]
        stdConc[stdIndex] <- get(paste(calName, "stds", sep = ""))[std]
    }
    
    area <- as.numeric(SumTable[,paste("Area_", compName, sep = "")])
    height <- as.numeric(SumTable[,paste("Height_", compName, sep = "")])
    calData <- data.frame(area, height, conc = stdConc)
    calData <- na.omit(calData)
    runOrder <- rownames(calData)
    
    
    if (nrow(calData) < 2 | length(unique(calData$conc)) < 2) {
      
      tkmessageBox(title = "No calibration generated", message = paste("Insufficient std data for", compName), icon="warning", type="ok")
      
      # Output results to HTML document
      
    } else {
      
      # Generate calibration for compound column i

      height.lm <- cooksD.lm(calData$conc, calData$height, "height", compName, calUnits, runOrder)
      area.lm <- cooksD.lm(calData$conc, calData$area, "area", compName, calUnits, runOrder)
      
      lmOptions <- c(area.lm, height.lm)
      
      lmChoice <- c(area.lm@cooks.rsq, height.lm@cooks.rsq)
      lmChoice <- which.max(lmChoice) 
      if (length(lmChoice) > 0) {
        bestFit <- lmOptions[[lmChoice]] 
      } else {
        bestFit <- lmOptions[[1]] # area.lm
      }

      SumTable <- fitSumTable(bestFit@cooks.lm, bestFit@var, compName, SumTable, conversionFactor[cMatches])
      conversionFactor <- NA

      # Generate RMarkdown results for compound compName
      
      outMarkdownFile <-  paste0(resultsDir, "\\", compName, " fit ", runDate, ".html")
      
      # This assumes that sparkyFit2410.Rmd is in the same directory as sparkyFit24.10.R and that
      # Session: Working Directory is set to Source File Location
      
      # tkmessageBox(title = "Starting Location", message = getwd(), icon="warning", type="ok")
      # tkmessageBox(title = "Markdown file exists?", message = file.exists('sparkyFit24.10.Rmd'), icon="warning", type="ok")
      
      rmarkdown::render(input = 'sparkyFit24.10.Rmd',
                        output_file = outMarkdownFile,
                        output_format = "html_document",
                        params = list(sequenceName = sequenceName,
                                      compName = compName,
                                      bestFit = bestFit,
                                      runDate = runDate,
                                      loopCounter = loopCounter,
                                      parserVersion = parserVersion,
                                      markdownVersion = markdownVersion,
                                      outfile = outfile
                        )
      )
      
      if (loopCounter) {
        loopCounter <- FALSE
      }
      
      newtxt <- readLines(outMarkdownFile, n = -1)
      fileConn <- file(sumReport, "a")
      writeLines(newtxt, fileConn)
      #write(newtxt, file = fileConn, append = TRUE)
      close(fileConn)
    }
  }
} #End of data processing

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Manual calibration selected. Output raw data only, to be processed by hand.
# Now write SumTable to .csv file
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

write.csv(SumTable, file = outfile, quote = FALSE, row.names = FALSE)

# Open HTML summary file

browseURL(sumReport)  

# Cleanup environment

rm(list = ls())

# To generate a list of the variables generated by a script do the following:

# 1) noquote(ls())
# 2) Copy/Paste output from Console
# 3) enclose text in double quotes and assign to variable a
# 4) Execute the following commands in order:
#       b <- gsub("\\[.+?\\]", "", a)  # Remove square bracket line numbering
#       c <- gsub("\\s+", " ", b)      # Remove extra spaces between variable names
#       d <- gsub(" ", ",", c)         # Replace spaces with comma
# 5) rm(...list of variables in d...)
