
### Script version 6
# Script implemented using 
# RStudio Version 1.2.5033 
# R version 3.6.3 (2021-08-18)


# UPDATES:
# Probit analysis is aborted when any dose has a percent mortality equal or lower than controls or 100% mortality  

# This script now includes LC1 and LC99 in the output matrix

# If input data contains 1 strain only, then the output is the summary of the probit analysis (text file)

# If input data contains >1 strain, then the output is a log-probit line (pdf) BUT NO summary of the probit analysis

# Current script does not depend on the MASS package to calculate lethal doses/concentrations. Instead, it calculates
# the different parameters (calculated in line 92 through 148) used by Finney (weigthing coefficients, Z, NWX, NWY, etc).
# Although the previous parameters are not displayed in the output file of the probit analysis, the user
# can inspect the parameters as stored in a matrix (please see lines 138 to 141).

### SELECT ALL THE SCRIPT BELOW TO RUN YOUR DATA
entryf<-file.choose()

dataf<-(read.table(entryf,header=T,stringsAsFactors=T)) # Added stringsAsFactors = T
######## ############## This was the original mortality correction when 0% and 100% mortality were obtained  mort<-ifelse(dataf$dead/dataf$total==0,0.0006,ifelse(dataf$dead/dataf$total==1,1-0.0006,dataf$dead/dataf$total)) #calculate the mortality in each dose and replicate from the data
mort<-dataf$dead/dataf$total
dataf<-cbind(dataf,mort) #add the mortality in the dataframe

#cont<-0 # if cont = 0 controls are not used
if (min(dataf$dose)==0) {
  datac<-dataf[dataf$dose==0,]
  dataf<-dataf[dataf$dose>0,]
  cont<-0} # create two data frames (one with the controls only and another with rest of the data)
dataf.output<-dataf
if (max(datac$dead)==0){cont<-0}
if (min(datac$dead)>0.00) {cont<-1}  # correction is applied to all treatment doses regardless of mortality in the controls
if (cont==1) {  # Apply Abbott's correction 
  datafc<-dataf # creates a dataframe identical to dataf
  strains<-levels(dataf$strain) #creates a vector with the name of the strains
  Cc<-1:length(strains) #creates a vector for the estimations of mortality in the controls, C
  Conv<-1:length(strains) #creates a vector for the convergence of the estimations of mortality in the controls
  for (i in 1:length(souch)) { # the procedure is done for each strain in turn 
    ii<-datafc$strain == souch[i]
    if (min(datafc$mort[ii])>0) { # estimates the mortality in the controls only if there are some dead individuals in the controls, otherwise C=0
      probit_C <- function(Cx){ 	# creates a function which applies the Abbott's correction to the data for a given C               
        datafc$mort[ii]<-(dataf$mort[ii]-Cx)/(1-Cx)
        datafc$dead[ii]<-datafc$mort[ii]*dataf$total[ii]
        
        y<-cbind(datafc$dead[ii],datafc$total[ii]-datafc$dead[ii])
        moda<-glm(y~log10(datafc$dose[ii]),family = quasibinomial(link=probit))                       
        L<--moda$deviance/2
        j<-datac$strain==souch[i]
        L<-L+sum(datac$dead[j]*log10(Cx)+(datac$total[j]-datac$dead[j])*log10(1-Cx))
        return(L)} # estimates conjointly the mortality in the controls (C) and the coefficients of the linear regression for the strain considered 
      
      param_low <- 0.000000000001
      param_up  <- min(datafc$mort[ii])
      param_ini<-NULL
      param_ini <- runif(1, min=param_low, max=param_up)
      result_opt <- optim(param_ini,probit_C,control=list(fnscale=-1,trace=1),method="L-BFGS-B",lower=param_low,upper=param_up) # Optimizing function to get the best estimates for C
      Cc[i]<-result_opt$par # gives the best C
      Conv[i]<-result_opt$convergence # gives the convergence of the optimizing function for the best C
      dataf$mort[ii]<-(dataf$mort[ii]-Cc[i])/(1-Cc[i]) # calculates the frequency of dead individuals in dataf taking the best C into account, using Abbott's formula
      dataf$dead[ii]<-dataf$mort[ii]*dataf$total[ii]	# updates the number of dead individuals in dataf
    }
    else Cc[i]<-0
  }
  fitC<-cbind(souch,Cc,Conv)
  colnames(fitC)<-c("Strain","MortalityControls","Convergence(OK if 0)") # creates a dataframe giving the results of the estimations of C, to be printed for the user in the results file
}

probmort<-sapply(dataf$mort,qnorm) # apply probit transformation to the probability of mortality
dataf<-cbind(dataf,probmort) # add these probit-transformed mortality in the dataframe

attach(dataf)
probit<-(5+dataf$probmort)
dataf<-cbind(dataf,probit) #add these probit values to dataframe


# Abbott's correction (without optimization) to include in matrix output
# Calculation of the LC1, 5, 10, 50, 90, 95 and 99 and their corresponding 95% confidence intervals

myprobit <- glm((dataf$mort) ~ log10(dataf$dose), family =quasibinomial(link = probit))
summary<-summary(myprobit)
summary
intercept<-5+myprobit$coefficients[1]## Intercept of the line equation
intercept<-round(intercept,digit=3)
slope<-summary$coefficients[2] ## Slope of the line equation
slope<-round(slope,digit=3)

X<-log10(dataf$dose)
Y<-(intercept + slope*X)
new.expected.mortality<-pnorm(Y-5)*dataf.output$total

Pi<-c(3.141593) # Create a vector with the value of pi. Pi is a constant included in the equation to calculate Z. Thus, the standard value used in R for this constant is 3.141593.
Y<-(intercept + slope*X)
P <- pnorm(Y-5)
Q <- 1-P

# Create a vector with total number of individuals tested per dose
N<-dataf$total


# Finney's calculation of Z (Equation 3.11 from P.33)
Z<-(1/sqrt(2*Pi))*exp(-0.5*(Y-5)^2)

# Working probit
Y0<-c(Y-P/Z)
Y0

# weighting is 1/Z according to Finney
k<-c(1/Z)/100
k

# calculation of weighting coefficient (w) - (Equation 3.10 from P.33)
Wt.coeff<-(Z^2)/(P*Q)
Wt.coeff

# calculation of weight to assigned dose/concentration (NW)
NW<-Wt.coeff*dataf$total

# calculation of weighting coefficient (NWX)
NWX<-NW*X
NWX
sum(NWX)

# calculation of parameter (NWY)
WY<-W*Y

NWY<-NW*Y
sum(NWY)
# calculation of parameter (NWXY)
NWXY<-NW*X*Y
sum(NWXY)

# Create a matrix with probit parameters. This matrix is meant to resemble Table 4.1 from Finney's book, P. 62.
probit.parameters<-matrix(0,nrow =length(dataf$dose) ,ncol = 12)
probit.parameters<-cbind(X,Y,Z,Y0,k,Wt.coeff,NW,NWX,NWY,NWXY,NW*(X^2),(NW*X)^2)
colnames(probit.parameters)<-c("X","Y","Z","Y0","k","Wt.coeff","NW","NWX","NWY","NWXY","NW(X^2)","(NWX)^2")
probit.parameters

## We now proceed to calculate the the parameters a and b from the linear regression (Y = aX + b)
# Calculation of the slope (constant b)
b<-(sum(NW)*sum(NWXY)-sum(NWX)*sum(NWY))/(sum(NW)*sum(NW*(X^2))-(sum(NW*X)^2))

# Calculation of the intersection to x axis (constant a)
a<- ((sum(NWY)/sum(NW))-(b*(sum(NWX)/sum(NW))))

# Calculation of X.hat. This parameter is used to estimate the variance for each dose (e.g. see formula on line 169).
X.hat<-sum(NWX)/sum(NW)

# Convert 1% mortality into probit (e.g. 50% mortality equals to a probit of 5)
probit_1<-round(qnorm(0.01)+5,digits = 3)

logLC1<-((probit_1-a)/b)
lc1<-10^c((probit_1-a)/b)
# Alternatively, the LD1 can be calculated using the function:


############# Calculation of LDs for output file #################


# Calculating LD5
probit_1<-round(qnorm(0.01)+5,digits = 3) # Convert percent mortality into probit (e.g. 90% mortality is equivalent to a probit value of 7.326)
logLC1<-((probit_1-a)/b)
lc1<-10^c((probit_1-a)/b)
Variance.logLC1 <- 1/b^2 * { 1/(sum(NW)) + ((logLC1-X.hat)^2 /(sum(NW*(X^2))-sum(NWX)^2/sum(NW))) } # Variance calculated as in equation 3.12 from Finney's book
Sm.LC1<- sqrt(Variance.logLC1) # Standard error calculated using equation 3.13 from Finney's book
lc1up <-10^c(logLC1 + c(1.96) * Sm.LC1) # Value of the upper 95% confidence interval
lc1low <-10^c(logLC1 + c(-1.96) * Sm.LC1) # Value of the lower 95% confidence interval

# Calculating LD5
probit_5<-round(qnorm(0.05)+5,digits = 3) # Convert percent mortality into probit (e.g. 90% mortality is equivalent to a probit value of 7.326)
logLC5<-((probit_5-a)/b)
lc5<-10^c((probit_5-a)/b)
Variance.logLC5 <- 1/b^2 * { 1/(sum(NW)) + ((logLC5-X.hat)^2 /(sum(NW*(X^2))-sum(NWX)^2/sum(NW))) } # Variance calculated as in equation 3.12 from Finney's book
Sm.LC5<- sqrt(Variance.logLC5) # Standard error calculated using equation 3.13 from Finney's book
lc5up <-10^c(logLC5 + c(1.96) * Sm.LC5) # Value of the upper 95% confidence interval
lc5low <-10^c(logLC5 + c(-1.96) * Sm.LC5) # Value of the lower 95% confidence interval

# Calculating LD10
probit_10<-round(qnorm(0.10)+5,digits = 3) # Convert percent mortality into probit (e.g. 90% mortality is equivalent to a probit value of 7.326)
logLC10<-((probit_10-a)/b)
lc10<-10^c((probit_10-a)/b)
Variance.logLC10 <- 1/b^2 * { 1/(sum(NW)) + ((logLC10-X.hat)^2 /(sum(NW*(X^2))-sum(NWX)^2/sum(NW))) } # Variance calculated as in equation 3.12 from Finney's book
Sm.LC10<- sqrt(Variance.logLC10) # Standard error calculated using equation 3.13 from Finney's book
lc10up <-10^c(logLC10 + c(1.96) * Sm.LC10) # Value of the upper 95% confidence interval
lc10low <-10^c(logLC10 + c(-1.96) * Sm.LC10) # Value of the lower 95% confidence interval

# Calculating LD50
probit_50<-round(qnorm(0.50)+5,digits = 3) # Convert percent mortality into probit (e.g. 90% mortality is equivalent to a probit value of 7.326)
logLC50<-((probit_50-a)/b)
lc50<-10^c((probit_50-a)/b)
Variance.logLC50 <- 1/b^2 * { 1/(sum(NW)) + ((logLC50-X.hat)^2 /(sum(NW*(X^2))-sum(NWX)^2/sum(NW))) } # Variance calculated as in equation 3.12 from Finney's book
Sm.LC50<- sqrt(Variance.logLC50) # Standard error calculated using equation 3.13 from Finney's book
lc50up <-10^c(logLC50 + c(1.96) * Sm.LC50) # Value of the upper 95% confidence interval
lc50low <-10^c(logLC50 + c(-1.96) * Sm.LC50) # Value of the lower 95% confidence interval

# Calculating LD90
probit_90<-round(qnorm(0.90)+5,digits = 3) # Convert percent mortality into probit (e.g. 90% mortality is equivalent to a probit value of 7.326)
logLC90<-((probit_90-a)/b)
lc90<-10^c((probit_90-a)/b)
Variance.logLC90 <- 1/b^2 * { 1/(sum(NW)) + ((logLC90-X.hat)^2 /(sum(NW*(X^2))-sum(NWX)^2/sum(NW))) } # Variance calculated as in equation 3.12 from Finney's book
Sm.LC90<- sqrt(Variance.logLC90) # Standard error calculated using equation 3.13 from Finney's book
lc90up <-10^c(logLC90 + c(1.96) * Sm.LC90) # Value of the upper 95% confidence interval
lc90low <-10^c(logLC90 + c(-1.96) * Sm.LC90) # Value of the lower 95% confidence interval

# Calculating LD95 
probit_95<-round(qnorm(0.95)+5,digits = 3) # Convert percent mortality into probit (e.g. 90% mortality is equivalent to a probit value of 7.326)
logLC95<-((probit_95-a)/b)
lc95<-10^c((probit_95-a)/b)
Variance.logLC95 <- 1/b^2 * { 1/(sum(NW)) + ((logLC95-X.hat)^2 /(sum(NW*(X^2))-sum(NWX)^2/sum(NW))) } # Variance calculated as in equation 3.12 from Finney's book
Sm.LC95<- sqrt(Variance.logLC95) # Standard error calculated using equation 3.13 from Finney's book
lc95up <-10^c(logLC95 + c(1.96) * Sm.LC95) # Value of the upper 95% confidence interval
lc95low <-10^c(logLC95 + c(-1.96) * Sm.LC95) # Value of the lower 95% confidence interval

# Calculating LD99 
probit_99<-round(qnorm(0.99)+5,digits = 3) # Convert percent mortality into probit (e.g. 90% mortality is equivalent to a probit value of 7.326)
logLC99<-((probit_99-a)/b)
lc99<-10^c((probit_99-a)/b)
Variance.logLC99 <- 1/b^2 * { 1/(sum(NW)) + ((logLC99-X.hat)^2 /(sum(NW*(X^2))-sum(NWX)^2/sum(NW))) } # Variance calculated as in equation 3.12 from Finney's book
Sm.LC99<- sqrt(Variance.logLC99) # Standard error calculated using equation 3.13 from Finney's book
# Upper 95% confidence intervals are calculated as LC +/- 1.96*Standard error
lc99up <-10^c(logLC99 + c(1.96) * Sm.LC99) # Value of the upper 95% confidence interval
lc99low <-10^c(logLC99 + c(-1.96) * Sm.LC99) # Value of the lower 95% confidence interval

r<-dataf.output$dead #dead observed
nP<-new.expected.mortality #dead expected
P<-pnorm(Y-5) #proportion of dead observed  
up<-r-nP
down<-nP*(1-P)
X2<-up^2/down
chitt<-(sum(X2)) #chi-square test between the observed dead numbers (data) and the dead numbers predicted by the regression (mods$fitted.values*total)
chitt<-round(chitt,digits = 3)
chitt.output<-(X2)
chi.sq<-pchisq(q=chitt,df=length(dataf$dead)-2,lower.tail=F)
chis.sq.round<-round(chi.sq,digits = 3)

chi.data<-paste("p(X2 =",chitt,",df = ",df=length(dataf$dead)-2,")=",chis.sq.round)
chi.data

extreme.values<-ifelse(new.expected.mortality<=5,"*",ifelse(new.expected.mortality/dataf.output$total*100>=95,"*"," "))

eq.line<-paste('Y =',intercept,'+',slope, '* X')
eq.line


# store results in array
names.strain.insecticide<-c("Strain:        ","Insecticide:","Date:        ")
byinsecticide<-levels(dataf$insecticide)
bystrain<-levels(dataf$strain)
bydate<-levels(dataf$date)
byins.bystrain.bydate<-c(bystrain,byinsecticide,bydate)
main.info<-matrix(0,nrow = 3,ncol = 2)
main.info<-cbind(names.strain.insecticide,byins.bystrain.bydate)
main.info

# Add a warning if control mortality exceeds 10%
if(datac$mort*100>10){warn.control.mort<-c("WARNING: CONTROL MORTALITY EXCEEDING 10%")}
if(datac$mort*100<10){warn.control.mort<-c(" ")}
warn.control.mort

control.mort<-matrix(0,nrow=1,ncol=2)
cmort.round<-round((datac$dead/datac$total)*100,digits = 3)
control.mort<-paste("Control Mortality (%) =",cmort.round,"(",datac$dead,"/",datac$total,")",warn.control.mort)
control.mort

raw.abbott.correction<-round((dataf.output$mort-datac$mort)/(1-datac$mort)*100,digits=3) #Abbott correction without optimization to include in output file
dataf$mort # Abbott's corrected percent mortality after optimization is applied (NOT INCLUDED IN THE OUTPUT FILE) 

# Warning message: No controls added
if (min(data_raw$dose)>0) {control.mort<-c("*** PLEASE ADD CONTROLS ***")}
if (min(data_raw$dose)==0) {control.mort<-matrix(0,nrow=1,ncol=2)
cmort.round<-round((datac$dead/datac$total)*100,digits = 3)
control.mort<-paste("Control Mortality (%) =",cmort.round,"(",datac$dead,"/",datac$total,")",warn.control.mort)
}

# Warning message: Extreme points expected mortality lower than 5% or higher than 95%
extreme.values<-ifelse(new.expected.mortality<=5,"*",ifelse(new.expected.mortality/dataf.output$total*100>=95,"*"," "))

# Warning message: Dose produced lower or equal mortality relative to control
lower.mort.than.control<-ifelse(dataf.output$mort<=datac$mort,"INVALID POINT: MORTALITY LOWER OR EQUAL THAN CONTROL","")
warning.100.mort<-ifelse(dataf.output$mort==1,"INVALID POINT: DOSE HAS 100% MORTALITY","")


# Create a matrix with input data (same data as input tab delimited text)
original.data<-matrix(0,nrow =length(dataf$dose) ,ncol = 6)
original.data<-cbind(dataf.output$dose,dataf.output$dead,dataf.output$total,round(dataf.output$mort*100,digits = 3),lower.mort.than.control,warning.100.mort)
colnames(original.data)<-c("Dose","Dead","Total","Observed Mortality (%)","","")
original.data

# Use a conditional, if there is 0%, 100% or % mortality lower than the control then the output will create a warning message BUT NO results will be provided
# Create a matrix with corrected mortality, probit, expected mortality, observed mortality and Chi-square 

if(min(dataf.output$mort)<=datac$mort) {
  output<-c(" ")

} else if (max(dataf.output$mort)==1) {
  output<-c(" ")

} else {  
  output<-matrix(0,nrow = length(dataf$dose),ncol = 9)
  output<-cbind(dataf$dose,raw.abbott.correction,"",round(dataf$probit,digits = 3),dataf$total,dataf.output$dead,round(new.expected.mortality,digits = 3),"",round(chitt.output,digits = 3),extreme.values)
  colnames(output)<-c("Dose","Mort.Corr(%)","Probit","Total","Dead","Expected","X2","","","")
}

#This is the mortality corrected after optimization

# Create matrix to output the intercept, slope and slope standard error

if(min(dataf.output$mort)<=datac$mort) {
  probit.line.summary<-c(" ")
} else if (max(dataf.output$mort)==1) {
  probit.line.summary<-c(" ") 
} else if (min(dataf.output$mort)==0) {
  probit.line.summary<-c(" ") 
  
} else {  
  probit.line.summary<-matrix(0,nrow=1,ncol=3)
  probit.line.summary<-cbind(intercept,slope,summary$coefficients[4])
  colnames(probit.line.summary)<-c("Interc.","Slope","SlopeSE")
  probit.line.summary
}

# Show the equation of the line
if(min(dataf.output$mort)<=datac$mort) {
  eq.mat<-c(" ")
} else if (max(dataf.output$mort)==1) {
  eq.mat<-c(" ") 
} else if (min(dataf.output$mort)==0) {
  eq.mat<-c(" ") 
  
} else {  
  eq.mat<-matrix(0,nrow=1,ncol=1)
  eq.mat<-cbind(eq.line)
  colnames(eq.mat)<-c("Equation of the line")
  eq.mat
}


# Add warning if data are represented by a line or not
if(min(dataf.output$mort)<=datac$mort) {
  warn.message<-c(" ")
} else if (max(dataf.output$mort)==1) {
  warn.message<-c(" ") 
} else if (min(dataf.output$mort)==0) {
  warn.message<-c(" ") 
  
} else {  
  if(chis.sq.round<0.05){warn.message<-c("DATA NOT REPRESENTED BY A LINE")}
  if(chis.sq.round>0.05){warn.message<-c("DATA REPRESENTED BY A LINE")}
  warn.message
}

if(min(dataf.output$mort)<=datac$mort) {
  chimatrix<-c(" ")
} else if (max(dataf.output$mort)==1) {
  chimatrix<-c(" ") 
} else if (min(dataf.output$mort)==0) {
  chimatrix<-c(" ") 
  
} else {  
  chimatrix<-matrix(0,nrow=1,ncol=2)
  chimatrix<-cbind(chi.data,warn.message)
  chimatrix
}


if(min(dataf.output$mort)<=datac$mort) {
  lcs<-c(" ")
} else if (max(dataf.output$mort)==1) {
  lcs<-c(" ") 
} else if (min(dataf.output$mort)==0) {
  lcs<-c(" ") 
  
} else {  
  lcs<-matrix(0,nrow=7,ncol=4)
  LCcol<-c(lc1,lc5,lc10,lc50,lc90,lc95,lc99)
  LCmin<-c(lc1low,lc5low,lc10low,lc50low,lc90low,lc95low,lc99low)
  LCmax<-c(lc1up,lc5up,lc10up,lc50up,lc90up,lc95up,lc99up)
  number<-c(1,5,10,50,90,95,99)
  lcs<-cbind(number,LCcol,LCmin,LCmax)
  colnames(lcs)<-c("LC","Value","Min","Max")
  lcs
}



strains<-levels(dataf$strain) # check what are the strains being tested

### Export probit analysis results in a tab delimited text
### If more than one strain included in input file then DO NOT create text file

if(length(strains)==1) {
  exporttab<-gsub(".txt","=Results.txt",entryf)
  write("                   SCOTT LAB - PROBIT ANALYSIS", exporttab)
  write("", exporttab, append=T)
  write("Script version: v5", exporttab, append=T)
  write.table(main.info, exporttab, sep = "\t", col.names=F,row.names = F,quote=F,append = T)
  write("", exporttab, append=T)
  write("", exporttab, append=T)
  write.table(control.mort,append = T ,exporttab, sep = "\t", col.names=F,row.names = F,quote=F)
  write("", exporttab, append=T)
  write("Original Data", exporttab, append=T)
  write("", exporttab, append=T)
  write.table(original.data,append = T ,exporttab, sep = "\t", col.names=T,row.names = F,quote=F)
  write("", exporttab, append=T)
  write("", exporttab, append=T)
  write.table(round(probit.line.summary,digits = 3),append = T ,exporttab, sep = "\t", col.names=T,row.names = F,quote=F)
  write("", exporttab, append =T)
  write.table(format(eq.mat,digits = 3),append = T, exporttab, sep = "\t", col.names=T,row.names = F,quote=F)
  write("", exporttab, append=T)
  write.table(format(chimatrix,digits = 3),append = T, exporttab, sep = "\t", col.names=F,row.names = F,quote=F)
  write("", exporttab, append=T)
  write("        Abbott                                  Dead               ", exporttab, append=T)
  write.table(output,append = T, exporttab, sep = "\t", col.names=T,row.names = F,quote=F)
  write("", exporttab, append=T)
  write("Extreme values are marked with an asterisk (*)", exporttab, append=T)
  write("", exporttab, append=T)
  write("", exporttab, append=T)
  write("                Conf. Interv. 95%", exporttab, append=T)
  write.table(round(lcs,digits = 3),append = T, exporttab, sep = "\t", col.names=T,row.names = F,quote=F)
 
} else {
  no.export<-vector() # If there are more than one strains then no probit output (.txt) is created
} 

###################################### Commands to plot the log-probit line #########################################

### Export probit plots (pdf)
### If only one strain included in input file then DO NOT create plot file

if(length(strains)>1) {
  
# assign a color to the different strains
if (length(strains)>10) palette(rainbow(length(strains))) else if (length(strains)>1) palette(rainbow(10)) 
if (is.null(dataf$color)==T) {
  color<-1:length(strain)
  dataf<-cbind(dataf,color)
  for (i in 1:length(strains)) {
    ii<-strain == strains[i]
    dataf$color[ii]<-i} # assign the color arbitrarily if the user did not provide any indication
}
colstrains<-1:length(strains)
for (i in 1:length(strains)) {
  ii<-strain == strains[i]
  colstrains[i]<-dataf$color[ii][1]}


if (is.null(dataf$symbol)==T) {
  symbol<-1:length(strain)
  dataf<-cbind(dataf,symbol)
  for (i in 1:length(strains)) {
    ii<-strain == strains[i]
    dataf$symbol[ii]<-i} # assign the data point symbol arbitrarily if the user did not provide any indication
}
symstrains<-1:length(strains)
for (i in 1:length(strains)) {
  ii<-strain == strains[i]
  symstrains[i]<-dataf$symbol[ii][1]}

dmin<-floor(log10(min(dose))) # determines the graph upper and lower log doses, from the data provided
dmax<-ceiling(log10(max(dose))) 
dose_min<- 10^(dmin) # lower log dose transformed to dose
dose_max<- 10^(dmax) # upper log dose transformed to dose
pmort_min<- qnorm(0.006)+5 # probability of minimal dose transformed to probit
pmort_max<- qnorm(1-0.006)+5 # probability of maximal dose transformed to probit

CIplot<-function(model){sum<-summary(model) #function to calculate the CI of the regressions for each strain
a<-sum$coefficient[1]+5    # model intercept
b<-sum$coefficient[2]      # model slope
minldose<-((pmort_min-a)/b) # predicted log(dose) for 0% mortality
maxldose<-((pmort_max-a)/b) # predicted log(dose) for 100% mortality
datalfit<-seq(minldose-0.2,maxldose+0.2,by=0.0153) # generates a set of log doses
datafit<-10^datalfit# transform log doses into doses
mortality_fit<-seq(0.0006,1-0.0006,0.01)# generates a set of mortality
probitfit<-qnorm(mortality_fit)+5 # transforms mortality into probits
pred_logLC <- ((probitfit-a)/b) # predicted log-probit line
pred_LC<-10^pred_logLC
pred.Variance.logLC <- 1/b^2 * { 1/(sum(NW)) + ((pred_logLC-X.hat)^2 /(sum(NW*(X^2))-sum(NWX)^2/sum(NW))) } # Variance calculated as in equation 3.12 from Finney's book
pred_Sm.LC<- sqrt(pred.Variance.logLC) # Standard error calculated using equation 3.13 from Finney's book
pred_lc_CI_up <-10^c(pred_logLC + c(1.96) * pred_Sm.LC)
pred_lc_CI_low <-10^c(pred_logLC + c(-1.96) * pred_Sm.LC)
dose_ci<-cbind(probitfit,pred_lc_CI_low,pred_lc_CI_up) # creat a matrix with predicted probit (X-axis) and 95% CI (Y-axis) to be plotted
return(dose_ci)
}

# build the graph in a pdf file with confidence intervals
exportgraph<-gsub(".txt","=PLOT.pdf",entryf) # create the pdf file from the original data file
pdf(exportgraph)

ii<-strain == strains[1] #plots the data points for the first strain
plot((dose[ii]),probit[ii],log = "x",xlim=c(dose_min,dose_max),ylim=c(floor(pmort_min),ceiling(pmort_max)),ylab="mortalit?",yaxt="n",xaxt="n", ann=FALSE ,col=colstrains[1],pch=symstrains[1])
for (i in 2:length(strains)) { 
  ii<-strain == strains[i]
  points(dose[ii],probit[ii],col=colstrains[i],pch=symstrains[i]) }  #plots the data points for the other strains 
for (i in 1:length(strains)) { # Compute the probit regression for each strain
 locdata<-dataf[strain == strains[i],]
 myprobit_per_strain <- glm(mort ~ log10(dose),data = locdata, family =quasibinomial(link = probit))
 CIfit<-CIplot(myprobit_per_strain) #plots the CI of the regressions for each strain
  lines(CIfit[,2],CIfit[,1],type="l", lty=3, col=colstrains[i])  # plots the lower 95% CI
  lines(CIfit[,3],CIfit[,1],type="l", lty=3, col=colstrains[i])# plots the upper 95% CI
  myprobit_per_strain$coefficients[1] <-myprobit_per_strain$coefficients[1]+5 # Add 5 to the intercept
abline(myprobit_per_strain, col=colstrains[i])  #plots the regressions for each strain
}


# graph format
#labels the y axis in probit scale
labely<-c(1,5,seq(10,90,10),95,99) 
axis(2, at=5+qnorm(labely/100),labels=labely,las=2, adj=0) # converts percent mortality into probit 
mtext("Mortality (%)", side=2, line=3)


#labels the x axis 
for (i in dmin:dmax) axis(1,at=10 ^i,labels=substitute(10^k,list(k=i))) 
axis.at <- 10 ^c(dmin:dmax)
axis(1, at = 2:9 * rep(axis.at[-1]/10, each = 8),
     tcl = -0.5, labels = FALSE)                   
mtext(expression( Log ( dose ("ng per individual")) ), side=1, line=3) 

# creates a legend for the graph
legend(dose_min, 7.5, strains, col = colstrains, pch=symstrains, cex = 0.8)

} else {
no.export.graph<-vector() # If there is one strain only then no graph (.pdf) is created 
}

remove(list=c("datac"))

dev.off()

##############################################################################################

### To assess the quality of this script there is a comparison between the results
### obtained with this script and David Finney's book using the same dataset (Example 1 from the book: Probit Analysis.1917, 3rd edition. Cambridge University Press)

### Example 1 results obtained by:

# THIS SCRIPT:
### Slope (Standard Error): 4.217	(0.366)
### Equation line:  Y = 2.11 + 4.217 * X
### X2 weight:      1.731
### LD50 (CI 95%):  4.845	(4.386 - 5.352)

# FINNEY:
# Slope (Standard Error): 4.176 (0.466) 
# Equation line:  Y = 2.134 + 4.176 * X
# X2 weight:      1.671
# LD50 (CI 95%):  4.85 (4.38 - 5.36) 


