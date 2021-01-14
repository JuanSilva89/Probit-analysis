
### Script version 4
# Script implemented using 
# RStudio Version 1.2.5033 
# R version 3.6.3 (2020-02-29)


# UPDATES:
# This script now includes LC1 and LC99 in the output matrix
# If input data contains 1 strain only, then the output is the summary of the probit analysis (text file)
# If input data contains >1 strain, then the output is a log-probit line (pdf) BUT NO summary of the probit analysis

### INSTALL "MASS" PACKAGE FOR THE FIRST TIME ONLY!!
install.packages("MASS")


### SELECT ALL THE SCRIPT BELOW TO RUN YOUR DATA
entryf<-file.choose()

dataf<-(read.table(entryf,header=T,stringsAsFactors=T)) # Added stringsAsFactors = T
mort<-ifelse(dataf$dead/dataf$total==0,0.0006,ifelse(dataf$dead/dataf$total==1,1-0.0006,dataf$dead/dataf$total)) #calculate the mortality in each dos and replicate from the data
dataf<-cbind(dataf,mort) #add the mortality in the dataframe

#cont<-0 # if cont = 0 controls are not used, if cont =1 they are used
if (min(dataf$dose)==0) {
  datac<-dataf[dataf$dose==0,]
  dataf<-dataf[dataf$dose>0,]
  cont<-0} # check wether there are controls and create to files, one with the controls one with rest of the data
dataf.output<-dataf
if (max(datac$dead)==0){cont<-0}
if (min(datac$dead)>0.00) {cont<-1}  # correction is applied to all treatment doses regardless of mortality in the controls
if (cont==1) {  #if mortality in the controls is higher than 0%, we use the Abbott's correction 
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
require(MASS)
myprobit <- glm((dataf$mort) ~ log10(dataf$dose), family =quasibinomial(link = probit))

logLC1 <- dose.p(myprobit, cf = 1:2, p = 0.01)
lc1 <- 10^c(logLC1 + c(0) * attr(logLC1, "SE"))
lc1up <-10^c(logLC1 + c(1.96) * attr(logLC1, "SE"))
lc1low <-10^c(logLC1 + c(-1.96) * attr(logLC1, "SE"))

logLC5 <- dose.p(myprobit, cf = 1:2, p = 0.05)
lc5 <- 10^c(logLC5 + c(0) * attr(logLC5, "SE"))
lc5up <-10^c(logLC5 + c(1.96) * attr(logLC5, "SE"))
lc5low <-10^c(logLC5 + c(-1.96) * attr(logLC5, "SE"))

logLC10 <- dose.p(myprobit, cf = 1:2, p = 0.1)
lc10 <- 10^c(logLC10 + c(0) * attr(logLC5, "SE"))
lc10up <-10^c(logLC10 + c(1.96) * attr(logLC10, "SE"))
lc10low <-10^c(logLC10 + c(-1.96) * attr(logLC10, "SE"))

logLC50 <- dose.p(myprobit, cf = 1:2, p = 0.5)
lc50 <- 10^c(logLC50 + c(0) * attr(logLC50, "SE"))
lc50up <-10^c(logLC50 + c(1.96) * attr(logLC50, "SE"))
lc50low <-10^c(logLC50 + c(-1.96) * attr(logLC50, "SE"))

logLC90 <- dose.p(myprobit, cf = 1:2, p = 0.9)
lc90 <- 10^c(logLC90 + c(0) * attr(logLC90, "SE"))
lc90up <-10^c(logLC90 + c(1.96) * attr(logLC90, "SE"))
lc90low <-10^c(logLC90 + c(-1.96) * attr(logLC90, "SE"))

logLC95 <- dose.p(myprobit, cf = 1:2, p = 0.95)
lc95 <- 10^c(logLC95 + c(0) * attr(logLC95, "SE"))
lc95up <-10^c(logLC95 + c(1.96) * attr(logLC95, "SE"))
lc95low <-10^c(logLC95 + c(-1.96) * attr(logLC95, "SE"))
lc95up

logLC99 <- dose.p(myprobit, cf = 1:2, p = 0.99)
lc99 <- 10^c(logLC99 + c(0) * attr(logLC99, "SE"))
lc99up <-10^c(logLC99 + c(1.96) * attr(logLC99, "SE"))
lc99low <-10^c(logLC99 + c(-1.96) * attr(logLC99, "SE"))
lc99up

summary<-summary(myprobit)
summary
intercept<-5+myprobit$coefficients[1]## Intercept of the line equation
intercept<-round(intercept,digit=3)
slope<-summary$coefficients[2] ## Slope of the line equation
slope<-round(slope,digit=3)

X<-log10(dataf$dose)
Y<-(intercept + slope*X)
new.expected.mortality<-pnorm(Y-5)/pnorm(dataf$probit-5)*dataf.output$dead

r<-dataf.output$dead #dead observed
nP<-new.expected.mortality #dead expected
P<-pnorm(Y-5) #dataf$mort #proportion of dead observed  
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
higher.control.mort<-ifelse(dataf.output$mort<=datac$mort,"INVALID POINT: MORTALITY LOWER OR EQUAL THAN CONTROL","")

# Create a matrix with input data (same data as input tab delimited text)
original.data<-matrix(0,nrow =length(dataf$dose) ,ncol = 5)
original.data<-cbind(dataf.output$dose,dataf.output$dead,dataf.output$total,round(dataf.output$mort*100,digits = 3),higher.control.mort)
colnames(original.data)<-c("Dose","Dead","Total","Observed Mortality (%)","")
original.data

# Create a matrix with corrected mortality, probit, expected mortality, observed mortality and Chi-square 
output<-matrix(0,nrow = length(dataf$dose),ncol = 9)
output<-cbind(dataf$dose,raw.abbott.correction,"",round(dataf$probit,digits = 3),dataf$total,dataf.output$dead,round(new.expected.mortality,digits = 3),"",round(chitt.output,digits = 3),extreme.values)
colnames(output)<-c("Dose","Mort.Corr(%)","Probit","Total","Dead","Expected","X2","","","")
output
#This is the mortality corrected after optimization

# Create matrix to output the intercept, slope and slope standard error
sortie<-matrix(0,nrow=1,ncol=3)
sortie<-cbind(intercept,slope,summary$coefficients[4])
colnames(sortie)<-c("Interc.","Slope","SlopeSE")
sortie

# Show the equation of the line
eq.mat<-matrix(0,nrow=1,ncol=1)
eq.mat<-cbind(eq.line)
colnames(eq.mat)<-c("Equation of the line")
eq.mat

# Add warning if data are represented by a line or not
if(chis.sq.round<0.05){warn.message<-c("DATA NOT REPRESENTED BY A LINE")}
if(chis.sq.round>0.05){warn.message<-c("DATA REPRESENTED BY A LINE")}
warn.message

chimatrix<-matrix(0,nrow=1,ncol=2)
chimatrix<-cbind(chi.data,warn.message)
chimatrix

lcs<-matrix(0,nrow=7,ncol=4)
LCcol<-c(lc1,lc5,lc10,lc50,lc90,lc95,lc99)
LCmin<-c(lc1low,lc5low,lc10low,lc50low,lc90low,lc95low,lc99low)
LCmax<-c(lc1up,lc5up,lc10up,lc50up,lc90up,lc95up,lc99up)
number<-c(1,5,10,50,90,95,99)
lcs<-cbind(number,LCcol,LCmin,LCmax)
colnames(lcs)<-c("LC","Value","Min","Max")
lcs


strains<-levels(dataf$strain) # check what are the strains being tested

### Export probit analysis results in a tab delimited text
### If more than one strain included in input file then DO NOT create text file

if(length(strains)==1) {
  exporttab<-gsub(".txt","=Results.txt",entryf)
  write("                   SCOTT LAB - PROBIT ANALYSIS", exporttab)
  write("", exporttab, append=T)
  write("Script version: v4", exporttab, append=T)
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
  write.table(round(sortie,digits = 3),append = T ,exporttab, sep = "\t", col.names=T,row.names = F,quote=F)
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
mortality_fit<-seq(0.0006,1-0.0006,0.01)# generates a set of mortalities
probitfit<-qnorm(mortality_fit)+5 # transforms mortalities into probits
pred_logLC <- dose.p(model, cf = 1:2, p = mortality_fit) # predicted log-probit line
pred_lc <- 10^c(pred_logLC + c(0) * attr(pred_logLC, "SE")) # predicted LC
pred_lc_CI_up <-10^c(pred_logLC + c(1.96) * attr(pred_logLC, "SE")) # predicted upper 95% CI
pred_lc_CI_low <-10^c(pred_logLC + c(-1.96) * attr(pred_logLC, "SE")) # predicted lower 95% CI
dose_ci<-cbind(probitfit,pred_lc_CI_low,pred_lc_CI_up) # predicted LC and 95% CI to be plotted
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
### LD50 (CI 95%):  4.845	(4.496 - 5.222)

# FINNEY:
# Slope (Standard Error): 4.176 (0.466) 
# Equation line:  Y = 2.134 + 4.176 * X
# X2 weight:      1.671
# LD50 (CI 95%):  4.85 (4.38 - 5.36) 


