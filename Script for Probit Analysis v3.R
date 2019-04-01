### INSTALL "MASS" PACKAGE FOR THE FIRST TIME ONLY!!
install.packages("MASS")





### SELECT ALL THE SCRIPT BELOW TO RUN YOUR DATA
entryf<-file.choose()

dataf<-(read.table(entryf,header=T))

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
  souch<-levels(dataf$strain) #creates a vector with the name of the strains
  Cc<-1:length(souch) #creates a vector for the estimations of mortality in the controls, C
  Conv<-1:length(souch) #creates a vector for the convergence of the estimations of mortality in the controls
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

probmort<-sapply(dataf$mort,qnorm) # apply probit transformation to the data
dataf<-cbind(dataf,probmort) # add these probit-transformed mortality in the dataframe

attach(dataf)
probit<-(5+dataf$probmort)
dataf<-cbind(dataf,probit) #add these probit values to dataframe
dataf
datac

 # Abbott's correction (without optimization) to include in matrix output

require(MASS)
myprobit <- glm((dataf$mort) ~ log10(dataf$dose), family =quasibinomial(link = probit))

logLC50 <- dose.p(myprobit, cf = 1:2, p = 0.5)
lc50 <- 10^c(logLC50 + c(0) * attr(logLC50, "SE"))
lc50up <-10^c(logLC50 + c(1.96) * attr(logLC50, "SE"))
lc50low <-10^c(logLC50 + c(-1.96) * attr(logLC50, "SE"))

logLC5 <- dose.p(myprobit, cf = 1:2, p = 0.05)
lc5 <- 10^c(logLC5 + c(0) * attr(logLC5, "SE"))
lc5up <-10^c(logLC5 + c(1.96) * attr(logLC5, "SE"))
lc5low <-10^c(logLC5 + c(-1.96) * attr(logLC5, "SE"))

logLC10 <- dose.p(myprobit, cf = 1:2, p = 0.1)
lc10 <- 10^c(logLC10 + c(0) * attr(logLC5, "SE"))
lc10up <-10^c(logLC10 + c(1.96) * attr(logLC10, "SE"))
lc10low <-10^c(logLC10 + c(-1.96) * attr(logLC10, "SE"))

logLC90 <- dose.p(myprobit, cf = 1:2, p = 0.9)
lc90 <- 10^c(logLC90 + c(0) * attr(logLC90, "SE"))
lc90up <-10^c(logLC90 + c(1.96) * attr(logLC90, "SE"))
lc90low <-10^c(logLC90 + c(-1.96) * attr(logLC90, "SE"))

logLC95 <- dose.p(myprobit, cf = 1:2, p = 0.95)
lc95 <- 10^c(logLC95 + c(0) * attr(logLC95, "SE"))
lc95up <-10^c(logLC95 + c(1.96) * attr(logLC95, "SE"))
lc95low <-10^c(logLC95 + c(-1.96) * attr(logLC95, "SE"))
lc95up
summary<-summary(myprobit)
summary
intercept<-5+myprobit$coefficients[1]## Intercept for the line equation
intercept<-round(intercept,digit=3)
slope<-summary$coefficients[2] ## Slope for the line equation
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


# store results in an array
names.strain.insecticide<-c("Strain:        ","Insecticide:","Date:        ")
byinsecticide<-levels(dataf$insecticide)
bystrain<-levels(dataf$strain)
bydate<-levels(dataf$date)
byins.bystrain.bydate<-c(bystrain,byinsecticide,bydate)
main.info<-matrix(0,nrow = 3,ncol = 2)
main.info<-cbind(names.strain.insecticide,byins.bystrain.bydate)
main.info


if(datac$mort*100>10){warn.control.mort<-c("WARNING: CONTROL MORTALITY EXCEEDING 10%")}
if(datac$mort*100<10){warn.control.mort<-c(" ")}
warn.control.mort

control.mort<-matrix(0,nrow=1,ncol=2)
cmort.round<-round((datac$dead/datac$total)*100,digits = 3)
control.mort<-paste("Control Mortality (%) =",cmort.round,"(",datac$dead,"/",datac$total,")",warn.control.mort)
control.mort

raw.abbott.correction<-round((dataf.output$mort-datac$mort)/(1-datac$mort)*100,digits=3) #Abbott correction without optimization to include in output file
dataf$mort # Abbott's corrected percent mortality after optimization is applied (NOT INCLUDED IN THE OUTPUT FILE) 

#Warning message: No controls added
if (min(data_raw$dose)>0) {control.mort<-c("*** PLEASE ADD CONTROLS ***")}
if (min(data_raw$dose)==0) {control.mort<-matrix(0,nrow=1,ncol=2)
cmort.round<-round((datac$dead/datac$total)*100,digits = 3)
control.mort<-paste("Control Mortality (%) =",cmort.round,"(",datac$dead,"/",datac$total,")",warn.control.mort)
}
control.mort

dataf$mort
datac$mort

extreme.values<-ifelse(new.expected.mortality<=5,"*",ifelse(new.expected.mortality/dataf.output$total*100>=95,"*"," "))

higher.control.mort<-ifelse(dataf.output$mort<=datac$mort,"INVALID POINT: MORTALITY LOWER OR EQUAL THAN CONTROL","")


original.data<-matrix(0,nrow =length(dataf$dose) ,ncol = 5)
original.data<-cbind(dataf.output$dose,dataf.output$dead,dataf.output$total,round(dataf.output$mort*100,digits = 3),higher.control.mort)
colnames(original.data)<-c("Dose","Dead","Total","Observed Mortality (%)","")
original.data

output<-matrix(0,nrow = length(dataf$dose),ncol = 9)
output<-cbind(dataf$dose,raw.abbott.correction,"",round(dataf$probit,digits = 3),dataf$total,dataf.output$dead,round(new.expected.mortality,digits = 3),"",round(chitt.output,digits = 3),extreme.values)
colnames(output)<-c("Dose","Mort.Corr(%)","Probit","Total","Dead","Expected","X2","","","")
output
#This is the mortality corrected after optimization

sortie<-matrix(0,nrow=1,ncol=3)
sortie<-cbind(intercept,slope,summary$coefficients[4])
colnames(sortie)<-c("Interc.","Slope","SlopeSE")
sortie

eq.mat<-matrix(0,nrow=1,ncol=1)
eq.mat<-cbind(eq.line)
colnames(eq.mat)<-c("Equation of the line")
eq.mat

if(chis.sq.round<0.05){warn.message<-c("DATA NOT REPRESENTED BY A LINE")}
if(chis.sq.round>0.05){warn.message<-c("DATA REPRESENTED BY A LINE")}
warn.message

chimatrix<-matrix(0,nrow=1,ncol=2)
chimatrix<-cbind(chi.data,warn.message)
chimatrix

lcs<-matrix(0,nrow=5,ncol=4)
LCcol<-c(lc5,lc10,lc50,lc90,lc95)
LCmin<-c(lc5low,lc10low,lc50low,lc90low,lc95low)
LCmax<-c(lc5up,lc10up,lc50up,lc90up,lc95up)
number<-c(5,10,50,90,95)
lcs<-cbind(number,LCcol,LCmin,LCmax)
colnames(lcs)<-c("LC","Value","Min","Max")
lcs

exporttab<-gsub(".txt","=Results.txt",entryf)
write("                   SCOTT LAB - PROBIT ANALYSIS", exporttab)
write("", exporttab, append=T)
write("Script version: v3", exporttab, append=T)
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

remove(list=c("datac"))
