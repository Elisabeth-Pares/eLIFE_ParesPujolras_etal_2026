###################################################################################
## Behavioural analysis code for: 
# "Perceptual glimpses are locally accumulated and globally maintained at distinct processing levels" 
# Includes Eq. 1 to Eq. 4 in the paper. 

# Author: Elisabeth Pares-Pujolras 
###################################################################################

source("C:/Users/elisa/Desktop/utils/elisR.R") # Load standard utils

#Initialise required packages
{
  require(dplyr)
  require(tidyr)
  require(lme4)
  require(ggplot2)
  require(cowplot)
  require(rstatix)
  require(formattable)
  require(ggpattern)
  require(sjPlot)
}

#Define data folder and create output directories				
{				
  dir <-dirname(rstudioapi::getSourceEditorContext()$path) #Get current directory	
  dir.data <- paste(dir,'/DATA/', sep = '')				
  dir.stats <- paste(dir,'/STATS/', sep = '')				
  dir.create(dir.stats)				
  dir.figures <- paste(dir,'/FIGURES/', sep = '')				
  dir.create(dir.figures)				
}		

#Plotting parameters
{
  lowHighColors = c('red2', 'dodgerblue2') #c('mediumpurple1', 'mediumpurple4')
  lowMidHighColors = c('red2', 'green2', 'dodgerblue')
  fullColors = c('red2', 'green2', 'green4', 'dodgerblue2')#rev(c('mediumpurple4','mediumpurple3','mediumpurple2','mediumpurple1'))#rev(grDevices::heat.colors(9))# colorspace::sequential_hcl( h = 126,
  greenCols = c('green4', 'green2')
  
  mycols <- colors()[c(76, 53,33,36)] #Saturated
  mycols <- colors()[c(411, 425,573,571)] #Pastel
  mycols <- colors()[c(57,58,59,60)] #Pastel
  
  mycols <- c('#ffa600', '#ff6361', '#bc5090', '#58508d')
 
  # New facet label names for supp variable
  supp.labs <- c("Experiment 1", "Experiment 2")
  names(supp.labs) <- c(1, 2)
}

# Load csv file with all data
#Preprocess
for (experiment in (c(1,2))){
  if (experiment == 1){ setwd('C:/Users/elisa/Desktop/GapsTask/Experiment 1/Data/behavioural data')}
  else if (experiment == 2 ){setwd('C:/Users/elisa/Desktop/GapsTask/Experiment 2/Data/behavioural data')}
  
  allfiles = list.files(pattern = "*.csv")
  unwanted <- list.files(pattern = "*procBeh")
  temp <- base::setdiff(allfiles, unwanted)
  
  myfiles = lapply(temp, read.csv)
  myData <- do.call(rbind,myfiles)
  
  myData <- filter(myData, npulses != 0) # Eliminates first row of data 
  
  
  myData$npulses <- as.factor(myData$npulses)
  
  myData$meanCoh <- as.factor((myData$p1coh + myData$p2coh)/2)
  myData[myData$npulses ==1,]$meanCoh <- myData[myData$npulses == 1,]$p1coh
  
  myData$meanCoh <- as.factor(myData$meanCoh)
  myData$gapdur <- as.factor(myData$gapdur)
  
  #Add block ID to all multiples of 70; all blocks were 70 trials long except for block 14 in Exp. 2
  myData$block <- 0
  myData[myData$trial <=70,]$block <- 1
  myData[myData$trial >70 & myData$trial<=70*2,]$block <- 2
  myData[myData$trial >70*2 & myData$trial<=70*3,]$block <- 3
  myData[myData$trial >70*3 & myData$trial<=70*4,]$block <- 4
  myData[myData$trial >70*4 & myData$trial<=70*5,]$block <- 5
  myData[myData$trial >70*5 & myData$trial<=70*6,]$block <- 6
  myData[myData$trial >70*6 & myData$trial<=70*7,]$block <- 7
  myData[myData$trial >70*7 & myData$trial<=70*8,]$block <- 8
  myData[myData$trial >70*8 & myData$trial<=70*9,]$block <- 9
  myData[myData$trial >70*9 & myData$trial<=70*10,]$block <-10
  
  if (experiment ==2) {
    myData[myData$trial >70*10 & myData$trial<=70*11,]$block <-11
    myData[myData$trial >70*11 & myData$trial<=70*12,]$block <-12
    myData[myData$trial >70*12 & myData$trial<=70*13,]$block <-13
    myData[myData$trial >70*13 & myData$trial<=70*14,]$block <-14
  }
  
  
  
  #Add label for high/low for each block 
  myData <- myData%>%
    group_by(ID,block)%>%
    mutate(thisLow = min(p1coh,p2coh)) %>%
    mutate(thisHigh = max(p1coh,p2coh))%>%
    arrange(ID, block)
  
  #Get mean high/low coherence values for each participant to use in Matlab for sample-wise regression 
  sumData <- myData%>%
    group_by(ID)%>%
    summarise(meanLow = mean(thisLow),
              meanHigh = mean(thisHigh))
  
  
  myData$cohLevel <- 'NA'
  myData[myData$meanCoh == myData$thisLow,]$cohLevel <- '0Low'
  myData[myData$meanCoh == myData$thisHigh,]$cohLevel <- '2High'
  myData[myData$npulses != 1 &  myData$meanCoh != myData$thisLow & myData$meanCoh != myData$thisHigh,]$cohLevel <- '1Med'
  
  
  myData$p1cohLev <- 'NA'
  myData[myData$p1coh == myData$thisLow,]$p1cohLev <- 'low_'
  myData[myData$p1coh == myData$thisHigh,]$p1cohLev <- 'high_'
  
  myData$p2cohLev <- 'NA'
  myData[myData$p2coh == myData$thisLow,]$p2cohLev <- 'low'
  myData[myData$p2coh == myData$thisHigh,]$p2cohLev <- 'high'
  
  myData$bothCoh <- 'NA'
  myData$bothCoh <- paste(myData$p1cohLev, myData$p2cohLev)
  
  myData$allCoh <- 'NA'
  myData[myData$bothCoh == 'high_ high',]$allCoh <- '6.HH'
  myData[myData$bothCoh == 'high_ low',]$allCoh <- '5.HL'
  myData[myData$bothCoh == 'low_ high',]$allCoh <- '3.LH'
  myData[myData$bothCoh == 'low_ low',]$allCoh <- '2.LL'
  
  myData[myData$p1cohLev == 'low_' & myData$npulses == 1,]$allCoh <- '1.L-'
  myData[myData$p1cohLev == 'high_' & myData$npulses == 1,]$allCoh <- '4.H-'
  
  
  #Binary response code 
  myData$resp <- NA
  myData[myData$pdir == 0,]$resp <- 1 #Right
  myData[myData$pdir == -180,]$resp <- 0 #Left
  
  #Signed response  
  myData$sign_resp <- NA
  myData[myData$pdir == 0,]$sign_resp <- 1 #Right
  myData[myData$pdir == -180,]$sign_resp <- -1 #Left
  
  #Add experiment number
  myData$exp <- experiment
  
  demoData <- myData%>%
    group_by(ID, exp)%>%
    summarise_all(.funs = c('first')) 
  
  
  #Scale coherences to use as numeric 
  myData$p1coh_s <- myData$p1coh/1000
  myData$p2coh_s <- myData$p2coh/1000
  
  myData[myData$npulses != 1 &  myData$meanCoh != myData$thisLow & myData$meanCoh != myData$thisHigh,]$cohLevel <- '1Med'
  
  #Average within participant
  gaData <- myData %>%
    group_by(ID, npulses, gapdur,cohLevel) %>%
    summarise(meanAcc = mean(accuracy), 
              meanRT = mean(RT))
  
  #Get only one row per trial
  sumData <- myData %>%
    group_by(ID, trial)%>%
    summarise_all(.funs = c('first')) 
  
  
  #Average within participant
  ga_sumData <- sumData %>%
    group_by(ID, block) %>%
    summarise(meanAcc = mean(accuracy), 
              meanCohs = mean(as.numeric(meanCoh)))
  
  #Condition codes are useful to specify location on plots 
  sumData$condCode  <- 0
  sumData[sumData$bothCoh == 'low_ low',]$condCode <- 3
  sumData[sumData$bothCoh == 'low_ high',]$condCode <- 4
  sumData[sumData$bothCoh == 'high_ low',]$condCode <- 5
  sumData[sumData$bothCoh == 'high_ high',]$condCode <- 6
  
  sumData[sumData$p1cohLev == 'low_' & sumData$npulses == 1,]$condCode <- 1
  sumData[sumData$p1cohLev == 'high_' & sumData$npulses == 1,]$condCode <- 2
  
  if (experiment == 1){e1data <- myData; e1sumData <- sumData; e1_gasumData <- ga_sumData; e1demo <- demoData; write.csv(e1sumData, 'C:/Users/elisa/Desktop/GapsTask/Experiment 1/Data/behavioural data/EP1_procBehData.csv')}
  else if (experiment == 2 ){e2data <- myData; e2sumData <- sumData; e2_gasumData <- ga_sumData; e2demo <- demoData;  write.csv(e2sumData, 'C:/Users/elisa/Desktop/GapsTask/Experiment 2/Data/behavioural data/EP2_procBehData.csv')}
  
}

#FIGURE S1- including excluded participants in Exp. 2
{
  sumData <-rbind(e1sumData, e2sumData)
  
  #Average within participant
  ga_sumData <- sumData %>%
    mutate(meanCohs = (thisLow+thisHigh)/2) %>%
    group_by(ID, exp,block) %>%
    summarise(meanAcc = mean(accuracy), 
              meanCohs = mean(meanCohs))
  
  #Plot coherence over blocks, overlay with accuracy 
  avCoh_Acc_exp1 <- ggplot(ga_sumData[ga_sumData$exp == 1,], aes(x = block, y = meanAcc*100)) +
    geom_bar(stat="summary", position = position_dodge(), alpha=1, width = 0.75, fill = 'grey50', show.legend= F) +
    stat_summary(fun.data = mean_se, geom = "errorbar", alpha=1, width = 0.2, position = position_dodge(0.75))+ 
    theme_elis() + theme(axis.text.y.right = element_text(color = 'red'),axis.title.y.right = element_text(color = 'red'),
                         axis.line.y.right = element_line(color = "red"), axis.ticks.y.right = element_line(color = "red")) +   
    geom_line(aes(y = as.numeric(meanCohs)/10, color = "red"), show.legend = F) +
    geom_hline(yintercept = c(50)) +
    scale_x_continuous(name = 'Block',  breaks = c(1,5,10))+
    scale_y_continuous(name = 'Accuracy %', breaks = c(0,50,100), sec.axis = sec_axis(~., name = "Coherence %", breaks = c(0,50,100),))  +
    labs(title = 'Experiment 1')+
    facet_wrap(ID~.)
  
  avCoh_Acc_exp2 <- ggplot(ga_sumData[ga_sumData$exp == 2,], aes(x = block, y = meanAcc*100)) +
    geom_bar(stat="summary", position = position_dodge(), alpha=1, width = 0.75, fill = 'grey50', show.legend= F) +
    stat_summary(fun.data = mean_se, geom = "errorbar", alpha=1, width = 0.2, position = position_dodge(0.75), inherit.aes = TRUE)+ 
    geom_line(aes(y = as.numeric(meanCohs)/10, color = "red"), show.legend = F) +
    theme_elis() + theme(axis.text.y.right = element_text(color = 'red'),axis.title.y.right = element_text(color = 'red'),
                         axis.line.y.right = element_line(color = "red"), axis.ticks.y.right = element_line(color = "red")) +   
    geom_hline(yintercept = c(50)) +
    scale_x_continuous(name = 'Block', breaks = c(1,5,10))+
    scale_y_continuous(name = 'Accuracy %', breaks = c(0,50,100), sec.axis = sec_axis(~., name = "Coherence %", breaks = c(0,50,100),))  +
    labs(title = 'Experiment 2')+
    facet_wrap(ID~.)
  
  FigS1 <- cowplot::plot_grid(avCoh_Acc_exp1,avCoh_Acc_exp2, nrow = 1, ncol = 2)
  FigS1
  ggsave(paste(dir.figures, 'Figure S1.tiff'), FigS1,				
         width = 16, height = 9, units = "cm", dpi = 600)
  
} 

#Now exclude bad participants for the main analysis
e2data<- filter(e2sumdata, (ID != 1 & ID != 21))
e2sumData <- filter(e2sumData, (ID != 1 & ID != 21))
e2demo <-filter(e2demo, , (ID != 1 & ID != 21))

#Update data after participant exclusions
myData <- rbind(e1sumData, e2sumData)
sumData <-rbind(e1sumData, e2sumData)
demoData <- rbind(e1demo, e2demo)
myData <- rbind(e1sumData, e2sumData)
allData <- rbind(e1data, e2data)

#Starting coherence values- take first block only 
cohdata <- sumData%>%
  group_by(ID,exp, p1cohLev)%>%#
  filter(block == 1)%>%
  summarise(meanCoh = mean(p1coh)/10)

cohdata <- cohdata %>%
  group_by(exp, p1cohLev)%>%#
  summarise(sdCoh = sd(meanCoh),
            meanCoh = mean(meanCoh))

#Demographics 
demographics <- demoData %>%
  group_by(exp) %>%
  summarise(meanAge = mean(age), 
            stdAge = sd(age), 
            sum(gender == 'f'))


################################################################################################################################################################
# Eq. 1: Investigate effect of coherence on accuracy in 1pulse trials
################################################################################################################################################################
myData$gapdur_numeric <- as.numeric(as.character(myData$gapdur))

myData <- myData %>%
  group_by(exp, ID) %>%
  mutate(zcoh1 = scale(p1coh), 
         zcoh2 = scale(p2coh), 
         zgap = scale(gapdur_numeric))

# Exp. 1
 eq1_1 <- glmer(accuracy ~ 1 + zcoh1 + (1 + zcoh1|ID), 
                data = myData[myData$exp == 1 & myData$npulses == 1,], 
                family = binomial(link = 'logit'),
                control=glmerControl(optimizer="bobyqa"))
 
 car::Anova(eq1_1, type = c(2)) #For type 2 tests & p-vals
 summary(eq1_1)  #To get coefficients
 

 eq1_2 <- glmer(accuracy ~ 1 + zcoh1 + (1 + zcoh1|ID), 
                data = myData[myData$exp == 2 & myData$npulses == 1,], 
                family = binomial(link = 'logit'),
                control=glmerControl(optimizer="bobyqa"))
 
 car::Anova(eq1_2, type = c(2)) #For type 2 tests & p-vals
 summary(eq1_2)  #To get coefficients

################################################################################################################################################################
## Eq. 2: Investigate effect of coherence and gap duration on accuracy in 2pulse trials
################################################################################################################################################################  
#Exp 1 - after R1, used bobyqa optimizer and full random effects structure
eq2_1 <- glmer(accuracy ~ 1 + zcoh1 + zcoh2 + gapdur_numeric+ zcoh1:gapdur_numeric + zcoh2:gapdur_numeric +(1+ zcoh1 + zcoh2 + gapdur_numeric+ zcoh1:gapdur_numeric + zcoh2:gapdur_numeric|ID),
               data = myData[myData$exp == 1 & myData$npulses == 2,],
               family = binomial(link = 'logit'),
               control=glmerControl(optimizer="bobyqa",
                                    optCtrl=list(maxfun=2e5)),
               cov.prior = wishart)

car::Anova(eq2_1, type = c(2))
summary(eq2_1)

#Produces singularity
require(blme)
#Exp 1 - Try using weakly informative prior for covariance matrices to push away from singularity
#Chung et al. 2013/2015 recommend wishart prior
eq2_1 <- bglmer(accuracy ~ 1 + zcoh1 + zcoh2 + gapdur_numeric+ zcoh1:gapdur_numeric + zcoh2:gapdur_numeric +(1+ zcoh1 + zcoh2 + gapdur_numeric+ zcoh1:gapdur_numeric + zcoh2:gapdur_numeric |ID),
                data = myData[myData$exp == 1 & myData$npulses == 2,],
                family = binomial(link = 'logit'),
                control=glmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e5)),
                cov.prior = wishart)

car::Anova(eq2_1, type = c(2))
summary(eq2_1)

#Resolves singularity, but  convergence issues persist. Fit  with all available optimizers to assess robustness (as recommended in ?convergence):
allFit(show.meth.tab=TRUE) 
gm_all_eq2_1 <- allFit(eq2_1)
sum_gm_all_eq2_1 <- summary(gm_all_eq2_1)

plot.fixef.allFit(gm_all_eq2_1)

save(sum_gm_all_eq2_1, 'C:/Users/elisa/Desktop/GapsTask/')
write.csv(sum_gm_all_eq2_1$fixef, 'C:/Users/elisa/Desktop/GapsTask/EP1_eq2_1_new.csv')

#Exp 2 - after R1, used bobyqa optimizer and full random effects structure
eq2_2 <- glmer(accuracy ~ 1 + zcoh1 + zcoh2 + gapdur_numeric+ zcoh1:gapdur_numeric + zcoh2:gapdur_numeric +(1+ zcoh1 + zcoh2 + gapdur_numeric+ zcoh1:gapdur_numeric + zcoh2:gapdur_numeric |ID),
               data = myData[myData$exp == 2 & myData$npulses == 2,],
               family = binomial(link = 'logit'), 
               control=glmerControl(optimizer="bobyqa"))

car::Anova(eq2_2, type = c(2))
summary(eq2_2)

#Produces singularity
require(blme)
#Exp 2 - Try using weakly informative prior for covariance matrices to push away from singularity
#Chung et al. 2013/2015 recommend wishart prior 
eq2_2 <- bglmer(accuracy ~ 1 + zcoh1 + zcoh2 + gapdur_numeric+ zcoh1:gapdur_numeric + zcoh2:gapdur_numeric +(1+ zcoh1 + zcoh2 + gapdur_numeric+ zcoh1:gapdur_numeric + zcoh2:gapdur_numeric |ID),
               data = myData[myData$exp == 2 & myData$npulses == 2,],
               family = binomial(link = 'logit'),
               control=glmerControl(optimizer="bobyqa",
                                    optCtrl=list(maxfun=2e5)),
               cov.prior = wishart)

car::Anova(eq2_2, type = c(2))
summary(eq2_2)

#Resolves singularity, but  convergence issues persist. Fit  with all available optimizers to assess robustness (as recommended in ?convergence):
allFit(show.meth.tab=TRUE) 
gm_all_eq2_2 <- allFit(eq2_2)
sum_gm_all_eq2_2 <- summary(gm_all_eq2_2)

plot.fixef.allFit(gm_all_eq2_2)

save(sum_gm_all_eq2_2, 'C:/Users/elisa/Desktop/GapsTask/')
write.csv(sum_gm_all_eq2_2$fixef, 'C:/Users/elisa/Desktop/GapsTask/EP2_eq2_2.csv')

################################################################################################################################################################
## Eq. 3: Investigate effect of pulse number on accuracy 
################################################################################################################################################################  
myData$npulses <- as.factor(myData$npulses)

#Exp. 1
eq3_1 <- glmer(accuracy ~ 1 + npulses + zcoh1 + zcoh1:npulses + (1 + npulses + zcoh1+ zcoh1:npulses |ID),
               data = myData[myData$exp == 1 & myData$cohLevel != '1Med',], #Take only same-coherence trials
               family = binomial(link = 'logit'),
               control=glmerControl(optimizer="bobyqa"))

car::Anova(eq3_1, type = c(2))
summary(eq3_1)

eq3_1 <- bglmer(accuracy ~ 1 + npulses + zcoh1 + zcoh1:npulses + (1 + npulses + zcoh1 + zcoh1:npulses|ID),
                data = myData[myData$exp == 1 & myData$cohLevel != '1Med',], #Take only different coherence trials
                family = binomial(link = 'logit'),
                control=glmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e5)),
                cov.prior = wishart)

car::Anova(eq3_1, type = c(2))
summary(eq3_1)

#If convergence issues persist, attempt fitting with all available optimizers 
allFit(show.meth.tab=TRUE) 
gm_all_eq3_1 <- allFit(eq3_1)
sum_gm_all_eq3_1 <- summary(gm_all_eq3_1)

plot.fixef.allFit(gm_all_eq3_1)

write.csv(sum_gm_all_eq3_1$fixef, 'C:/Users/elisa/Desktop/GapsTask/EP2_eq3_1.csv')


#Follow-up tests 
eq3_11 <- glmer(accuracy ~ 1 + npulses  + (1 + npulses |ID),
               data = myData[myData$exp == 1 & myData$cohLevel == '0Low',], #Run on 0Low/2High separately
               family = binomial(link = 'logit'))

car::Anova(eq3_11, type = c(2))
summary(eq3_11)


#Exp. 2
eq3_2 <- glmer(accuracy ~ 1 + npulses + zcoh1 + zcoh1:npulses + (1 + npulses + zcoh1 + zcoh1:npulses|ID),
               data = myData[myData$exp == 2 & myData$cohLevel != '1Med',], #Take only different coherence trials
               family = binomial(link = 'logit'),
               control=glmerControl(optimizer="bobyqa",
                         optCtrl=list(maxfun=2e5)))
car::Anova(eq3_2, type = c(2))
summary(eq3_2)
#--> singular fit: need to regularize with BLME package (Chung 2013)
eq3_2 <- bglmer(accuracy ~ 1 + npulses + zcoh1 + zcoh1:npulses + (1 + npulses + zcoh1 + zcoh1:npulses|ID),
               data = myData[myData$exp == 2 & myData$cohLevel != '1Med',], #Take only different coherence trials
               family = binomial(link = 'logit'),
               control=glmerControl(optimizer="bobyqa",
                                    optCtrl=list(maxfun=2e5)),
               cov.prior = wishart)

car::Anova(eq3_2, type = c(2))
summary(eq3_2)

#If convergence issues persist, attempt fitting with all available optimizers 
allFit(show.meth.tab=TRUE) 
gm_all_eq3_2 <- allFit(eq3_2)
sum_gm_all_eq3_2 <- summary(gm_all_eq3_2)

plot.fixef.allFit(gm_all_eq3_2)

write.csv(sum_gm_all_eq3_2$fixef, 'C:/Users/elisa/Desktop/GapsTask/EP2_eq3_2.csv')

allFit(show.meth.tab=TRUE) 
gm_all <- allFit(eq3_2)
summary(gm_all)

################################################################################################################################################################
## Eq. 4: Investigate effect of pulse order on accuracy 
################################################################################################################################################################  
#Use normalised p1 & p2 coherence values to compute sums & differences 
myData <- myData %>%
  group_by(exp, ID) %>%
  mutate(c1plusc2 = p1coh+p2coh, 
         c2minc1 = p2coh-p1coh, 
         zc1plusc2 = scale(c1plusc2),
         zc2minc1 = scale(c2minc1))

#Exp. 1
eq4_1 <- glmer(accuracy ~ 1 + zc1plusc2 + zc2minc1 + (1+zc1plusc2 + zc2minc1|ID), 
               data = myData[myData$exp == 1  & myData$npulses == 2 & myData$cohLevel == '1Med',], # Take only two-pulse, mixed coherence trials. 
               family = binomial(link = 'logit'))
car::Anova(eq4_1, type = c(2))

summary(eq4_1)

#Exp. 2
eq4_2 <- glmer(accuracy ~ 1 + zc1plusc2 + zc2minc1 + (1+zc1plusc2 + zc2minc1|ID), 
               data = myData[myData$exp == 2  & myData$npulses == 2 & myData$cohLevel == '1Med',], # Take only two-pulse, mixed coherence trials. 
               family = binomial(link = 'logit'))

car::Anova(eq4_2, type = c(2))
summary(eq4_2)
