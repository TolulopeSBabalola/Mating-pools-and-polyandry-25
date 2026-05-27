

#Female dance flies return to the mating pool sooner than males: the role of female ornaments in mating rates and polyandry 

##StartHere---- 

###load packages 

library(dplyr)
library(tidyr)
library(ggplot2)
library(lmPerm)
library(readr)

###load data 
rawdata<-(read.csv("rawdata.csv"))
Lekreturnrates <- read_csv("Lekreturnrates.csv")
View(Lekreturnrates)

##filter out swarms between of <1

rr<-Lekreturnrates%>%filter(SwarmsBetween>1)
view(rr)

#isolate Female data
rrf<-rr%>%filter(Sex=="f")
summary(rrf)

#isolate Male data
rrm<-rr%>%filter(Sex=="m")
summary(rrm)

###GLM probability of recapture
prorecap<-glm(status~Sex, data=rawdata, family=binomial)
summary(prorecap)
simulationoutpu<-simulateResiduals(fittedModel = prorecap)
plot(simulationoutpu)

##figure 1 : OSR per day

summary_df <- rawdata %>%
  group_by(Date, Sex, SwarmNo) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = Sex,
    values_from = count,
    values_fill = 0
  ) %>%
  mutate(total = rowSums(across(where(is.numeric))))

# Define the desired order
order1 <- c("June 8 2020", "June 9 2020", "June 10 2020", "June 11 2020",
                   "June 12 2020","June 13 2020","June 14 2020","June 15 2020",
                   "June 16 2020","June 17 2020","June 18 2020")
# Reorder the factor levels in  data frame
summary_df$Date <- factor(summary_df$Date, levels = order1)
summary_df <- summary_df %>%
  arrange(Date) %>%        
  mutate(cum_total = cumsum(total))
summary_df <- summary_df %>%
  mutate(cum_prop = cum_total / max(cum_total))
summary_df <- summary_df %>%
  mutate(day_num = 1:n())

ggplot(data=summary_df, aes(x=SwarmNo))+
  geom_point(aes(y=m/total, size = total, colour = m/total))+
  geom_step(aes(y=cum_prop, group = 1))+
  theme_classic()+
  scale_colour_gradient(low = "orange", high = "maroon") +
  scale_size(range = c(5, 25))+
  scale_x_continuous(breaks = 1:12)+
  scale_y_continuous("Operational Sex Ratio", sec.axis = sec_axis(~ . * 387, name="Total Capture"))+
  theme_classic()+
  guides(
    colour = guide_colourbar(title = "", barwidth = 0.5, barheight = 8, position = "left"),
    size="none")+
  labs(y="Cumulative % Captured", x="Swarm")+
  theme( axis.title.x = element_text(size=16,face ="bold"),
         axis.text.x = element_text(size = 12, face = "bold"),
         axis.text.y = element_text(size = 12, face = "bold"),
         axis.title.y = element_text(size=16,face ="bold"),
         legend.text= element_blank())

###Differences in means and Significance testing-----
###means of M and F return 
mean_se(rrf$SwarmsBetween)
mean_se(rrm$SwarmsBetween)

se_rrm <- sd(rrm$SwarmsBetween)/sqrt(19)
se_rrf <- sd(rrf$SwarmsBetween)/sqrt(29)

isTRUE(mean(rrm$SwarmsBetween) > mean(rrf$SwarmsBetween)
       + (2 *se_rrm)) 

isTRUE(mean(rrf$SwarmsBetween) < mean(rrm$SwarmsBetween)
       - (2 *se_rrm)) 

##ordinary least squares 
diffmodel <- glm(SwarmsBetween ~ Sex, data = rr, family = poisson(link ="log"))
summary(diffmodel)

exp(coef(diffmodel))          # rate ratios
confint(diffmodel)            # CI on log scale
exp(confint(diffmodel))       # CI on RR scale
summary(diffmodel)$coefficients


##unpaired t-test 
t.test(rrm$DaysBetween,rrf$DaysBetween)
t.test(rrm$SwarmsBetween,rrf$SwarmsBetween)
t.test(rrm$SwarmsBetween, mu = mean(rrf$SwarmsBetween), alternative = "greater")

####Figure2-----
###Boxplot with over layed jitter 
ggplot(data=rr, aes(x=Sex, y=SwarmsBetween))+
  # make sure outlier points are not being repeated by boxplot
  geom_boxplot(aes(fill=Sex), outlier.shape = NA, alpha=0.2, show.legend = FALSE)+ 
  geom_jitter(aes(color=SwarmsBetween), size=10, alpha=0.7, width = 0.2)+
  #ggbeeswarm::geom_quasirandom(aes(color=SwarmsBetween), size=7, alpha=0.7)+
  # seeing if I can get the boxplots to have color
  scale_fill_manual(name = "", values = c("#930333", "#B05300")) + 
  scale_color_gradient(low="maroon", high="orange", space="Lab",guide="colorbar")+
  theme(panel.background = element_blank())+ 
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=12,face="bold"),legend.title = element_text(size=12),legend.text=element_text(size=10),
        axis.line=element_line("black"))+ 
  ylab("Interval Between Mating Pool Returns")+
  scale_x_discrete(breaks=c("f", "m"),labels=c("Female", "Male"))+labs(color = "Number of Swarms)")


###Permutation testing-----
##resample data and configure data frame of resampled means ( for figure 3)

#resampled F means x1000
fresample<-tibble(num=1:1000)%>%
  group_by(num)%>%
  mutate(means = mean(sample(rrf$SwarmsBetween, replace = TRUE)))

#resampled M means x1000
mresample<-tibble(num=1:1000)%>%
  group_by(num)%>%
  mutate(means=mean(sample(rrm$SwarmsBetween,replace=TRUE)))

##rename columns
fresample<-fresample%>%rename(F=means)
mresample<-mresample%>%rename(M=means)

#combining M and F resampled means into a data frame and converting it from wide to long
comboresample<-cbind(fresample$F,mresample$M)
comboresample<-as_tibble(comboresample)
comboresample<-comboresample%>%rename(Male=V2,Female=V1)
combo<-gather(comboresample,sex, means, Female,Male, factor_key = TRUE)

ggplot(data=rr,aes(fill=Sex))+geom_boxplot(aes(y=SwarmsBetween))

#####permutation test for swarms between return

#filter by sex for swarms between
allswarms<-rr%>%pull(SwarmsBetween)
fswarm<-rr%>%filter(Sex=="f")%>%pull(SwarmsBetween)
mswarm<-rr%>%filter(Sex=="m")%>%pull(SwarmsBetween)
length(fswarm)
length(mswarm)

##difference in means
observed<-mean(mswarm)-mean(fswarm)

#Permuation
N<-100000
result<-numeric(N)
for (i in 1:N) { index<-sample(length(allswarms),
                                size=length(fswarm),
                                replace = FALSE)
result[i]<-mean(allswarms[index])-mean(allswarms[-index])
}

##Plot 
ggplot()+
  geom_histogram(aes(result),color="blue", fill="lightblue")+
  geom_vline(xintercept=observed,color="red")+
  xlab("Male-Female")+
  ggtitle("Swarm Perm")

##pvalue
(sum(result >= observed) + 1)/(N + 1)

####Figure 3----
##1000x permutation of meals histogram, grouped by sex 
ggplot(data=combo,aes(x=means,fill=sex, color = sex))+
  geom_histogram(position = "identity", alpha = 0.6, bins = 100)+
  scale_color_manual(
    name = "Sex:",
    values = c("#5A032A","#B05300"),
    labels = c("Female", "Male"))+  
  scale_fill_manual(
    name = "Sex:",
    values = c("maroon","orange"),
    labels = c("Female", "Male"))+
  labs(x="Mean Return Rate after Resampling",y="Count")+
  theme(panel.background=element_rect(fill="white", color="white"))+ 
  # family is how you set font type
  theme(axis.title= element_text(color="black",family = "Arial", size=16,face="bold"),
        legend.title =  element_text(color="black",family = "Arial", size=16,face="bold"),
        axis.text = element_text(size=14),
        # adding in axis lines to make them more visible
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(limits = c(0, 75))







#BBnet model-----

#SetWD and load packages 
library(bbnet) #assuming it is installed

rm(list = ls(all = TRUE))

setwd("~/Library/CloudStorage/OneDrive-UniversityofToronto/Mating Pools/BabalolaMurray2025-data and code/BBModBabalolaMurray2025")

###Interaction Matrix-----
int.matrix<-read.csv("bbmodparameters.csv", header=T) 

####Scenarios-----

###Scenario1 is a lineage with nuptial gifts, of high value, females are highly dependent on these gifts and there is strong sensory bias for large females 

Scenario1<-read.csv('bbmodscenario1.csv', header=T) # read in files - should be an n x 2 matrix 

###Scenario2 is everything at 0, idk if this actiaually means anything

Scenario2<-read.csv('bbmodscenario2.csv', header=T) # can run up to 6 scenarios at once

##scenario 3 is a lineage with nuptial gifts, but gifts are not vauble, no strong sensory bias, mating rarely occurs without a gift.  

Scenario3<-read.csv('bbmodscenario3.csv', header=T)

###scenario 4 Just strong sensory bias
Scenario4<-read.csv('bbmodscenario4.csv', header=T)

###scenario5 Nuptial gifts are Not present and Mating is not dependent 
Scenario5<-read.csv('bbmodscenario5.csv', header=T)

###Predicted changes Figure S3-----
bbn.predict(bbn.model=int.matrix, priors1 = Scenario1, priors2 = Scenario2, priors3=Scenario3, priors4=Scenario4, priors5=Scenario5, figure=2, boot_max = 100, font.size = 12)

###Figure 1-----
##visualizing The whole network based on parameters 
##color code: 
#red (priors that are assumed preexisting and/or binary)
#grey (inputs that are variable and continuous, simplified for the model), 
#white(primarily outcomes, however one does effect the other)
#yellow (mating pool variables)

my_network <- read.csv('bbmoddiagram.csv', header=T) 

bbn.network.diagram(bbn.network = my_network, font.size = 2, arrow.size = 8, arrange = layout.circle)








