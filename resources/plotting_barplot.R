library(ggplot2)
library(viridis)
library(hrbrthemes)
library(dplyr)
library(tidyr)


#### Plotting weekly data

Hospital_data <- read.delim("DUKE_CMB_consolidated_monthly_results_updated_7_19_2021_cleaned.txt")

head(Hospital_data)
dim(Hospital_data)
summary(Hospital_data)

Hospital_week_data <- subset(Hospital_data, Week != "NA" )

head(Hospital_week_data)
dim(Hospital_week_data)
summary(Hospital_week_data)


#totals <- Hospital_week_data %>%
#  group_by(pango_tab, week) %>% mutate(value = n())  %>% ungroup() %>%   
#  mutate('relative'=unlist(by(data = value, INDICES = week,
#                              FUN = function(x) round(x/sum(x)*100, digits = 1))))


Duke_CMB_totals <- Hospital_week_data %>% 
  group_by(pango_tab, Week) %>% 
  count(pango_tab, Week)  %>% 
  group_by(Week) %>% 
  mutate(pct = n / sum (n), total= sum(n)) 

head(Duke_CMB_totals)

Duke_CMB_totals$pango_tab
Duke_CMB_totals$Week

Duke_CMB_totals$pango_tab <- factor(Duke_CMB_totals$pango_tab,levels = c("B.1.1.7" , "B.1.351",  "B.1.617.2", "P.1" , 
                                                                         "B.1.427", "B.1.429", "B.1.525","B.1.526" , "B.1.617.1" , "B.1.617.3" , "P.2", 
                                                                         "B.1.2",  "others"  ))
Duke_CMB_totals$Week <- factor(Duke_CMB_totals$Week,levels = c("week_1","week_2","week_3","week_4","week_5",
                                               "week_6","week_7","week_8","week_9","week_10","week_11",
                                               "week_12","week_13","week_14","week_15",
                                               "week_16", "week_17", "week_18", "week_19", 
                                               "week_20", "week_21", "week_22", "week_23",
                                               "week_24", "week_25",  "week_26", "week_27", 
                                               "week_28"))


summary(Duke_CMB_totals)


# Variants of Concern
levels(Duke_CMB_totals$pango_tab)[levels(Duke_CMB_totals$pango_tab)=="B.1.1.7"] <- "B.1.1.7 (Alpha)"
levels(Duke_CMB_totals$pango_tab)[levels(Duke_CMB_totals$pango_tab)=="B.1.351"] <- "B.1.351 (Beta)"
levels(Duke_CMB_totals$pango_tab)[levels(Duke_CMB_totals$pango_tab)=="B.1.617.2"] <- "B.1.617.2 (Delta)"
levels(Duke_CMB_totals$pango_tab)[levels(Duke_CMB_totals$pango_tab)=="P.1"] <- "P.1 (Gamma)"

# Variants of Interest
levels(Duke_CMB_totals$pango_tab)[levels(Duke_CMB_totals$pango_tab)=="B.1.427"] <- "B.1.427 (Epsilon)"
levels(Duke_CMB_totals$pango_tab)[levels(Duke_CMB_totals$pango_tab)=="B.1.429"] <- "B.1.429 (Epsilon)"
levels(Duke_CMB_totals$pango_tab)[levels(Duke_CMB_totals$pango_tab)=="B.1.525"] <- "B.1.525 (Eta)"
levels(Duke_CMB_totals$pango_tab)[levels(Duke_CMB_totals$pango_tab)=="B.1.526"] <- "B.1.526 (Iota)"
levels(Duke_CMB_totals$pango_tab)[levels(Duke_CMB_totals$pango_tab)=="B.1.617.1"] <- "B.1.617.1 (Kappa)"
levels(Duke_CMB_totals$pango_tab)[levels(Duke_CMB_totals$pango_tab)=="B.1.617.3"] <- "B.1.617.3 (--)"
levels(Duke_CMB_totals$pango_tab)[levels(Duke_CMB_totals$pango_tab)=="P.2"] <- "P.2 (Zeta)"



levels(Duke_CMB_totals$pango_tab)

levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_1"] <- "Week 1 (Jan 1, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_2"] <- "Week 2 (Jan 8, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_3"] <- "Week 3 (Jan 15, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_4"] <- "Week 4 (Jan 22, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_5"] <- "Week 5 (Jan 29, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_6"] <- "Week 6 (Feb 5, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_7"] <- "Week 7 (Feb 12, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_8"] <- "Week 8 (Feb 19, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_9"] <- "Week 9 (Feb 26, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_10"] <- "Week 10 (Mar 5, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_11"] <- "Week 11 (Mar 12, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_12"] <- "Week 12 (Mar 19, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_13"] <- "Week 13 (Mar 26, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_14"] <- "Week 14 (Apr 2, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_15"] <- "Week 15 (Apr 9, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_16"] <- "Week 16 (Apr 16, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_17"] <- "Week 17 (Apr 23, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_18"] <- "Week 18 (Apr 30, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_19"] <- "Week 19 (May 7, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_20"] <- "Week 20 (May 14, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_21"] <- "Week 21 (May 21, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_22"] <- "Week 22 (May 28, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_23"] <- "Week 23 (Jun 4, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_24"] <- "Week 24 (Jun 11, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_25"] <- "Week 25 (Jun 18, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_26"] <- "Week 26 (Jun 25, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_27"] <- "Week 27 (Jul 2, 2021)"
levels(Duke_CMB_totals$Week)[levels(Duke_CMB_totals$Week)=="week_28"] <- "Week 28 (Jul 9, 2021)"



levels(Duke_CMB_totals$Week)


library(RColorBrewer)
# Define the number of colors you want in case we need more colors
nb.cols <- 13
mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(nb.cols)

library(ggplot2)
ggplot(Duke_CMB_totals, aes(x= Week, fill = pango_tab, y = n)) + scale_fill_manual(values = mycolors,  drop = FALSE) +
  geom_bar(position="fill", stat="identity") +
  geom_text(aes(label=paste0(100*round(pct,digits=2),"%")), vjust=1.5, hjust=0.5,color = "black", size =2,  position=position_fill())+
  stat_summary(aes(label = total, y = 1.02), fun = mean, geom = "text", size = 3, vjust = 0) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +  theme(axis.text.x=element_text(angle=60,hjust=1)) +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank()) +
  ggtitle("Evolution of Coronavirus Variants (Duke-CMB)") +
  xlab("Week of Collection") + ylab( "Percentage of Variants per Week") + guides(fill=guide_legend(title="Lineage"))
  
  
########################
## plotting monthly data
  
head(Hospital_data)
dim(Hospital_data)
summary(Hospital_data)

Hospital_month_data <- subset(Hospital_data )

head(Hospital_month_data)
dim(Hospital_month_data)
summary(Hospital_month_data)


head(Hospital_month_data)

#totals <- Hospital_week_data %>%
#  group_by(pango_tab, week) %>% mutate(value = n())  %>% ungroup() %>%   
#  mutate('relative'=unlist(by(data = value, INDICES = week,
#                              FUN = function(x) round(x/sum(x)*100, digits = 1))))


totals_month <- Hospital_month_data %>% 
  group_by(pango_tab, Month) %>% 
  count(pango_tab, Month)  %>% 
  group_by(Month) %>% 
  mutate(pct = n / sum (n), total= sum(n)) 


levels(totals_month$Month)


totals_month$pango_tab <- factor(totals_month$pango_tab,levels = c("B.1.1.7" , "B.1.351",  "B.1.617.2", "P.1" , 
                                                                   "B.1.427", "B.1.429", "B.1.525","B.1.526" , "B.1.617.1" , "B.1.617.3" , "P.2", 
                                                                   "B.1.2",  "others"  ))

totals_month$Month <- factor(totals_month$Month,levels =c("march_2020", "april_2020", "may_2020", "june_2020", "july_2020", "august_2020", 
                                                          "september_2020", "october_2020", "november_2020", "december_2020",
                                                          "january_2021", "february_2021", "march_2021", "april_2021", "may_2021", "june_2021",
                                                          "july_2021") )
                             

totals_month
summary(totals_month)


# Variants of Concern
levels(totals_month$pango_tab)[levels(totals_month$pango_tab)=="B.1.1.7"] <- "B.1.1.7 (Alpha)"
levels(totals_month$pango_tab)[levels(totals_month$pango_tab)=="B.1.351"] <- "B.1.351 (Beta)"
levels(totals_month$pango_tab)[levels(totals_month$pango_tab)=="B.1.617.2"] <- "B.1.617.2 (Delta)"
levels(totals_month$pango_tab)[levels(totals_month$pango_tab)=="P.1"] <- "P.1 (Gamma)"

# Variants of Interest
levels(totals_month$pango_tab)[levels(totals_month$pango_tab)=="B.1.427"] <- "B.1.427 (Epsilon)"
levels(totals_month$pango_tab)[levels(totals_month$pango_tab)=="B.1.429"] <- "B.1.429 (Epsilon)"
levels(totals_month$pango_tab)[levels(totals_month$pango_tab)=="B.1.525"] <- "B.1.525 (Eta)"
levels(totals_month$pango_tab)[levels(totals_month$pango_tab)=="B.1.526"] <- "B.1.526 (Iota)"
levels(totals_month$pango_tab)[levels(totals_month$pango_tab)=="B.1.617.1"] <- "B.1.617.1 (Kappa)"
levels(totals_month$pango_tab)[levels(totals_month$pango_tab)=="B.1.617.3"] <- "B.1.617.3 (--)"
levels(totals_month$pango_tab)[levels(totals_month$pango_tab)=="P.2"] <- "P.2 (Zeta)"
levels(totals_month$pango_tab)


# Reoganize by month

levels(totals_month$Month)[levels(totals_month$Month)=="march_2020"] <- "March 2020"
levels(totals_month$Month)[levels(totals_month$Month)=="april_2020"] <- "April 2020"
levels(totals_month$Month)[levels(totals_month$Month)=="may_2020"] <- "May 2020"
levels(totals_month$Month)[levels(totals_month$Month)=="june_2020"] <- "June 2020"
levels(totals_month$Month)[levels(totals_month$Month)=="july_2020"] <- "July 2020"
levels(totals_month$Month)[levels(totals_month$Month)=="august_2020"] <- "August 2020"
levels(totals_month$Month)[levels(totals_month$Month)=="september_2020"] <- "September 2020"
levels(totals_month$Month)[levels(totals_month$Month)=="october_2020"] <- "October 2020"
levels(totals_month$Month)[levels(totals_month$Month)=="november_2020"] <- "November 2020"
levels(totals_month$Month)[levels(totals_month$Month)=="december_2020"] <- "December 2020"
levels(totals_month$Month)[levels(totals_month$Month)=="january_2021"] <- "January 2021"
levels(totals_month$Month)[levels(totals_month$Month)=="february_2021"] <- "February 2021"
levels(totals_month$Month)[levels(totals_month$Month)=="march_2021"] <- "March 2021"
levels(totals_month$Month)[levels(totals_month$Month)=="april_2021"] <- "April 2021"
levels(totals_month$Month)[levels(totals_month$Month)=="may_2021"] <- "May 2021"
levels(totals_month$Month)[levels(totals_month$Month)=="june_2021"] <- "June 2021"
levels(totals_month$Month)[levels(totals_month$Month)=="july_2021"] <- "July 2021"

levels(totals_month$Month)

library(RColorBrewer)
# Define the number of colors you want in case we need more colors
nb.cols <- 13
mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(nb.cols)


library(ggplot2)
ggplot(totals_month, aes(x= Month, fill = pango_tab, y = n)) + scale_fill_manual(values = mycolors,  drop = FALSE) +
  geom_bar(position="fill", stat="identity") +
  geom_text(aes(label=paste0(100*round(pct,digits=2),"%")), vjust=1.5, hjust=0.5,color = "black", size =2,  position=position_fill())+
  stat_summary(aes(label = total, y = 1.02), fun = mean, geom = "text", size = 3, vjust = 0) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +  theme(axis.text.x=element_text(angle=60,hjust=1)) +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank()) +
  ggtitle("Evolution of Coronavirus Variants (Duke-CMB)") +
  xlab("Month of Collection") + ylab( "Percentage of Variants per Month") + guides(fill=guide_legend(title="Lineage"))


  
