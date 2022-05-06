###################################################################################
########   Building LDA: Using Groundtruthed Carcasses                     ########
###################################################################################
#This code uses a subset of the real data used for analysis as an example of how to run code, and so results differ

library(MASS)
library(pROC)
library(tidyverse)
library(caret)

load("subset_carcassdata.Rdata") #subset of carcass data cleaned and covariates added, and Carcass Values assigned

#Carcass values here represent: 
# 1.	Very Large carcasses weighing over 2000kg
# 2.	Large carcasses between 900kg-1900kg
# 3.	Medium sized carcasses, including most ungulates (100kg-900kg)
# 4.	Small carcasses including small antelopes (below 100kg)
# 5.	Non-carcass


#########################
# Create Model
#########################

# Split the data into training (80%) and test set (20%)
theme_set(theme_classic())

set.seed(123)
training.samples <- mydata$Carcass.Value %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data <- mydata[training.samples, ]
test.data <- mydata[-training.samples, ]

#########################
model <- lda(Carcass.Value~Total_Hours + Total_Birds + Dist.to.Rivers + Dist.to.Rangers + Tree.cover, data = train.data)
model

# Check model
predictions <- model %>% predict(test.data)
mean(predictions$class==test.data$Carcass.Value)  # Model accuracy


# Model ROC/AUC
roc_obj<-roc(test.data$Carcass.Value, rowSums(predictions$posterior))
auc(roc_obj)
plot(roc_obj, print.auc=TRUE)

#Class ROC/AUC distinctions 
roc_obj2<-roc(test.data$Carcass.Value==1, predictions$posterior[,1],  quiet = FALSE)
auc(roc_obj2)
plot(roc_obj2, print.auc=TRUE)

roc_obj3<-roc(test.data$Carcass.Value==2, predictions$posterior[,2],  quiet = FALSE)
auc(roc_obj3)
plot(roc_obj3, print.auc=TRUE)

roc_obj4<-roc(test.data$Carcass.Value==3, predictions$posterior[,3],  quiet = FALSE)
auc(roc_obj4)
plot(roc_obj4, print.auc=TRUE)

roc_obj4<-roc(test.data$Carcass.Value==4, predictions$posterior[,3],  quiet = FALSE)
auc(roc_obj4)
plot(roc_obj4, print.auc=TRUE)

roc_obj4<-roc(test.data$Carcass.Value==6, predictions$posterior[,3],  quiet = FALSE)
auc(roc_obj4)
plot(roc_obj4, print.auc=TRUE)


###########################################################################
###### Using LDA on unchecked cluster points to predict carcass type ######
###########################################################################
# We used our validated LDA to predict carcass class type for each of the clusters that had not been verified by ground-truthing

load("subset_uncheckeddata.Rdata") #subset of unchecked cluster points with covariates added

# Make carcass class predictions

predictions2 <- model %>% predict(unchecked_points)  #for tool output, same format as test data 
table(predictions2$class)

#Add predictions to dataframe
pred <- as.data.frame(predictions2$posterior) 
unchecked_points <- merge(unchecked_points, predictions2, by = "row.names")
unchecked_points <- unchecked_points[-c(1,13:16)]

################################################
### Pick class threshold cut-off/inclusion   ###
################################################

#How many points per cut-off?
sum(unchecked_points$posterior.5>0.2)
sum(unchecked_points$posterior.5>0.5)
sum(unchecked_points$posterior.5>0.6)
sum(unchecked_points$posterior.5>=.4)

#Delete at chosen cut-off
unchecked_points <- unchecked_points[!(unchecked_points$posterior.5 >= .4),]


##############################################
###  Add back in Groundtruthed carcasses   ###
##############################################

# Match Columns
unchecked_points <- unchecked_points[-c(7:11)]
colnames(unchecked_points)[6] <- "Carcass.Value"

#Merge
Carcass.dataset <- rbind(mydata, unchecked_points)



### Add additional variables for Capture-Recapture analyses: Lat/Long data, spatial variables, region and protected area
###  This code not provided as it contains sensitive data



