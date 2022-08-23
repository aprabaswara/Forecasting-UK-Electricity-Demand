## Loading required library and defined required function
library(mgcv)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(data.table)
library(lubridate)
library(dplyr)

mape_compute <- function(ytrue,ypred){
  # Function to compute MAPE
  # Input: ytrue = actual values, ypred = predicted values
  # Output: MAPE
  
  error_ratio <- (ypred-ytrue)/ytrue
  mape <- mean(abs(error_ratio),na.rm = TRUE)
  return(mape)
}

rmse_compute <- function(ytrue,ypred){
  # Function to compute RMSE
  # Input: ytrue = actual values, ypred = predicted values
  # Output: RMSE
  
  mse <- mean((ypred-ytrue)^2,na.rm = TRUE)
  rmse <- sqrt(mse)
  return(rmse)
}

season_convert <- function(month_name){
  # Function to identify season in UK based on month
  # Input: month_name = month name
  # Output: name of season in UK based on month
  
  # Convert integer (1-12) to month name
  if(is.numeric(month_name)==T|is.integer(month_name)==T){
    month_name <- month.name[month_name]
  }
  
  # Identify winter season
  if(month_name %in% c(month.name[c(1,2,12)],month.abb[c(1,2,12)])){
    season <- "Winter"
  }
  
  # Identify spring season
  else if (month_name %in% c(month.name[c(3,4,5)],month.abb[c(3,4,5)])){
    season <- "Spring"
  }
  
  # Identify summer season
  else if (month_name %in% c(month.name[c(6,7,8)],month.abb[c(6,7,8)])){
    season <- "Summer"
  }
  
  # Identify autumn season
  else{
    season <- "Autumn"
  }
  return(season)
}

error_data <- function(ape,square_error,list_var){
  # Function to produces MAPE and RMSE summary table
  # Input:
  # ape = absolute percentage error, square_error = square error,
  # list_var = list of factor or categorical data
  # Output: MAPE and RMSE summary table
  error_df <- aggregate(cbind(mape=ape,mse=square_error), list_var, FUN=mean)
  error_df$rmse <- sqrt(error_df$mse)
  error_df$mse <- NULL
  return(error_df)
}

##---------------------------------------------------------------------------
## Loading the data and perform data preparation steps

# Load UK grid load data and convert date into julian day
load("UKL.RDA") 
day <- julian(UKL$date, origin = min(UKL$date))

# Add full name of day of week
UKL$day_desc <- wday(UKL$date, label=TRUE,abbr=FALSE)

# Add month name data
UKL$month_desc <- factor(month.name[UKL$month],levels=month.name)

# Convert day of week to integer (1 to 7)
UKL$day_num <- as.numeric(UKL$dow)

#categorize day based on weekdays and weekends
weekend_day <- c("Saturday","Sunday")
UKL$day_type <- ifelse(UKL$day_desc %in% weekend_day,"Weekends","Weekdays")
UKL$day_type <- factor(UKL$day_type)

#Add season data based on month
UKL$season <- unlist(lapply(UKL$month,season_convert))
UKL$season <- factor(UKL$season)

##--------------------------------------------------------------------------------
## Exploratory data analysis 

# Plot load against day
plot(day,UKL$load,type='l',xlab='Days',ylab='Load (MW)',main="Daily Load")

# Average load against time of day based on day type
p1 <- ggplot(UKL, aes(x=tod, y=load)) + aes(colour = day_type) + stat_summary(fun = mean, geom="line")+
  labs(x="Half Hour Period", y = "Average Load (MW)",colour="Day Type")+
  theme_classic()+
  ggtitle("Half Hourly Average Load (Based on Day Type)")+
  theme(plot.title = element_text(hjust = 0.5))

# Average load against time of day based on season
p2 <- ggplot(UKL, aes(x=tod, y=load)) + aes(colour = season) + 
  stat_summary(fun = mean, geom="line")+
  labs(x="Half Hour Period", y = "Average Load (MW)",colour="Season")+
  theme_classic()+
  ggtitle("Half Hourly Average Load (Based on Season)")+
  theme(plot.title = element_text(hjust = 0.5))

p1/p2

# Plot average daily temperature against day
temp_text <- expression("Average Daily Temperature ("*~degree*C*")")
plot(day,UKL$temp,type='l',xlab="Days",ylab=temp_text,
     main="Average Daily Temperature vs Days")


# Mean average daily temperature against time of day based on day type
tempavg_text <- expression("Mean Average Daily Temperature ("*~degree*C*")")
p3 <- ggplot(UKL, aes(x=tod, y=temp)) + aes(colour = day_type) + stat_summary(fun = mean, geom="line")+
  labs(x="Half Hour Period", y = tempavg_text,colour="Day Type")+
  theme_classic()+ggtitle("Weekdays and Weekends Mean Average Daily Temperature")+
  theme(plot.title = element_text(hjust = 0.5))

# Mean average daily temperature against time of day based on season
p4 <- ggplot(UKL, aes(x=tod, y=temp)) + aes(colour = season) + 
  stat_summary(fun = mean, geom="line")+
  labs(x="Half Hour Period", y = tempavg_text,colour="Season")+
  theme_classic()+ggtitle("Mean Average Daily Temperature (Based on Season)")+
  theme(plot.title = element_text(hjust = 0.5))

p3/p4

# Plot load against average daily temperature
plot(UKL$temp,UKL$load,xlab=temp_text,ylab="Load (MW)",
     main="Load Vs Average Daily Temperature")

##------------------------------------------------------------------------------
## GAM model for electricity forecasting

# Split the data into training and test data
train.index <- which(UKL$year != max(UKL$year))
train.dt <- data.table(UKL[train.index,])
test.dt <- data.table(UKL[-train.index,])

# Initialize list to store GAM model, RMSE, MAPE
model_list <- list()
mape_test <- list()
rmse_test <- list()
mape_train <- list()
rmse_train <- list()

# Initialize list to store training and test data
test_list <- list()
train_list <- list()

# Fit proposed GAM model at each half hour
for (i in unique(UKL$tod)){
  tod_train <- which(train.dt$tod==i)
  tod_test <- which(test.dt$tod==i)
  model_list[[i+1]] <- bam(load ~ s(timeCount, bs = "cr",k=10)+
                             s(load48,by=day_num,bs = "cr",k=10)+
                             s(toy, bs = "cc", k = 10)+
                             s(day_num, bs = "cc", k = 7)+
                             s(month, bs = "cc", k = 12)+
                             s(temp, bs = "cr",k=10)+
                             s(temp95, bs = "cr",k=10)+
                             s(year,bs="cr",k=3)+
                             ti(day_num,load48,bs=c("cr","cc"),k=c(7,10))+
                             ti(day_num,temp,bs=c("cr","cc"),k=c(7,10))+
                             ti(day_num,temp95,bs=c("cr","cc"),k=c(7,10)),
                           select=TRUE,data = train.dt[tod_train,],family = gaussian)
  
  # Compute test data RMSE and MAPE at each half hour
  rmse_test[[i+1]] <- rmse_compute(test.dt$load[tod_test],predict.gam(model_list[[i+1]],test.dt[tod_test,]))
  mape_test[[i+1]] <- mape_compute(test.dt$load[tod_test],predict.gam(model_list[[i+1]],test.dt[tod_test,]))
  
  # Compute training data RMSE and MAPE at each half hour
  rmse_train[[i+1]] <- rmse_compute(train.dt$load[tod_train],model_list[[i+1]]$fitted.values)
  mape_train[[i+1]] <- mape_compute(train.dt$load[tod_train],model_list[[i+1]]$fitted.values)
  
  # Store the training data at each half hour
  train_list[[i+1]] <- train.dt[tod_train,]
  
  # Add the training data predicted load and residual
  train_list[[i+1]]$resid <- model_list[[i+1]]$residuals
  train_list[[i+1]]$pred_load <- model_list[[i+1]]$fitted.values
  
  # Store the test data at each half hour
  test_list[[i+1]] <- test.dt[tod_test,]
  
  # Add the test data predicted load
  test_list[[i+1]]$pred_load <- predict.gam(model_list[[i+1]],test.dt[tod_test,])
}

# Concatenate the training data at each half hour
train_set <- bind_rows(train_list, .id = "column_label")

# Order the training data by date and half hour period
train_set <- train_set[order(date,tod),]

# Concatenate the test data at each half hour
test_set <- bind_rows(test_list, .id = "column_label")

# Compute residual, squared error, absolute percentage error
test_set$resid <- test_set$load-test_set$pred_load
test_set$square_error <- test_set$resid^2
test_set$ape <- 100*abs(test_set$resid/test_set$load)

# Order the test data by date and half hour period
test_set <- test_set[order(date,tod),]

# Model residual diagnostic plots
par(mfrow=c(2,2))
qqnorm(train_set$resid, ylab="Residuals",main = "Normal QQ-Plot")
qqline(train_set$resid,col="red")
plot(train_set$pred_load, train_set$resid, xlab = "Fitted values", 
     ylab = "Residuals", main = "Residuals vs Fitted values")
hist(train_set$resid, xlab = "Residuals", main = "Histogram of residuals")
plot(train_set$pred_load, train_set$load, xlab = "Fitted values", 
     ylab = "Response", main = "Response vs Fitted values")

##--------------------------------------------------------------------------
## Forecasting result table and plot

# Model fit comparison table
rmse_vec <- c(mean(unlist(rmse_train)),mean(unlist(rmse_test)))
mape_vec <- c(mean(unlist(mape_train)),mean(unlist(mape_test)))
error_df <- data.frame(mape = mape_vec, rmse = rmse_vec)
colnames(error_df) <- c("MAPE", "RMSE (MW)")
rownames(error_df) <- c("Training data", "Test data")
error_df

# Plot daily actual and predicted load in test data
test_index <- which(UKL$date %in% test_set$date)
plot(day[test_index],test_set$load,xlab="Days",ylab="Load (MW)",
     col="grey",type="l",main="Daily Actual Load and Predicted Load")
lines(day[test_index],test_set$pred_load,col="black")
legend("topright", legend=c("Actual Load", "Predicted Load"),
       col=c("grey", "black"),lty=1,bty='n',cex=0.7)

# Plot residual vs day in test data
plot(day[test_index],test_set$resid,xlab="Days",ylab="Residual (MW)",
     type="l",col="red",main="Daily Residuals")

# Plot residual vs half hour period
plot(test_set$tod,test_set$resid,xlab="Half Hour Period",
     ylab="Residual (MW)",type="l",col="grey",main="Half Hourly Residual")

# Plot MAPE vs half hour period based on day of week
list_var <- list(tod=test_set$tod,dow=test_set$day_desc)
dow_error <- error_data(ape,square_error,list_var)
ggplot(aes(x=tod,y=mape,group=dow,color=dow),data=dow_error)+
  geom_line()+theme_classic()+
  labs(x="Half Hour Period",y="MAPE (%)",colour="Day of Week")+
  ggtitle("Half Hourly MAPE (Based on Day of Week)")+
  theme(plot.title = element_text(hjust = 0.5))

# Plot RMSE vs half hour period based on day of week
ggplot(aes(x=tod,y=rmse,group=dow,color=dow),data=dow_error)+
  geom_line()+theme_classic()+
  labs(x="Half Hour Period",y="RMSE (MW)",colour="Day of Week")+
  ggtitle("Half Hourly RMSE (Based on Day of Week)")+
  theme(plot.title = element_text(hjust = 0.5))

# Plot MAPE vs half hour period based on day type
list_var <- list(tod=test_set$tod,day_type=test_set$day_type)
day_error <- error_data(ape,square_error,list_var)
p7 <- ggplot(aes(x=tod,y=mape,group=day_type,color=day_type),data=day_error)+
  geom_line()+theme_classic()+
  labs(x="Half Hour Period",y="MAPE (%)",colour="Day Type")+
  ggtitle("Half Hourly MAPE")+
  theme(plot.title = element_text(hjust = 0.5)) 

# Plot RMSE vs half hour period based on day of week
p8 <- ggplot(aes(x=tod,y=rmse,group=day_type,color=day_type),data=day_error)+
  geom_line()+theme_classic()+
  labs(x="Half Hour Period",y="RMSE (MW)",colour="Day Type")+
  ggtitle("Half Hourly RMSE")+
  theme(plot.title = element_text(hjust = 0.5)) 

ggarrange(p7,p8,nrow=1,ncol=2,common.legend=TRUE,legend="bottom")

# Plot MAPE vs half hour period based on season
list_var <- list(tod=test_set$tod,season=test_set$season)
season_error <- error_data(ape,square_error,list_var)
p9 <- ggplot(aes(x=tod,y=mape,group=season,color=season),data=season_error)+
  geom_line()+theme_classic()+
  labs(x="Half Hour Period",y="MAPE (%)",colour="Season")+
  ggtitle("Half Hourly MAPE")+
  theme(plot.title = element_text(hjust = 0.5)) 

# Plot RMSE vs half hour period based on season
p10 <- ggplot(aes(x=tod,y=rmse,group=season,color=season),data=season_error)+
  geom_line()+theme_classic()+
  labs(x="Half Hour Period",y="RMSE (MW)",colour="Season")+
  ggtitle("Half Hourly RMSE")+
  theme(plot.title = element_text(hjust = 0.5)) 

ggarrange(p9,p10,nrow=1,ncol=2,common.legend=TRUE,legend="bottom")