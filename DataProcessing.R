##LIBRARY INVOKATION##
library(tidyverse) # ggplot, dplyr, and other functionality
library(lubridate) # date/timestamp data handling
library(qlcal) # trading calendar tools
library(zoo) # rollmean()
start <- Sys.time()

##SETTING PARAMETERS##
setwd('..')

maxTTM <- 252
#Creating ICE London Business Calendar
setCalendar("UnitedKingdom")

##DEFINING NA FUNCTIONS##
na.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
na.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)
na.first <- function(x) ifelse( !all(is.na(x)), first(x, na_rm=T), NA)
na.last <- function(x) ifelse( !all(is.na(x)), last(x, na_rm=T), NA)
na.sum <- function(x) ifelse( !all(is.na(x)), sum(x, na.rm=T), NA)
na.mean <- function(x) ifelse( !all(is.na(x)), mean(x, na.rm=T), NA)


##Reading and Processing Raw ContractDetails.csv##
ContractDetails <- read.csv('Data/Raw/ContractDetails.csv') %>%
  mutate(Exchange.Code = as.factor(Exchange.Code), Trading.Status = as.factor(Trading.Status)) %>%
  select(-c(IdentifierType, Identifier, Company.Name, Exercise.Style, Domicile, Strike.Price, Ticker)) %>%
  arrange(RIC)

ExpiryMapping <- ContractDetails %>% filter(
  (nchar(RIC) == 6 & substr(Last.Trading.Day, 3, 3) == dplyr::if_else(substr(RIC, 6, 6) > 2, 2, 3)) # ensure that a 6 digit ric actually ends in the year range of 2023-2029
  | substr(RIC, 8, 8) == substr(as.character(Last.Trading.Day), 3, 3), # ensure that an 8 digit ric ends in the year as shown in last.trading.day
  First.Trading.Date != '', Last.Trading.Day != '', as.integer(substr(First.Trading.Date, 1, 4)) > 2005) %>%
  select(RIC, Security.Description, First.Trading.Date, Last.Trading.Day) %>%
  group_by(RIC) %>%
  summarise(First.Trading.Date = first(First.Trading.Date), Last.Trading.Day = last(Last.Trading.Day))

##Reading and Processing Raw AllContracts.csv##
AllContracts <- read.csv('Data/Raw/AllContracts.csv') %>%
  mutate(Date.Time = as.POSIXct(Date.Time, format = '%Y-%m-%dT%H:%M:%OSZ', tz = 'UTC'), RIC = X.RIC) %>% #renaming some columns + formatting DT
  mutate(TimeInterval = floor_date(Date.Time, unit = "15 minutes")) %>% #form 15 minute intervals
  group_by(TimeInterval, RIC) %>% #grouping by RIC and 15 min intervals
  summarise(Date.Time = as.POSIXct(min(Date.Time)), Open = na.first(Open), Last = na.last(Last), High = na.max(High), Low = na.min(Low), Num_Trades = na.sum(No..Trades), .groups = "keep") %>% #aggregating the numeric columns. removing NA values to ensure correct calcs
  dplyr::select(-Date.Time) %>%
  mutate(ExpiryYear = as.integer(substr(RIC, 6, 6)) + if_else(TimeInterval > paste0('201', substr(RIC,6,6), '-12-31'), 2020, 2010)) %>% # add row "ExpiryYear" = if(year(date.time) > 201x, then 202x, else 201x). don't need to concern with 2020 vs 2030 as 2030 contract not yet traded
  mutate(RIC = paste0(RIC, if_else(ExpiryYear < 2023, paste0('^', substr(as.character(ExpiryYear), 3, 3)), ''))) %>% # modify RICs to include decade number
  group_by(RIC) %>%
  arrange(RIC, TimeInterval)
AllContracts$Num_Trades[is.na(AllContracts$Num_Trades)] = 0

AllReturns <- AllContracts %>% 
  inner_join(ExpiryMapping, by = 'RIC') %>%
  mutate(RIC = as.factor(RIC), DaystoMaturity = businessDaysBetween(as.Date(TimeInterval), as.Date(Last.Trading.Day)), DaysTraded = businessDaysBetween(as.Date(First.Trading.Date), as.Date(TimeInterval)), absRet = Last - lag(Last), logRet = log(Last) - log(lag(Last))) %>%
  ungroup() %>%
  select(RIC, Time = TimeInterval, Last, absRet, logRet, Num_Trades, DaystoMaturity, DaysTraded)

AllRV <- AllReturns %>%
  mutate(Date = as.Date(Time)) %>%
  group_by(Date, RIC) %>%
  summarise(RV = na.sum(logRet^2), Last = na.last(Last), absRet = na.sum(absRet), logRet = na.sum(logRet), Volume = na.sum(Num_Trades), DaystoMaturity = first(DaystoMaturity), DaysTraded = first(DaysTraded), NumBlocksTraded = na.sum(ifelse(Num_Trades > 0, 1, 0)))

##Reading and Processing Raw DecContracts.csv##
DecContracts <- read.csv('Data/Raw/DecContracts.csv') %>%
  mutate(Date.Time = as.POSIXct(Date.Time, format = '%Y-%m-%dT%H:%M:%OSZ', tz = 'UTC'), RIC = X.RIC) %>% #renaming some columns + formatting table
  mutate(TimeInterval = floor_date(Date.Time, unit = "15 minutes")) %>% #form 15 minute intervals
  group_by(TimeInterval, RIC) %>% #grouping by RIC and 15 min intervals
  summarise(Date.Time = as.POSIXct(min(Date.Time)), Open = na.first(Open), Last = na.last(Last), High = na.max(High), Low = na.min(Low), Num_Trades = na.sum(No..Trades), .groups = "keep") %>% #aggregating the numeric columns. removing NA values to ensure correct calcs
  dplyr::select(-Date.Time) %>%
  mutate(ExpiryYear = as.integer(substr(RIC, 6, 6)) + if_else(TimeInterval > paste0('201', substr(RIC,6,6), '-12-31'), 2020, 2010)) %>% # add row "ExpiryYear" = if(year(date.time) > 201x, then 202x, else 201x). don't need to concern with 2020 vs 2030 as 2030 contract not yet traded
  mutate(RIC = paste0(RIC, if_else(ExpiryYear < 2023, paste0('^', substr(as.character(ExpiryYear), 3, 3)), ''))) %>% # modify RICs to include decade number
  group_by(RIC) %>%
  arrange(RIC, TimeInterval)
DecContracts$Num_Trades[is.na(DecContracts$Num_Trades)] = 0

DecReturns <- DecContracts %>% 
  inner_join(ExpiryMapping, by = 'RIC') %>%
  mutate(RIC = as.factor(RIC), DaystoMaturity = businessDaysBetween(as.Date(TimeInterval), as.Date(Last.Trading.Day)), DaysTraded = businessDaysBetween(as.Date(First.Trading.Date), as.Date(TimeInterval)), absRet = Last - lag(Last), logRet = log(Last) - log(lag(Last))) %>%
  ungroup() %>%
  select(RIC, Time = TimeInterval, Last, absRet, logRet, Num_Trades, DaystoMaturity, DaysTraded)

DecRV <- DecReturns %>%
  mutate(Date = as.Date(Time)) %>%
  group_by(Date, RIC) %>%
  summarise(RV = na.sum(logRet^2), m3 = na.sum(logRet^3), m4 = na.sum(logRet^4), RS = sqrt(sum(!is.na(logRet))) * m3/RV^(3/2), RK = sum(!is.na(logRet)) * m4/RV^2, Last = na.last(Last), absRet = na.sum(absRet), logRet = na.sum(logRet), Volume = na.sum(Num_Trades), DaystoMaturity = first(DaystoMaturity), DaysTraded = first(DaysTraded), NumBlocksTraded = sum(ifelse(Num_Trades > 0, 1, 0))) %>%
  ungroup()

TradingDates <- data.frame(cbind(ExpiryMapping, Last.Recorded.Day = as.Date(unlist(lapply(as.Date(ExpiryMapping$Last.Trading.Day), advanceDate, days = -10))))) %>%
  filter(substr(RIC, 5, 5) == 'Z') %>% #filter just December contracts
  arrange(Last.Recorded.Day) %>%
  mutate(First.Recorded.Day = lag(Last.Recorded.Day), Last.Recorded.Day.2nd = lag(Last.Recorded.Day), First.Recorded.Day.2nd = lag(First.Recorded.Day))

FirstDecNearby <- DecRV %>%
  inner_join(TradingDates, by = 'RIC') %>%
  filter(Date > First.Recorded.Day & Date <= Last.Recorded.Day)

SecondDecNearby <- DecRV %>%
  inner_join(TradingDates, by = 'RIC') %>%
  filter(Date > First.Recorded.Day.2nd & Date <= Last.Recorded.Day.2nd)

##WRITING OUTPUTS## 
write_csv(ExpiryMapping, 'Data/Processed/ExpiryMapping.csv')
write_csv(AllReturns, 'Data/Processed/AllReturns.csv')
write_csv(AllRV, 'Data/Processed/AllRV.csv')
write_csv(DecReturns, 'Data/Processed/DecReturns.csv')
write_csv(DecRV, 'Data/Processed/DecRV.csv')
write_csv(FirstDecNearby, 'Data/Processed/FirstDecNearby.csv')
write_csv(SecondDecNearby, 'Data/Processed/SecondDecNearby.csv')
print(Sys.time() - start)


