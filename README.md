Welcome to the repository for Joel Clement's Honours Thesis Project, examining realised moment dynamics in Phase III and IV EUA futures contracts!
Please see the below list of files contained within this repository, and their respective purposes:

- DataProcessing.R: Processes raw data in 'Data/Raw/' and outputs processed data in 'Data/Processed/'. This includes calculating time until contract maturity, futures rolling, calculating log-returns and realised moments.
- MomentModelBuilding.R: Contains code creating functions for all aspects of analysis. This includes some common data transformations, model fitting, log-likelihood functions and model forecasting.
- EDA.Rmd: A high level analysis of processed data, particularly focused on descriptive statistics and trading volumes of futures contracts. Consists of mostly summary statistics and many high level plots.
- MomentModelCalibration.Rmd: Uses MomentModelBuilding.R functions to calibrate the TTM-HAR-RM and TTM-HAR-eps functions by finding the K that minimises AIC (also examining SIC and Adjusted R Squared).
- MomentModelAnalysis.Rmd: Uses calibrated models to analyse trends and perform the realised moments forecasting exercise.

Please see below for a description of data files referenced throughout the code (these are not uploaded nor available in this repository).

Data in 'Data/Raw/' contains raw .csv files after extraction from Refinitiv Datascope. The files are listed below:

- AllContracts.csv: A raw table containing a series of prices at 5-minute intervals for all traded EUA futures contracts.
- ContractDetails.csv: A raw table containing key contract details for each futures contract, including details of contract settlement and expiry dates.
- DecContracts.csv: A raw table containing the same data as AllContracts.csv, but just for December-expiring futures contracts.

Data in 'Data/Processed/' contains processed .csv files created by DataProcessing.R. The files are listed below:

- ExpiryMapping.csv: Uses data from ContractDetails.csv to create a mapping table of start and end contract dates for each RIC.
- AllReturns.csv: Table of log-returns and calculated number of days until maturity (using ExpiryMapping.csv and AllContracts.csv) for all futures contracts over all dates.
- AllRV.csv: Same as AllReturns.csv, with a realised variance series calculated.
- DecReturns.csv: Same as AllReturns.csv, but filtered to only include December-expiring contracts. Calcuated using the input DecContracts.csv.
- DecRV.csv: Modified version of DecReturns.csv, calculating realised moments for each contract.
- FirstDecNearby.csv: Rolled series (returns, RV, RS, RK) of nearest December contracts. Used for the bulk of subsequent analysis
- SecondDecNearby.csv: Rolled series (returns, RV, RS, RK) of second-nearest December contracts. Not directly used for any analysis, but was briefly investigated out of interest. Briefly used in EDA.Rmd.
