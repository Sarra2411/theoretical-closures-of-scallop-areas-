# Theoretical closures of scallop fishing areas
## Description
The aim of this R code is to simulate the closures (number and duration) of scallops fishing areas due to shellfish contaminations by harmful algal blooms. This R script uses the _in situ_ dataset from the French monitoring network
REPHYTOX (French Monitoring Network for Phycotoxins in marine organisms) and information related to the monitoring and the management of scallops fishery in the Eastern English Channel to generate ***theoretical closures*** of the production areas. The processing consists of two main steps: 
1. The creation of the _closures calendar_ using a binary coding based on the comparison of the phycotoxin values in scallop with the regulatory thresholds.
2. The simulation of the number and the duration of bans according to regulations. 
## Input data
Input data used in this code are associated to this code. It concerns data on phycotoxins concentrations in scallops (REPHYTOX data), data on regulatory thresholds for each toxin familly, and data about the fishing season calendar of scallops in the French coasts of the Eastern Channel.  
## How to use the code
To run this R code you need to download the input data in your desktop and import it into your R space. The package _'dplyr'_ is required to run the code, if you don't have it you need to install it using ``install.packages("dplyr")``.  
## Output data
The output data of this code consists of two tables containing the number and duration of closures for each production area per week (1st table) and per fishing season (2nd table). 
# Licence
This project is licensed under the terms of the GNU General Public License v3.0 license. 
