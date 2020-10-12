# adjusted-typhoid-incidence

This project includes the code used for the simulated portion of the manuscript "A Bayesian approach for estimating typhoid fever incidence from large-scale facility-based passive surveillance data" (https://www.medrxiv.org/content/10.1101/2020.10.05.20206938v1). In this study, we developed an approach using a Bayesian framework that adjusts the count of reported blood-culture-positive cases for healthcare seeking, blood culture collection, and blood culture detection while combining information from prior published studies. For more details about the methods, please refer to the full manuscript.

## Files included in this project

### Data
*NOTE:* Data provided in this GitHub project are NOT the site-specific data. These data have been simulated to evaluate the model formulation. 
`strataa_datsim.mat`-- This dataset contains 4 simulated scenarios for each of 5 variables. Each scenario has 100000 "individuals" simulated for each of the 5 variables, if applicable.
The 4 simulated scenarios are:
(1) low prob of seeking care, high prob of being tested, low prior antibiotic usage
(2) low prob of seeking care, high prob of being tested, high prior antibiotic usage
(3) high prob of seeking care, low prob of being tested, low prior antibiotic usage
(4) high prob of seeking care, low prob of being tested, high prior antibiotic usage
The 5 variables are : 
- `TFpositive`: Whether an individual had a positive blood culture result (if they were tested)
- `bctest`: Whether an individual received a blood culture test
- `fever`: Whether an individual had a fever
- `riskfactor`: Whether an individual had the typhoid fever risk factor
- `soughtcare`: Whether an individual sought care for the fever
Each of the 5 variables is an indicator where 1=yes and 0=no.

### R code
- `sim model full.R`: this code fits the simulated data using the FULL model. It formats the data, runs the model, then outputs some model diagnostics and results. You will need to specify how many "samples" from the healthcare utilization survey (HUS) you would like to use (lines 25-27). 
- `sim model simple approach.R`: tthis code fits the simulated data using the SIMPLE APPROACH model. It formats the data, runs the model, then outputs some model diagnostics and results. You will need to specify how many "samples" from the healthcare utilization survey (HUS) you would like to use (lines 23-25). 

## Running the code
To run this project, first specify how many "samples" from the healthcare utilization survey (HUS) you would like to use (lines 25-27 in full model, 23-25 in simpler model R code). Next, set the seed (if desired). Then, run the model all the way through. You may need to run the model several times (with different seeds) for convergence of all 4 scenarios. This is due to the normal-mixture model, which often poses difficulties with convergence.

## Authors

* **Maile T Phillips**
* **Joshua L Warren**
* **Virginia E Pitzer**
