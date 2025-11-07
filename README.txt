This repository contains the R scripts to reproduce the numerical experiments from the paper "Learning extremal graphical structures in high dimensions" by S. Engelke, M. Lalancette and S. Volgushev (Ann. Statist., 2025+).

In the file simulations/simulation_paper.R, line 22, change the value of sim_setting to the name of any of the JSON files found in simulations/config/ and source the file. This will produce a .RDS file in simulations/output/, the name of which contains information on the simulation configurations. Once all these output files are produced, copy them to figures/data and source figures/plots_paper.R to produce all the simulation figures in figures/.

For the data applications, first source the file applications/generate_currency_data.R to create the pre-processed data for the financial data application. Then, source the files applications/application_danube26.R and applications/application_exchange.R to produce all the figures in applications/figures/.
