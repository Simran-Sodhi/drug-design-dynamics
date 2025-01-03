# drug-design-dynamics

## Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/Simran-Sodhi/drug-design-dynamics.git

## Project Structure
- `PLAS20K/`: Contains datasets related to the PLAS20K project  
- `Rshiny/`: Scripts and resources for the R Shiny application  
- `SPICE/`: Data and parsing scripts pertaining to the SPICE project  // we are not using this data anymore
- `machineLearningMethods/`: Machine learning models and training scripts  
- `metropolisMethod/`: Implementation of the Metropolis algorithm

## To download input data for Metropolis simulation from PLAS20K
#### Input: extended_PLAS20K.csv
- usePDBnames.go to extract all the pdb_id of protein-ligand complex used in PLAS20K (stored in PLAS20K_pdb_ids.txt)
- use batch_download.sh to grab all the pdb files from RSCB
- use splitPDB.go to seperate proteins and ligands (output two pdb files for proteins and ligands)
- use convert_pdb_to_mol2.sh (calls Open Babel) to convert all pdb files to mol2 files.

## Running the metropolis simulation from the go code
- You can use the metropolisMethod/main.go to run the metropolis simulation
- You need to provide data in metropolisMethod/Data. Some sample data is present there
- In main.go there are three options: one to simulate multiple ligands RunMultipleLigands(), one to get RMSD values: TestMethodRMSD() and the third for the R Shiny app: RShinyAppMain(args []string)
- All the outputs go into the metropolisMethod/Output folder


## R shiny

Interactive web app for predicting protein-ligand interations by evaluating their binding energies using Metropolis and machine learning simulations method.

## Requirements

- Set your Python path at Line 241 in app.R to ensure that the ML python script can be executed.
- Put all external data under folder either ./MCdata or ./MLdata for data to be accessible (like the example data).

## Usage

- Metropolis:<br />
Tha app takes one protein file (.pdb/.mol2) and multiple ligand files (.mol2 files uploaded together in a directory) as input,<br />after successfully uploaded the data, hit "Run Simulation" button, and you will get:<br />
  <br />
  (1) the plot of all protein-ligand pairs' binding energies<br />
  (2) the structure of protein-ligand pair with minimum binding energy shown<br />
  (need to preprocess the output protein and ligand file (.mol2) using Chimera or PyMOL to generate the video and save as ./output/results.mp4)<br />
   <br />
      
- Machine Learning:<br />
The app takes protein-ligand interaction dataset (.csv file) as input,<br />you can then select features you are interested in from ["electrostatic", "polar_solvation", "non_polar_solvation", "vdW"],<br />and select models you would like to use from ["RandomForest", "DecisionTree", "XGBoost", "LightGBM", "SVM"],<br />after successfully uploaded the data and selected the parameters, hit "Run Simulation" button, and you will get:<br />
  <br />
  (1) the evalution matrix of all models (evalution method includes: MSE, RMSE, MAE, R2, MAPE)<br />
  (2) feature importance tables generated by models you selected<br />

## Recorded Code demo available at: 

https://drive.google.com/drive/folders/1g_GTiWV2_l0lUbO9OTQ9euYVeeyWyaBW?usp=drive_link

