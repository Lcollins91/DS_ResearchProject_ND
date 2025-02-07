# Training_Data Information

The training data used in this project has been generated by MATLAB, saved as .csv files and .mat files and organized into folders here depending on the model for which it was created. These folders are:
  * _Line_Cut_: This model uses simulated line cut spectroscopy data perpendicular to a solid wall of CO molecules.
  * _Circle_Spec_: This model uses simulated spectroscopy data at the center of a solid-walled circle of CO molecules. For this system we do not have corresponding experimental data. 
  * _Hexagon_: This model uses simulated spectroscopy data at the center of a solid-walled hexagon of CO molecules built with the same proportions as the one we assembled within the STM. 
  * _Graphene_: This model uses simulated spectroscopy data at the "top" and "bond" site of a synthetic graphene lattice, which we compare to experimental data taken when Ken worked at Stanford. 


## Line_Cut Folder

* LineCutExpData180530.csv 
* LineCutTrainingData050318.csv/.mat
* LineCutTrainingData051418.csv/.mat
* LineCutTrainingData052118_fixed.mat
* LineCutTrainingData052118_fixedfit.csv


## Circle_Spec Folder

* CircleSpec_TrainingData052218.csv/.mat


## Hexagon Folder

Subheading do not represent subfolders, just subcategories of data within the Hexagon folder. 

#### CO Coordinates Files
* **hexagon.ngef**: This file is the Nanonis Grid Experiment File that was used to build the synthetic solid-walled hexagon.
* **NewHexagonNGEFPoints.csv**: This file is the conversion of the points in hexagon.ngef to a .csv file and was generated using the "Convert hexagon.ngef to csv.ipynb" notebook, since the .ngef file cannot be read directly by MATLAB.

#### Bias Files
* **HexagonBias.csv**: Bias ranging from -0.4V to 0.5V with 201 points total.
* **HexagonBias_v2.csv**: Bias ranging from -0.25V to 0.25V with 201 points total.
* **HexagonBias_v3.csv**: Bias ranging from -0.4V to 0.5V with 451 points total.
* **HexagonBias_v4.csv**: Bias ranging from -0.4V to 0.5V with 451 points total.

#### Experimental Data
* **HexagonExperimentalData053118_specPoints.csv**: 
* **HexagonExperimentalData053118_v2.csv**:
* **HexagonExperimentalData060618_v3.csv**:

#### Training Data
* **HexagonTrainingData052818_specPoints.csv/.mat**:
* **HexagonTrainingData053118_v2.csv/.mat**:
* **HexagonTrainingData060518_c3.csv/.mat**:
* **HexagonTrainingData060518_v4.csv/.mat**:
* **HexagonTrainingData060618_v5.csv**:


## Graphene Folder

Subheading do not represent subfolders, just subcategories of data within the Hexagon folder. 

#### Experimental Data
* ES_AG_Exp_data_SCALED.csv/.mat (the scale factor in this file has been hacked, multiplied by 0.964 to account for the faulty dispersion)
* ES_AG_Exp_data.csv/.mat
* Spec10a.mat

#### Training Data
* ES_AG_Spec_data_SCALED.csv/.mat (this uses the Scaled training points) 
* ES_AG_Spec_data.csv/.mat

