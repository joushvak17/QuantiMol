<div align="center">
  <h1>EGFR QSAR Modeling Software 1.0 🧬</h1>
  <p align="center"><strong>Quantitative structure-activity relationship (QSAR) modeling software for molecular description calculations</strong></p>
</div>

# ![alt text](https://github.com/joushvak17/QSAR-Modeling-Software/blob/master/assets/ui_image.jpg)

## 📑 File Outline

### [Data Collection and Model Development](https://github.com/joushvak17/EGFR-QSAR-Modeling-Software/tree/master/Data%20Collection%20and%20Model%20Development)

- DataCollectionAndPreperation.ipynb: Data collection file for EGFR dataset from ChEMBL. Output file is the "EGFR_Data_Preprocessed.csv"
- FeatureExtraction.ipynb: Feature calculation and extraction through RDKit. Output file is the "EGFR_Feature_Extraction.csv"
- ModelDevelopment.ipynb: XGBClassifier model development
- TestFileGeneration.ipynb: A different EGFR dataset from ChEMBL that can be used as a test dataset for the software. Output file is the "EGFR_Data_TestFile.csv"

### [Main Files](https://github.com/joushvak17/EGFR-QSAR-Modeling-Software/tree/master/Main%20Files)

- Interface.py: Main interface file that runs the Tkinter window application
- InputFrame.py: Tkinter frame that allows for data to be uploaded, deleted, computed, exported, analyzed and visualized
- TableFrame.py: Tkinter frame that shows a dataframe of calculated descriptor values and activity predictions
- OutputFrame.py: Tkinter frame that allows for user to navigate through the uploaded data, visualize 2D molecular images, and see individual values
- Interface.spec: File that defines the PyInstaller configurations
- hook -xgboost.py: File that defines additional operations, for XGBoost library, that need to be performed at runtime when using PyInstaller

### [Main Files/Image Files](https://github.com/joushvak17/EGFR-QSAR-Modeling-Software/tree/master/Main%20Files/Image%20Files)

- Icon.ico: Main interface file icon image. Made by Wanicon
- Analysis.png: Image file for data analysis button. Made by RaftelDesign
- Compute.png: Image file for compute button. Made by RaftelDesign
- CSV.png: Image file for export button. Made by Mpanicon
- Delete.png: Image file for delete button. Made by Futuer
- Next.png: Image file for next button. Made by Pixel Perfect
- Previous.png: Image file for previous button. Made by Pixel Perfect
- Search.png: Image file for search button. Made by Maxim Basinski Premium
- Upload.png: Image file for upload button. Made by Freepik
- Visualization.png: Image file for data visualization button. Made by Fajri Rama

### [Main Files/ML Models](https://github.com/joushvak17/EGFR-QSAR-Modeling-Software/tree/master/Main%20Files/MLModels)

- XGBClassifierEGFR.joblib: Machine learning classifier for binary activity class prediction (82% accuracy).

## ⚗️ Usage

**Notes: The visualization/analysis code that is implemented in the InputFrame.py is not the most optimized. It works fine for now, but will/needs be updated later so it doesn't affect negative performance/usability and gives the user more options towards exploring the data.**

This software tool is an Quantitative structure–activity relationship modeling software that works with .csv data files. The user will want to define a .csv file that has one column and each row is a string format of SMILE data. The user will start the software where they can hit the "Upload Data" button to upload their .csv data file. The user is also given the ability to delete data files that they have uploaded by clicking on the file so it is highlighted and then hitting the "Delete Data" button.

Once the user uploads a file that they want, they can hit the "Compute Data" button. This process will take some time and there is a chance that the software will enter a Not Responding mode. Once all the data has been computed, a dataframe will show on the top right of the software with calculated descriptor values and the binary (active/inactive) prediction value at the last column of the dataframe. Below that dataframe, the user will see a 2D image of what the SMILE data looks like with the calculated values. The user can then navigate to the SMILE that they want to see using the "Next", "Previous", or the "Search" button which allows a user to enter the SMILE number. The software indexes the numbers from 1 till the length of the uploaded data, so users will need to find what the index is of a specific SMILE that they are looking for.

Finally a user can perform some data visualization using the "Data Visualize" button, where they will recieve a prompt in the lower left corner of the software about the type of visualizations they would like to perform. As of right now only boxplot, histogram, and pairwise scatter plots are supported. A user can also do a data analysis using the "Data Analysis" button, where they will recieve a prompt in the same location. The data analysis as of right now only supports a data summary.

## 💻 Installation

The current version of the software is v1.0.0 and can be downloaded through the releases. You will be prompted to select a location to install the folder and you can create a desktop shortcut if you want. About ~500 MB of free disk space is required. Source code and executable can be found in releases.

## 📜 License

The EGFR QSAR Modeling Software is under the MIT License. See the [LICENSE](LICENSE) file for details.