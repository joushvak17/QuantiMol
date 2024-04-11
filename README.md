___
Quantitative structure-activity relationship (QSAR) modeling software for molecular description calculations using RDKit and binary (active, inactive) activity prediction through machine learning classification method.

**Notes: The visualization/analysis code that is implemented in the InputFrame.py is not the most optimized. It works fine for now, but will/needs be updated later so it doesn't affect negative performance/usability and gives the user more options towards exploring the data.** 

Data Collection and Model Development:
- DataCollectionAndPreperation.ipynb: Data collection file for EGFR dataset from ChEMBL. Output file is the "EGFR_Data_Preprocessed.csv" 
- FeatureExtraction.ipynb: Feature calculation and extraction through RDKit. Output file is the "EGFR_Feature_Extraction.csv"
- ModelDevelopment.ipynb: XGBClassifier model development
- TestFileGeneration.ipynb: A different EGFR dataset from ChEMBL that can be used as a test dataset for the software. Output file is the "EGFR_Data_TestFile.csv"

Main Files:
- Interface.py: Main interface file that runs the Tkinter window application
- InputFrame.py: Tkinter frame that allows for data to be uploaded, deleted, computed, exported, analyzed and visualized
- TableFrame.py: Tkinter frame that shows a dataframe of calculated descriptor values and activity predictions
- OutputFrame.py: Tkinter frame that allows for user to navigate through the uploaded data, visualize 2D molecular images, and see individual values

Image Files (Main Files): All images are acquired from Flaticon.com
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

MLModels (Main Files):
- XGBClassifierEGFR.joblib: Machine learning classifier for binary activity class prediction (82% accuracy).
___
