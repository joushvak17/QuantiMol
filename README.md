___
Quantitative structure-activity relationship (QSAR) modeling software for molecular description calculations and activity prediction through machine learning methods.

**Notes: The visualization/analysis code that is implemented in the InputFrame.py is not the most optimized. It works fine for now, but will/needs be updated later so it doesn't affect negative performance/usability and gives the user more options towards exploring the data.** 

Data Collection and Model Development:
- DataCollectionAndPreperation.ipynb: Data collection file for EGFR dataset from CHEMBL
- FeatureExtraction.ipynb: Feature calculation and extraction through RDKit
- ModelDevelopment.ipynb: XGBClassifier model development

Main Files:
- Interface.py: Main interface file that runs the Tkinter window application
- InputFrame.py: Tkinter frame that allows for data to be uploaded, deleted, computed, exported, analyzed and visualized
- TableFrame.py: Tkinter frame that shows a dataframe of calculated descriptor values and activity predictions
- OutputFrame.py: Tkinter frame that allows for user to navigate through the uploaded data, visualize 2D molecular images, and see individual values

Image Files:
- Icon.ico: Main interface file icon image
- Analysis.png: Image file for data analysis button
- Compute.png: Image file for compute button
- CSV.png: Image file for export button
- Delete.png: Image file for delete button
- Next.png: Image file for next button
- Previous.png: Image file for previous button
- Search.png: Image file for search button
- Upload.png: Image file for upload button
- Visualization.png: Image file for data visualization button

MLModels:
- XGBClassifierEGFR.joblib: Main machine learning classifier for binary activity class
___
