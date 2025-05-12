<div align="center">
  <h1>QuantiMol 🧪</h1>
  <p align="center">
    <strong>
      <a href="https://github.com/joushvak17/QuantiMol/releases"><img src="https://img.shields.io/github/v/release/joushvak17/QuantiMol" alt="Release Version"></a>
      <a href="https://github.com/joushvak17/SeqCraft/issues"><img src="https://img.shields.io/github/issues/joushvak17/SeqCraft" alt="Issues"></a>
      <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue.svg" alt="License: MIT"></a>
      <br>
      Quantitative structure-activity relationship (QSAR) modeling software for molecular description calculations and activity prediction using machine learning.
    </strong>
  </p>
</div>

## 📑 Table of Contents
- [Usage](#️-usage)
- [Installation](#️-installation)
- [License](#-license)

## ⚗️ Usage
> [!NOTE]
> The visualization/analysis code that is implemented in the InputFrame.py is not the most optimized. It works fine for now, but will/needs be updated later so it doesn't affect negative performance/usability and gives the user more options towards exploring the data.

This software tool is an Quantitative structure–activity relationship modeling software that works with .csv data files. The user will want to define a .csv file that has one column and each row is a string format of SMILE data. The user will start the software where they can hit the "Upload Data" button to upload their .csv data file. The user is also given the ability to delete data files that they have uploaded by clicking on the file so it is highlighted and then hitting the "Delete Data" button.

Once the user uploads a file that they want, they can hit the "Compute Data" button. This process will take some time and there is a chance that the software will enter a Not Responding mode. Once all the data has been computed, a dataframe will show on the top right of the software with calculated descriptor values and the binary (active/inactive) prediction value at the last column of the dataframe. Below that dataframe, the user will see a 2D image of what the SMILE data looks like with the calculated values. The user can then navigate to the SMILE that they want to see using the "Next", "Previous", or the "Search" button which allows a user to enter the SMILE number. The software indexes the numbers from 1 till the length of the uploaded data, so users will need to find what the index is of a specific SMILE that they are looking for.

Finally a user can perform some data visualization using the "Data Visualize" button, where they will recieve a prompt in the lower left corner of the software about the type of visualizations they would like to perform. As of right now only boxplot, histogram, and pairwise scatter plots are supported. A user can also do a data analysis using the "Data Analysis" button, where they will recieve a prompt in the same location. The data analysis as of right now only supports a data summary.

## ⚙️ Installation
The current version of the software is v1.0.0 and can be downloaded through the [releases page](https://github.com/joushvak17/QuantiMol/releases), which also includes the source code. 

> [!IMPORTANT]
> The version that is included in the releases page is for the old software and a newer one is currently in development.

## 📜 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
