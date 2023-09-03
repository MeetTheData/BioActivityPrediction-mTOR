[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://bioactivityprediction-mtor.streamlit.app)
# BioActivityPrediction-mTOR (Drug Discovery)
Dysregulation of the mTOR pathway has been linked to the development and progression of cancer. This app allows you to predict the bioactivity of compounds towards inhibting the mTOR enzyme to find potential leads in Drug Discovery.

## Data Collection and Data Preprocessing 
The data represnts mtor single proteins in Homo Sapiens. The data was extracted from `ChemBL Database` using `chembl_webresource_client`. Handled missing values, duplicate values and standardized values with varied scales.

## Feature Extraction and Feature Selection
Performed feature extraction using `RDKit` library and PaDel-Descriptors software to extract Lipinksi Descriptors and PubChem fingerprints. Implement Feature Selection to drop low variance features using `Variance Threshold`.

## Model Building and Results
Having 179 features, `RandomForestRegressor` was modelled resulting in `Rsquare` of 0.68 and `Mean Square Error` of 0.55. 
