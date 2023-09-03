import streamlit as st
import altair as alt
import pandas as pd
from PIL import Image
import subprocess
import os
import base64
import pickle
import bz2file as bz2
from sklearn.feature_selection import VarianceThreshold

with open('app2.css') as f:
    st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)

# Molecular descriptor calculator
def desc_calc():
    #Performs the descriptor calculation
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.remove('molecule.smi')

# File download
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href

def decompress_pickle(file):
    data = bz2.BZ2File(file, 'rb')
    data = pickle.load(data)
    return data

# Model building
def build_model(input_data):
    # Reads in saved regression model
    #load_model = pickle.load(open('model.pkl', 'rb'))
    load_model = decompress_pickle('model.pbz2')
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    st.header('**Prediction output**')
    st.markdown(""" 
    We use a RandomForestRegressor to predict the bioactivity (pIC50). The model explains `69.29%` variabilty in target data. the model has a Mean Square Error of `0.55`, indicating an average difference of 0.55 from actual values exists.
    """)
    prediction_output = pd.Series(prediction, name='Predicted')
    molecule_name = pd.Series(desc['Name'], name='Chembl_ID')
    df = pd.concat([molecule_name, prediction_output], axis=1)
    st.write(df)
    alt.themes.register('black_marks', black_marks)
    alt.themes.enable('black_marks')
    st.header('**Actual vs Predicted**')
    st.markdown(""" 
    The graph shows difference in the actual and predicted bioactivity values for our sample compounds.
    """)
    graph_dat = df.merge(st.session_state.df2, on='Chembl_ID', how='inner')
    graph_dat = graph_dat.melt('Chembl_ID', var_name='Legend', value_name='pIC50')
    chart = alt.Chart(graph_dat).mark_line().encode(
        x=alt.X('Chembl_ID:N'),
        y=alt.Y('pIC50:Q'),
        color=alt.Color('Legend:N'))
    st.altair_chart(chart, use_container_width=True)

    
def black_marks():
    return {
        "config": {
            "axis": {
              "labelFontSize": 14,
              "labelFontWeight": 400,
              "labelColor": "#e6eaf1",
              "titleFontWeight": 600,
              "titleFontSize": 18,
              "gridColor": "#4DD5F3",
              "titleColor": "#e6eaf1",
              "labelPadding": 16,
            },
            "range": {
              "category": [
                "#ffffff",
                "#414FEA"],
            },
        },
     "padding": {"bottom": 20, "top": 40, "right": 20, "left": 20}
    }

def add_bg_from_url():
    st.markdown(
         f"""
         <style>
         .stApp {{
             background-image: url("https://scitechdaily.com/images/DNA-Change-Concept.gif");
             background-attachment: fixed;
             background-size: cover
         }}
         [data-testid='stHeader'] {{
             background-color: rgba(0,0,0,0);
         }}
         </style>
         """,
         unsafe_allow_html=True
     )

add_bg_from_url() 

def refresh():
    if  st.session_state.count_s > 0:
        del st.session_state.sample_clicked
        del st.session_state.df

def callback():
    st.session_state.sample_clicked = True
    


# Logo image
#image = Image.open('img.jpeg')

#st.image(image, use_column_width=True)

# Page title
st.markdown("""
# Bioactivity Prediction App

Dysregulation of the `mTOR` pathway has been linked to the development and progression of cancer. This app allows you to predict the bioactivity of compounds towards inhibting the `mTOR` enzyme to find potential leads in Drug Discovery.
""")

if 'df' not in st.session_state and 'sample_clicked' not in st.session_state:
    sample_data = pd.read_csv('sample_data.csv')
    data = sample_data[['canonical_smiles_new','molecule_chembl_id','pIC50']].sample(5)
    st.session_state.df = data[['canonical_smiles_new','molecule_chembl_id']]
    st.session_state.df2 = data[['molecule_chembl_id','pIC50']]
    st.session_state.df3 = sample_data[['molecule_chembl_id', 'MW', 'LogP', 'NumHDonors', 'NumHAcceptors']]
    st.session_state.df2.rename(columns={'molecule_chembl_id': 'Chembl_ID', 'pIC50':'Actual'}, inplace=True)
    st.session_state.sample_clicked = False
    st.session_state.count_s = 0
    
st.sidebar.header('Sample Data')
st.sidebar.caption('Click button to generate SMILES input data.')
sample = st.sidebar.button('Sample', on_click=refresh)

if sample or st.session_state.sample_clicked:
    st.session_state.count_s += 1
    st.header('**Generated Sample data**')
    st.markdown(""" 
    Canonical smiles is a standard notation for representing the molecular structure of a chemical compound using a short and human-readable string of characters, making it easier to compare and analyze chemical structures.
    """)
    st.write(st.session_state.df)
    st.sidebar.header('Predict Data')
    st.sidebar.caption('Click button to predict bioactivity.')
    predict = st.sidebar.button('Predict', on_click=callback)

    if predict:
        st.session_state.sample_clicked = False
        st.session_state.df.to_csv('molecule.smi', sep = '\t', header = False, index = False)
        st.header('**Lipinski Descriptors**')
        st.markdown(""" 
           We calculate 4 Lipinski Descriptors using `RDKit` library. These descriptors help us evaluate molecules with certain sufficient properties making them potential features to predict Bioactivity. The features are Molecular Weight (MW), Partition Coefficient (LogP), Number of Hydrogen Donors (NumHDonors), Number of Hydrogen Acceptors (NumHAcceptors).
    """)
        lip = st.session_state.df.merge(st.session_state.df3, on='molecule_chembl_id', how='inner')
        st.write(lip.iloc[:,1:])
        #st.header('**Original input data**')
        #st.write(st.session_state.df)

        with st.spinner("Calculating descriptors..."):
            desc_calc()

        # Read in calculated descriptors and display the dataframe
        st.header('**Calculated PaDel descriptors**')
        st.markdown(""" 
           We use [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) [[Read the Paper]](https://doi.org/10.1002/jcc.21707) to calculate molecule fingerprints. PaDel descriptors are numerical values that capture various characteristics of chemical compounds. They help predict bioactivity by quantifying structural and molecular features.
    """)
        desc = pd.read_csv('descriptors_output.csv')
        st.write(desc)
        st.write(desc.shape)

        # Read descriptor list used in previously built model
        st.header('**Feature Selection Using Variance Threshold**')
        st.markdown(""" 
           We use Variance Threshold to eliminate features that carry very little useful information or provide minimal discrimination between data points. We reduce 882 features to 179 features.
    """)
        dat = pd.read_csv('descriptor_data2.csv')
        dat = dat.set_index('Name')
        desc_li = desc['Name'].values.tolist()
        dat = dat.loc[desc_li,:]
        st.write(dat)
        #dat = dat.iloc[:,:-4]
        st.write(dat.shape)
        #desc_subset = desc.loc[:,dat.columns]
        #st.write(desc_subset)
        #st.write(desc_subset.shape)

        # Apply trained model to make prediction on query compounds
        build_model(dat)
        

else:
    st.info('Click Sample button in the sidebar to start!')
        

st.markdown("""
**Credits**
- Descriptor calculated using [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) [[Read the Paper]](https://doi.org/10.1002/jcc.21707).
- Project inspired from DataProfessor.
---
""")