"""
Created on Mon Apr 10 19:12:04 2023

@author: Dhruv

For processing experimental data - making graphs

# V1
Use Plotly for advanced ploting
Optimize zero_extensometer()

"""
print("Importing packages stage")
# Import libraries
import os
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
import itertools
import matplotlib.pyplot as plt
from loess.loess_1d import loess_1d
import plotly.express as px
import plotly.graph_objects as go
from statsmodels.nonparametric.smoothers_lowess import lowess
print("Importing packages stage finish") 

""" 
The comments are updated to provide clearer explanations of what each section of the code does.  
"""

class Specimen:
    def __init__(self, home_dir, material_type, condition_type, notch_type, specimen_name, extensometer_plot):
        '''
        Parameters
        ----------
        home_dir : str
            Global folder for saving data.
        material_type : str
            WP1_BM1, ......
        condition_type : str
            AIR, H2_HC1, H2_HC2
        notch_type : str
            SRB, NRB_R2, NRB_R6
        specimen_name : str
            specimen name
        extensometer : str
            Which extensometer has full data

        Returns
        -------
        Object Specimen
        '''
        print('Creating object for '+specimen_name)
        self.home_dir = home_dir
        self.material_type = material_type
        self.condition_type = condition_type
        self.notch_type = notch_type
        self.specimen_name = specimen_name
        self.extensometer_plot = extensometer_plot
        print(self.extensometer_plot)
        if extensometer_plot == 'A':
            self.extensometer_plot = 'Analog In 1'
        else:
            self.extensometer_plot = 'Analog In 2'
        print("Home dir", home_dir)
        print("Condition type", condition_type)
        print("Specimen name", specimen_name)
        self.mts_data_file = f"{home_dir}/{condition_type}/MTS/{specimen_name}/specimen.dat"
        
    def read_csv(self):
        '''
        Returns
        -------
        Creates a self Dataframe containing raw data
        '''
        self.data_df = pd.read_csv(self.mts_data_file, sep='\t', skiprows=[0,1,2,4])
        print(self.specimen_name + ' reading done')
        
    def load_clip(self, load_col='ESH B Force', clip_value=1):
        '''
        Parameters
        ----------
        load_col : str, optional
            Which column has the load value. The default is 'ESH B Force'.
        clip_value : integer, optional
            The value of load (kN) below which data is to ignored. The default is 1 kN.

        Returns
        -------
        Clipped dataframe.
        '''
        index_flag = int(0.5 * len(self.data_df))
        top_df = self.data_df[(self.data_df.index <= index_flag)]
        check_df = self.data_df[(self.data_df.index > index_flag) & (self.data_df.index <= int(len(self.data_df)))]
        eval_df = check_df[check_df[load_col] >= clip_value]
        self.clipped_df = pd.concat([top_df, eval_df])
        print('Load clip executed for ' + self.specimen_name)
        
    def loess_smooth(self, frac=0.05):
        """
        Apply Locally Weighted Scatterplot Smoothing (LOWESS) to smooth the force data using an extensometer as the predictor.

        Parameters
        ----------
        frac : float, optional
            The fraction of data points to use for the local regression, also known as the smoothing factor.
            Default is 0.05, which corresponds to 5% of the data points.
    
        Returns
        -------
        None

        Notes
        -----
        The LOESS method is a non-parametric regression technique used to smooth noisy data. It fits a weighted regression model
        locally to a subset of data points around each point of interest, resulting in a smoothed curve. In this method, the
        'ESH B Force' data is smoothed using the values from the extensometer plot as the predictor. The resulting smoothed
        load data is stored as a new column in the clipped DataFrame, with the column name 'Smoothed load'.
        """
        
        # converting to numpy arrays
        np_load = self.clipped_df['ESH B Force'].values
        np_extensometer = self.clipped_df[self.extensometer_plot].values

        xout, smoothed_load, wout = lowess(np_load, np_extensometer, frac=frac, return_weights=True)
    
        self.clipped_df['Smoothed load'] = smoothed_load
              
    def zero_extensometer(self):
        """ 
        Zero the extensometer data relative to the first data point.

        This method adjusts the values of the extensometer data in the clipped DataFrame relative to the first data point.
        If the first data point is positive, all subsequent values are subtracted by the first value. If the first data point
        is negative, all subsequent values are added by the absolute value of the first value.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        The extensometer data is adjusted to have a reference point of zero, relative to the first data point. This method modifies the extensometer data in place, without creating a new DataFrame. The name of the specimen being processed is printed as part of the output message for reference. 
        """
        
        extensometer_col = self.clipped_df[self.extensometer_plot]
        extensometer_first_val = extensometer_col.iloc[0]
    
        if extensometer_first_val >= 0:
         self.clipped_df[self.extensometer_plot] -= extensometer_first_val
        elif extensometer_first_val < 0:
         self.clipped_df[self.extensometer_plot] += abs(extensometer_first_val)
    
        print('Extensometer zeroed for ' + self.specimen_name)


    def plotly_line_plot(self, linewidth, linestyle, color):
        '''
        Create a line plot using Plotly with advanced features.

        Parameters
        ----------
        linewidth : int
            Line width of the plot.
        linestyle : str
            Line style of the plot.
        color : str
            Color of the plot.

        Returns
        -------
        None
        '''
        
        # Create a figure
        fig = go.Figure()

        # Add a line trace to the figure
        fig.add_trace(go.Scatter(
            x=self.clipped_df[self.extensometer_plot],
            y=self.clipped_df['Smoothed load'],
            mode='lines',
            line=dict(width=linewidth, dash=linestyle, color=color),
            name=self.material_type + "_" + self.condition_type + "_" + self.notch_type + "_" + self.specimen_name + "_" + str(self.test_duration_min) + " min"
        ))

        # Update the figure layout
        fig.update_layout(
            title=self.material_type + ": " + self.notch_type + ": Force-displacement",
            title_font=dict(size=12, color='black'),
            title_pad=dict(t=20),
            xaxis_title=self.extensometer_plot,
            yaxis_title='Smoothed load',
            showlegend=True,
            legend=dict(font=dict(size=10)),
            margin=dict(t=50, r=20, b=50, l=50),
            plot_bgcolor='white'
        )

        # Show the figure
        fig.show()


    def sav_gol_smooth(self, window_size=101):
        '''
        Apply Savitzky-Golay smoothing to 'ESH B Force' column in the clipped dataframe.

        Parameters
        ----------
        window_size : int, optional (default=101)
            Window size for Savitzky-Golay smoothing.

        Returns
        -------
        None
        '''
        # Convert window_size to integer
        window_size = int(window_size)

        # Extract 'ESH B Force' column as numpy array
        np_load = self.clipped_df['ESH B Force'].to_numpy()

        # Apply Savitzky-Golay smoothing to 'ESH B Force' column
        self.clipped_df['Smoothed load'] = savgol_filter(np_load, window_size, 2)

    def testing_time(self):
        self.test_duration_sec = round(self.clipped_df['Time'].iloc[-1],2)
        self.test_duration_min = round(self.test_duration_sec/60,2)

# Example usage of loess_smooth function
#specimen1 = Specimen(home_dir, material_type, condition_type, notch_type, specimen_name, extensometer)
#specimen1.loess_smooth()

# Reading data stage 
  
def process_specimens(specimens_df, home_dir):
    """
    Process specimens data from a DataFrame and create Specimen objects.

    Args:
        specimens_df (pd.DataFrame): DataFrame containing specimens data.
        home_dir (str): Directory path for Specimen objects.

    Returns:
        set: A set of Specimen objects.
    """
    # Create an empty set to store Specimen objects
    specimens = set()

    # Loop through each row in specimens_df
    for _, row in specimens_df.iterrows():
        # Extract the required column values
        args_dict = {'home_dir': home_dir,
                     'material_type': row['material_type'],
                     'condition_type': row['condition_type'],
                     'notch_type': row['notch_type'],
                     'specimen_name': row['specimen_name'],
                     'extensometer_plot': row['extensometer_plot']}

        # Create a Specimen object and add it to the set
        value_i = Specimen(**args_dict)
        specimens.add(value_i)

    # Read in data for each Specimen object in the set
    for specimen in specimens:
        specimen.read_csv()

    return specimens

# Clip load values stage

def clip_load(specimens, load_col, clip_value):
        """
        Clips the load values of Specimen objects in a set to a specified threshold value.

        Args:
            specimens (set): Set of Specimen objects to clip the load values for.
            load_col (str): Name of the column containing the load values in the Specimen objects.
            clip_value (float): Threshold value for clipping the load values.

        Returns:
            None
        """
        for specimen in specimens:
            specimen.load_clip(load_col, clip_value)

# Extensometer stage 

def zero_extensometers(specimens):
    """
    Zeros the extensometer values of Specimen objects in a set.

    Args:
        specimens (set): Set of Specimen objects to zero the extensometer values for.

    Returns:
        None
    """
    for specimen in specimens:
        specimen.zero_extensometer()

    # Calculate total testing time
    for specimen in specimens:
        specimen.testing_time()

    # Apply Smoothing
    o = 0
    for specimen in specimens:
        o = o+1
        print("Applying smoothing",o)
        specimen.sav_gol_smooth(181)


# Plot stage
# One of the following dash styles:
#         ['solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot']
#   - A string containing a dash length list in pixels or percentages
#         (e.g. '5px 10px 2px 2px', '5, 10, 2, 2', '10% 20% 40%', etc.)

def plot_with_plotly(specimens):
    # Define colors and linestyle of curves depending on condition type
    color_dict = {
        'AIR': '#1e64c8',
        'H2_HC1': '#d81e56',
        'H2_HC2': 'royalblue'
    }
    linestyle_dict = {
        'AIR': 'solid',
        'H2_HC1': 'dashdot',
        'H2_HC2': 'dash'
    }

    # Create plotly figure
    fig = go.Figure()

    for specimen in specimens:
        color = color_dict.get(specimen.condition_type, 'blue')
        linestyle = linestyle_dict.get(specimen.condition_type, 'solid')
        fig.add_trace(go.Scatter(
            x=specimen.clipped_df[specimen.extensometer_plot],
            y=specimen.clipped_df['ESH B Force'],
            mode='lines',
            name=f"{specimen.material_type}_{specimen.condition_type}_{specimen.notch_type}_{specimen.specimen_name}",
                line=dict(color=color, dash=linestyle)
        ))

    # Update layout
    fig.update_layout(
        legend=dict(
            xanchor="right",
            yanchor="top",
            x=1,
            y=-0.2,
            bgcolor='white',
            font=dict(size=22),
        ),
        xaxis_title='Elongation [mm]',
        xaxis_title_font=dict(size=22),
        yaxis_title='Tensile force [kN]',
        yaxis_title_font=dict(size=22),
        xaxis=dict(
            tickmode = 'linear',
            tick0 = 0.5,
            tickfont=dict(size=22)),
        yaxis=dict(tickfont=dict(size=22)),
        showlegend=True,
        margin=dict(l=0, r=0, t=0, b=0),
        paper_bgcolor='white',
        plot_bgcolor='white',
    )

    return fig.show()
#############################################################################################################################################
# INPUT PARAMETERS TO BE GIVEN FOR VARIOUS FUNCTIONS

print("Input filename stage") 
# Input data: Specify the file path and sheet name for the Excel file containing the data
# CHANGE
#home_dir = r"S:/shares/flx_lsms_hydrogen/EXPERIMENTAL_DATA/LABO_SOETE/WP1_TENSILE"
home_dir = r"C:/Users/Dhruv/Documents/Jubica_work/shares/flx_lsms_hydrogen/EXPERIMENTAL_DATA/LABO_SOETE/WP1_TENSILE"

excel_file_name = 'Input_Graphs.xlsx'
#excel_sheet_name = '20230322'
#print("Home directory is selected as",home_dir)
print("Input filename stage finish") 


print("Read excel stage")
# Read the Excel file into a DataFrame
# CHANGE
#specimens_excel_df = pd.read_excel(home_dir + '\\' + excel_file_name, sheet_name=excel_sheet_name)
specimens_excel_df_SRB = pd.read_excel('InputGraphsSRB.xlsx')
specimens_excel_df_R2 = pd.read_excel('InputGraphsR2.xlsx')
specimens_excel_df_R6 = pd.read_excel('InputGraphsR6.xlsx')


print("\n Specimens SRB\n",specimens_excel_df_SRB)
#print("\n Specimens R2",specimens_excel_df_R2)
#print("\n Specimens R6",specimens_excel_df_R6)

print("Read excel stage finish") # For debugging


# Set up graph layout: Specify titles and font sizes for the graph
Title_x = 'Elongation [mm]'
Title_y = 'Tensile force [kN]'

print("Reading data stage - process specimens") 
# Function inputs for process_specimens() function
# Extracting required columns from specimens_excel_df
specimens_df_SRB = specimens_excel_df_SRB[['specimen_name', 'material_type', 'notch_type', 'condition_type', 'extensometer_plot']]

#print("\nSRB dataframe\n",specimens_df_SRB) 
specimens_df_R2 = specimens_excel_df_R2[['specimen_name', 'material_type', 'notch_type', 'condition_type', 'extensometer_plot']]
specimens_df_R6 = specimens_excel_df_R6[['specimen_name', 'material_type', 'notch_type', 'condition_type', 'extensometer_plot']]

#home_dir = r"S:/shares/flx_lsms_hydrogen/EXPERIMENTAL_DATA/LABO_SOETE/WP1_TENSILE"

# Function call for process_specimens() function
specimens_SRB = process_specimens(specimens_df_SRB,home_dir)
specimens_R2 = process_specimens(specimens_df_R2,home_dir)
specimens_R6 = process_specimens(specimens_df_R6,home_dir)

print("Reading data stage finish - specimens processed") 

print("Clip load values stage start") 
# Function inputs for clip_load() function
# Clip data to remove data at the end of the curve
load_col = 'ESH B Force'
clip_value = 1
# specimens in Line 447 -> value is process_specimens(specimens_df,home_dir)
# Function call for clip_load() function
clip_load(specimens_SRB,load_col,clip_value) 
clip_load(specimens_R2,load_col,clip_value)
clip_load(specimens_R6,load_col,clip_value)

print("Clip load values stage finish") 


print("Extensometer stage start") 
# Function inputs for zero_extensometers() function
# specimens in Line 447 -> value is process_specimens(specimens_df,home_dir)

# Function call for zero_extensometers() function
zero_extensometers(specimens_SRB) #fine till here
zero_extensometers(specimens_R2)
zero_extensometers(specimens_R6)

print("Extensometer stage finish") 


print("Plot stage")
# Function inputs for plot_with_plotly() function
# specimens in Line 447 -> value is process_specimens(specimens_df,home_dir)
Size_title = 32
Size_title_x = 28
Size_title_y = 28
Size_axis_x = 22
Size_axis_y = 22

# Function call for plot_with_plotly() function
plot_with_plotly(specimens_SRB)
plot_with_plotly(specimens_R2)
plot_with_plotly(specimens_R6)