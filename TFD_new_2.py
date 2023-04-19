# Import libraries
import os
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from loess.loess_1d import loess_1d
import itertools
import plotly.graph_objects as go

# Clearing the Screen
os.system('cls')

print("Input filename stage") # For debugging
# Input data
#Where to find the raw data and which specimens you want to plot?
home_dir = r"S:/shares/flx_lsms_hydrogen/EXPERIMENTAL_DATA/LABO_SOETE/WP1_TENSILE"
print("Home directory is selected as",home_dir)
print("Input filename stage finish") # For debugging

print("Read excel stage") # For debugging
#excel_file_name = '/Input_GraphsJ.xlsx' #it is assumed it is in the home_dir
#excel_sheet_name = '20230322' #sheetname of the excel
specimens_excel_df = pd.read_excel('Input_GraphsJ.xlsx')
#print (specimens_excel_df) # For debugging
print("Read excel stage finish") # For debugging

#Lay-out graph
Title_x = 'Elongation [mm]'
Title_y = 'Tensile force [kN]'
Size_title = 32
Size_title_x = 28
Size_title_y = 28
Size_axis_x = 22
Size_axis_y = 22

# Class of specimens
class Specimen:
    def __init__(self, home_dir, material_type, condition_type, notch_type, specimen_name, extensometer):
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
        self.extensometer = extensometer
        if extensometer == 'A':
            self.extensometer_plot = 'Analog In 1'
        else:
            self.extensometer_plot = 'Analog In 2'
        self.mts_data_file = home_dir+"/"+condition_type+"/MTS/"+specimen_name+"/specimen.dat"
        print(self.mts_data_file) # For debugging
        #######################################

    def read_csv(self):
        '''
        Returns
        -------
        Creates a self Dataframe containing raw data

        '''
        self.data_df = pd.read_csv(self.mts_data_file, sep='\t', skiprows=[0,1,2,4])
        print(self.specimen_name+' reading done')
        #######################################

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
        None

        '''
        # Filter data based on clip_value
        self.load_col = load_col
        self.clip_value = clip_value
        self.filtered_data_df = self.data_df[self.data_df[load_col] >= clip_value]
        print('Data filtered for '+self.specimen_name)
        #######################################

    def smooth_data(self, window=25, polyorder=3):
        '''
        Parameters
        ----------
        window : int, optional
            Window size for Savitzky-Golay filter. The default is 25.
        polyorder : int, optional
            Polynomial order for Savitzky-Golay filter. The default is 3.

        Returns
        -------
        None

        '''
        # Apply Savitzky-Golay filter to smooth data
        self.window = window
        self.polyorder = polyorder
        self.filtered_data_df['Smoothed Load'] = savgol_filter(self.filtered_data_df[self.load_col], window, polyorder)
        print('Data smoothed for '+self.specimen_name)
        #######################################

    def loess_regression(self, frac=0.1):
        '''
        Parameters
        ----------
        frac : float, optional
            Fraction of data for LOESS regression. The default is 0.1.

        Returns
        -------
        None

        '''
        # Perform LOESS regression to obtain a trendline
        self.frac = frac
        loess = loess_1d(self.filtered_data_df['Smoothed Load'].values, self.filtered_data_df['Elongation [mm]'].values, frac=frac)
        self.trendline = loess.eval(self.filtered_data_df['Elongation [mm]'].values)
        print('LOESS regression performed for '+self.specimen_name)
        #######################################

    def plot_data(self):
        '''
        Returns
        -------
        None

        '''
        # Create plot using Plotly
        fig = go.Figure()

        # Add raw data
        fig.add_trace(go.Scatter(x=self.data_df['Elongation [mm]'], y=self.data_df[self.load_col],
                                 mode='lines',
                                 name='Raw Data',
                                 line=dict(color='blue')))

        # Add filtered data
        fig.add_trace(go.Scatter(x=self.filtered_data_df['Elongation [mm]'], y=self.filtered_data_df[self.load_col],
                                 mode='lines',
                                 name='Filtered Data',
                                 line=dict(color='green')))

        # Add smoothed data
        fig.add_trace(go.Scatter(x=self.filtered_data_df['Elongation [mm]'], y=self.filtered_data_df['Smoothed Load'],
                                 mode='lines',
                                 name='Smoothed Data',
                                 line=dict(color='orange')))

        # Add trendline
        fig.add_trace(go.Scatter(x=self.filtered_data_df['Elongation [mm]'], y=self.trendline,
                                 mode='lines',
                                 name='Trendline',
                                 line=dict(color='red')))

        # Update plot layout
        fig.update_layout(
            title=self.specimen_name,
            xaxis_title=Title_x,
            yaxis_title=Title_y,
            title_font=dict(size=Size_title),
            xaxis=dict(title_font=dict(size=Size_title_x), tickfont=dict(size=Size_axis_x)),
            yaxis=dict(title_font=dict(size=Size_title_y), tickfont=dict(size=Size_axis_y)),
            showlegend=True
        )

        # Show plot
        fig.show()
        #################################################################################################

def plot_all_data():
    # List of specimen names
    specimens = ['SRB', 'R6', 'R2']

    ## Clip data to remove data at the end of the curve
    load_col = 'ESH B Force'
    clip_value = 1
    
    # Loop through all specimens in the excel sheet and plot data
    for specimen in specimens:
        spec = Specimen(specimen)
        spec.load_data(home_dir)
        specimen.load_clip(load_col,clip_value)
        specimen.zero_extensometer()
        specimen.testing_time()
        specimen.sav_gol_smooth(181)
        spec.filter_data()
        spec.smooth_data()
        spec.loess_regression()
        spec.plot_data()

    # Call the function to plot data for all specimens

    plot_all_data()