"""
@author: jubica

Graphs

"""
print("Importing packages stage") # For debugging
# Import libraries
import os
# Clearing the Screen
os.system('cls')
# import re
import numpy as np
import pandas as pd
# import sys
from scipy.signal import savgol_filter
# from scipy import interpolate
# import warnings
# from math import pi
import itertools

import plotly.express as px

# Plot related packages
import matplotlib.pyplot as plt

# from FLUXYS_experimental_data_FUNCTIONS import load_clip
from loess.loess_1d import loess_1d
print("Importing packages stage finish") # For debugging


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
Size_title   = 12 # 32
Size_title_x = 12 # 28
Size_title_y = 12 # 28
Size_axis_x  = 12 # 22
Size_axis_y  = 12 # 22

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
        Clipped dataframe.

        '''
        index_flag = int(0.5*len(self.data_df))
        top_df = self.data_df[(self.data_df.index<=index_flag)]
        check_df = self.data_df[(self.data_df.index>index_flag) & (self.data_df.index<=int(len(self.data_df)))]
        eval_df = check_df[check_df[load_col]>=clip_value]
        self.clipped_df = pd.concat([top_df,eval_df])
        print('Load clip executed for '+self.specimen_name)
        #######################################
        
    def loess_smooth(self, frac=0.05):
        #TODO: Fix this- this is taking too long. ☺
        #converting to numpy arrays
        np_load = self.clipped_df['ESH B Force'].to_numpy()
        np_extensometer = self.clipped_df[self.extensometer_plot].to_numpy()

        xout, self.clipped_df['Smoothed load'], wout = loess_1d(np_extensometer, np_load, degree=1, frac=0.05)
        
    def zero_extensometer(self):
        if self.clipped_df[self.extensometer_plot][0] >= 0:
            self.clipped_df[self.extensometer_plot] = self.clipped_df[self.extensometer_plot]-self.clipped_df[self.extensometer_plot][0]
        elif self.clipped_df[self.extensometer_plot][0] < 0:
            self.clipped_df[self.extensometer_plot] = self.clipped_df[self.extensometer_plot]+abs(self.clipped_df[self.extensometer_plot][0])
        
        print('Extensometer zeroed for '+self.specimen_name)
    
    def basic_plt_line_plot(self, linewidth, linestyle, color):
        '''      
        Parameters
        ----------
        linewidth : int
            line width of plot.
        linestyle : str
            line style of plot.
        color : str
            colour of plot.

        Returns
        -------
        Adds plot to the matplotlib figure.
        '''
        # ax.plot(self.clipped_df[self.extensometer_plot], self.clipped_df['ESH B Force'],color = color, lw=linewidth, ls=linestyle, label = self.material_type+"_"+self.condition_type+"_"+self.notch_type+"_"+self.specimen_name+"_"+str(self.test_duration_min)+" min")
        ax.plot(self.clipped_df[self.extensometer_plot], self.clipped_df['Smoothed load'],color = color, lw=1.5, ls=linestyle, label = self.material_type+"_"+self.condition_type+"_"+self.notch_type+"_"+self.specimen_name+"_"+str(self.test_duration_min)+" min")
        ax.set_title(self.material_type+": "+self.notch_type+": Force-displacement",size=12, pad = 20)

    def sav_gol_smooth(self, window_size = 101):
        #converting to numpy arrays
        window_size = int(window_size)
        np_load = self.clipped_df['ESH B Force'].to_numpy()
        # np_extensometer = self.clipped_df[self.extensometer_plot].to_numpy()

        self.clipped_df['Smoothed load'] = savgol_filter(np_load, window_size, 2)

    def testing_time(self):
        self.test_duration_sec = round(self.clipped_df['Time'].iloc[-1],2)
        self.test_duration_min = round(self.test_duration_sec/60,2)
    
# Reading in data
# Creating another df with only the columns that are required for further analysis.    
print("Reading data stage") # For debugging
specimens_df = specimens_excel_df[['specimen_name', 'material_type', 'notch_type', 'condition_type', 'extensometer_plot']]
specimens = set()

#print("Specimens df",specimens_df) # For debugging
#print("specimens",specimens) # for debugging

#Reading in each of the specimens from the dictionary
for i in range(len(specimens_df)):
    specimen_i = specimens_df.loc[i, 'specimen_name']
    material_i = specimens_df.loc[i,'material_type']
    condition_i = specimens_df.loc[i,'condition_type']
    notch_i = specimens_df.loc[i,'notch_type']
    extensometer_i = specimens_df.loc[i, 'extensometer_plot']
    value_i = Specimen(home_dir, material_i, condition_i, notch_i, specimen_i, extensometer_i)
    specimens.add(value_i)

for specimen in specimens:
    specimen.read_csv() 
    
print("Reading data stage finish") # For debugging

# Remove load below a certain value
print("Clip data stage") # For debugging

## Clip data to remove data at the end of the curve
load_col = 'ESH B Force'
clip_value = 1

for specimen in specimens:
    specimen.load_clip(load_col,clip_value)
print("Clip data stage finish") # For debugging

print("Extensometer stage") # For debugging
# Zero extensometers
for specimen in specimens:
    specimen.zero_extensometer()

print("Extensometer stage finish") # For debugging

print("Calculate total testing time stage") # For debugging
# Calculate total testing time
for specimen in specimens:
    specimen.testing_time()
print("Calculate total testing time stage finish") # For debugging

print("Smoothing stage") # For debugging

# Apply Smoothing
for specimen in specimens:
    specimen.sav_gol_smooth(181)
print("Smoothing stage finish") # For debugging


print("Plot stage") # For debugging

# Plot
# list of markers to iterate

linestyles = itertools.cycle(('solid', 'dotted', 'dashed'))

# fig, ax = plt.subplots(figsize=(12,8))
fig, ax = plt.subplots(figsize=(6,4))      # VIKI

#Define colors and linestyle of curves depending on condition type
air_count = 0
hc1_count = 0
hc2_count = 0
linestyles_list = ['solid', 'dashed']
#linestyles_list = ['solid', 'dotted', 'dashed']
for specimen in specimens:
    if specimen.condition_type == 'AIR':
        color = '#1e64c8'
        # color = 'red'   # VIKI
        linestyle = linestyles_list[air_count]
        air_count += 1 
    elif specimen.condition_type == 'H2_HC1':
        color = '#d81e56'
        # color = 'navy'  # VIKI
        linestyle = linestyles_list[hc1_count]
        hc1_count += 1 
    #elif specimen.condition_type == 'H2_HC2':
        #color = 'royalblue'
        #linestyle = linestyles_list[hc2_count]
        #hc2_count += 1 
    else:
        color = 'pink'
        print('Condition type not known: ' + specimen.condition_type)
    specimen.basic_plt_line_plot(3, linestyle, color)


# ax.plot(specimen_1.clipped_df[specimen_1.extensometer_plot], specimen_1.clipped_df['ESH B Force'],color = "black", lw=4, ls='solid', label = specimen_1.material_type+"_"+specimen_1.condition_type+"_"+specimen_1.notch_type+"_"+specimen_1.specimen_name)
# ax.plot(clipped_df[extensometer_plot], clipped_df['Smoothed load'],color = "red", lw=2, ls='solid', label = 'Smoothed')

numbercond = 0
for cond in [air_count, hc1_count, hc2_count]:
    if cond > 0:
        numbercond += 1

# ax.legend(loc='best', bbox_to_anchor=(1, -0.2), fancybox=True, shadow=True, ncol=numbercond, fontsize=9) #location of legend
ax.legend(loc='best', fontsize=12) #location of legend (VIKI)
ax.set_xlabel('Elongation [mm]',size=Size_title_x)
ax.set_ylabel('Tensile force [kN]',size=Size_title_y) #link for scientific notations below
# ax1.set_title(specimen_id+" Validation Attempt: "+simulation_id)

ax.set_facecolor('white')

# Define axis limits depending on notch type
if specimen.notch_type == 'SRB':
    ax.set_ylim([0,20])
    ax.set_xlim([0,8])
else: 
    ax.set_ylim([0,30])
    ax.set_xlim([0,3])
    
ax.xaxis.set_tick_params(labelsize=Size_axis_x)
ax.yaxis.set_tick_params(labelsize=Size_axis_y)

plt.tight_layout()
plt.grid(visible=False)

plt.savefig('SRB.svg', dpi=1200)  

# Show Plot
#plt.show()


#print("Plot stage finish") # For debugging


def plot_force_displacement(specimen):
    '''
    Parameters
    ----------
    specimen : Specimen
        The Specimen object for which the graph is to be plotted.

    Returns
    -------
    None.

    '''
    # Create a scatter plot with the extensometer displacement on the x-axis and the force on the y-axis
    fig = px.scatter(specimen.clipped_df, x=specimen.extensometer_plot, y='ESH B Force')

    # Set the graph title and axis labels
    fig.update_layout(
        title=specimen.specimen_name,
        xaxis_title='Elongation [mm]',
        yaxis_title='Tensile force [kN]'
    )

    # Show the graph
    fig.show()