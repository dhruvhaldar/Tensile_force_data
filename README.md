# Tensile_force_data using Plotly
1. Images are generated and saved as png in the same folder.
2. Go through the docstrings to know what to change (mostly trace for changing dash and color)

## plot_with_plotly() function
1. This sets the color variable to the value associated with the condition_type key in the color_dict dictionary for the current specimen. If there is no value associated with the condition_type key, the default value 'blue' is used.
2. This sets the linestyle variable to the value associated with the condition_type key in the linestyle_dict dictionary for the current specimen. If there is no value associated with the condition_type key, the default value 'solid' is used.
3. This adds a new trace to the fig figure object, using the go.Scatter() constructor. 
4. This sets the x-values of the trace to the extensometer_plot column of the clipped_df DataFrame for the current specimen.  
5. This sets the y-values of the trace to the 'ESH B Force' column of the clipped_df DataFrame for the current specimen. 
6. This sets the mode of the trace to 'lines', indicating that it should be plotted as a line.   
7. This sets the name of the trace to a formatted string that includes the material_type, condition_type, notch_type, and specimen_name attributes of the current specimen.   
8. This sets the line color and dash style of the trace using the color and linestyle variables, respectively.  
9. The resulting code creates a figure with a separate trace for each specimen in the specimens list, where each trace represents a line plot of the #'ESH B Force' column of the clipped_df DataFrame against the extensometer_plot column. The color and dash style of each line are determined by the # color_dict and linestyle_dict dictionaries based on the condition_type attribute of each specimen.
