# GUI-for-stride-by-stride-EMG-visualisation-or-analysis

This GUI can be used to visualize electromyography (EMG) data. It is meant for EMG data collected during gait for which you have the gait events. This GUI can be used to visualize electromyography (EMG) data. It is meant for EMG data collected during gait for which you have the gait events. Based on e.g. the left heel strikes, one can visualize the data per stride or per several strides and inspect the data. The GUI also contains rectification and filtering options, as well as the function to replace the EMG values of (a) stride(s) by NaNs. This can be convenient if only a few strides contain noise clearly observable during visual inspection (e.g. for me this happened because one of the EMG wires was jammed). You can save the processed EMG data as well as the indices of any NaN values by pressing "save EMG" before closing the GUI.

Use:
- Use the 'GUI_EMG.m' function in a main script
- Use the four inputs as described in the function's help in the input order as described

Note:
The main purpose of this GUI is data visualisation. Analysis by means of the GUI is not transparent and therefore not recommended.
