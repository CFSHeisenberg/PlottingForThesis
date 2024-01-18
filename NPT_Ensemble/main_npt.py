import tkinter as tk
from tkinter import filedialog
import os.path
from readerAndLoader_npt import ReaderAndLoader
from calculateAndPlot_npt import DistancePlotter
from calculateAndPlot_npt import plotFromSavedDistances
import tkinter as tk
from tkinter import *
import os

DEFAULT_FILE_PREFIX = "MOF5"
DEFAULT_DIRECTORY = "/home/mfi/Desktop/mfi/filesManuelOtt"
DEFAULT_INDEX_FILENAME = "indicessorted.dat"  # Set default index filename

ReaderAndLoader = ReaderAndLoader()
distancePlotter = None
centroid_vars = []

def read_and_save_data():
    directory = directory_entry.get()
    file_prefix = prefix_entry.get()
    index_filename = index_entry.get()
    ReaderAndLoader.read_and_save_data(directory, file_prefix, index_filename)
    load_button.config(bg="green")  # Change button color to green
    check_saved_data()
    #plot_button.config(state=tk.NORMAL)  # Enable the plot button

def browse_directory():
    directory = filedialog.askdirectory()
    directory_entry.delete(0, tk.END)
    directory_entry.insert(tk.END, directory)
    check_saved_data()

def browse_distances_file():
    file_path = filedialog.askopenfilename()
    distances_entry.delete(0, tk.END)
    distances_entry.insert(tk.END, file_path)
    check_saved_data_distances()
    #plotfromdistance_button.config(state=tk.NORMAL)  # Enable the plot button

def check_saved_data_distances():
    file_path = distances_entry.get()
    if os.path.isfile(file_path):
        saved_distance_data_label.config(text="Saved distances found!", fg="green")
        plotfromdistance_button.config(state=tk.NORMAL)  # Enable the plot button
        directory = directory_entry.get()
        index_filename = index_entry.get()
        
        numOfCentroids = ReaderAndLoader.getNumCentroids(directory, index_filename)
                    
        # Create tick boxes for each centroid
        global centroid_vars
        centroid_vars = []
        for i in range(numOfCentroids):
            var = tk.BooleanVar()
            centroid_vars.append(var)
            centroid_checkbox = tk.Checkbutton(window, text=f"Centroid {i+1}", variable=var)
            centroid_checkbox.pack()
                    
    else:
        saved_distance_data_label.config(text="No saved distances found.", fg="black")
        plotfromdistance_button.config(state=tk.DISABLED)  # Disable the plot button

def check_saved_data():
    file_path = os.path.join(directory_entry.get(), "sourcedata.npy")
    file_path_lattices = os.path.join(directory_entry.get(), "sourcedata_lattices.npy")
    if os.path.isfile(file_path) and os.path.isfile(file_path_lattices):
        load_saved_data_button.config(state=tk.NORMAL)
        saved_data_label.config(text="Source data and lattice-data found in chosen directory!", fg="green")
        plot_button.config(state=tk.DISABLED)  # Disable the plot button

    elif os.path.isfile(file_path) and not os.path.isfile(file_path_lattices):
        load_saved_data_button.config(state=tk.DISABLED)
        saved_data_label.config(text="No lattice-data (sourcedata_lattices.npy) found in chosen directory. Try reading in the trajectories.", fg="black")
        plot_button.config(state=tk.DISABLED)  # Disable the plot button
        
    elif not os.path.isfile(file_path) and os.path.isfile(file_path_lattices):
        load_saved_data_button.config(state=tk.DISABLED)
        saved_data_label.config(text="No source data (sourcedata.npy) found in chosen directory. Try reading in the trajectories.", fg="black")
        plot_button.config(state=tk.DISABLED)  # Disable the plot button
    else:
        load_saved_data_button.config(state=tk.DISABLED)
        saved_data_label.config(text="No source data (sourcedata.npy) and no lattice-data (sourcedata_lattices.npy) found in chosen directory. Try loading the data using the button above.", fg="black")
        plot_button.config(state=tk.DISABLED)  # Disable the plot button

def plot_distances():
    global distancePlotter
    global centroid_vars
    #data, numAtoms, centroid_indices_flat, guest_indices, numOfCentroids, lattice_values = ReaderAndLoader.load_data(directory_entry.get(), index_entry.get(), os.path.join(directory_entry.get(), "sourcedata.npy"))
    #distancePlotter = DistancePlotter(merged_df_values=data, numAtoms=numAtoms, centroid_indices_j=centroid_indices_flat, guest_indices=guest_indices, numOfCentroids=numOfCentroids, lattice_values=lattice_values)
    stepsToPlot = int(stepsToPlot_entry.get())
    distancePlotter.plot_distances(stepsToPlot, centroid_vars)
    
    # Start the main event loop
def load_saved_data():
    global distancePlotter
    directory = directory_entry.get()
    index_filename = index_entry.get()
    file_path = os.path.join(directory_entry.get(), "sourcedata.npy")
    file_path_lattices = os.path.join(directory_entry.get(), "sourcedata_lattices.npy")
    
    data, numAtoms, centroid_indices_flat, guest_indices, numOfCentroids, lattice_values, centroid_indices, lat_data = ReaderAndLoader.load_data(directory, index_filename, file_path, file_path_lattices)  # Call the function to load the saved data
    distancePlotter = DistancePlotter(merged_df_values=data, numAtoms=numAtoms, centroid_indices_j=centroid_indices, guest_indices=guest_indices, numOfCentroids=numOfCentroids, merged_lattice_values = lat_data)
    
    num_atoms_label.config(text=f"Number of Atoms: {numAtoms}")
    num_centroids_label.config(text=f"Number of Centroids: {numOfCentroids}")
    plot_button.config(state=tk.NORMAL)  # Enable the plot button

    # Create tick boxes for each centroid
    global centroid_vars
    centroid_vars = []
    for i in range(numOfCentroids):
        var = tk.BooleanVar()
        centroid_vars.append(var)
        centroid_checkbox = tk.Checkbutton(window, text=f"Centroid {i+1}", variable=var)
        centroid_checkbox.pack()
        
# Function to ask for directory
def ask_directory():
    directory = filedialog.askdirectory()
    if directory:
        calculate_distances(directory)
        
def plot_from_distances_file():
    file_path_distance = distances_entry.get()
    index_filename = index_entry.get()
    file_path = os.path.join(directory_entry.get(), "sourcedata.npy")
        
    plotFromSavedDistances(centroid_vars, file_path_distance)

# Function to calculate distances and save data
def calculate_distances(save_directory):
    
    steps = int(stepsToCalc_entry.get())
    index_filename = index_entry.get()
    directory = directory_entry.get()
    file_path = os.path.join(directory_entry.get(), "sourcedata.npy")
    file_path_lattices = os.path.join(directory_entry.get(), "sourcedata_lattices.npy")

    
    data, numAtoms, centroid_indices_flat, guest_indices, numOfCentroids, lattice_values, centroid_indices, lat_data = ReaderAndLoader.load_data(directory, index_filename, file_path, file_path_lattices)  # Call the function to load the saved data
    distancePlotter = DistancePlotter(merged_df_values=data, numAtoms=numAtoms, centroid_indices_j=centroid_indices, guest_indices=guest_indices, numOfCentroids=numOfCentroids, merged_lattice_values=lat_data)
    file_path = distancePlotter.calcAndSaveAllDistancesForNSteps(steps, save_directory)
    saved_distance_data_label.config(text=f"Distances saved to: {file_path}")
    

# Create the main window
window = tk.Tk()
window.title("Trajectory Analyzer for centroid distances")

# Create a frame for the directory, prefix, and index filename widgets
input_frame = tk.Frame(window, highlightbackground="red", highlightthickness=4)
input_frame.pack(side="left", padx=20)

# Create a frame for the load and load saved data buttons
button_frame = tk.Frame(window, highlightbackground="blue", highlightthickness=4)
button_frame.pack(side="left", padx=20)

# Create a new frame for the stepstocalc entry and calculate distances button
calc_frame = tk.Frame(window, highlightbackground="green", highlightthickness=4)
calc_frame.pack(side="left", padx=20)

# Create a frame for the load precalculated distance file and plot chosen centroids section
load_distances_frame = tk.Frame(window, highlightbackground="orange", highlightthickness=4)
load_distances_frame.pack(side="left", padx=20)

# Create a frame for the properties of simulation
properties_frame = tk.Frame(window, highlightbackground="purple", highlightthickness=4)
properties_frame.pack(side="left", padx=20)

# Create a label for the frame
label = tk.Label(input_frame, text="MANDATORY: File locations", fg="red")
label.pack()

# Create the directory label and entry
directory_label = tk.Label(input_frame, text="Directory containing XYZ files and indices file:")
directory_label.pack()
directory_entry = tk.Entry(input_frame)
directory_entry.insert(tk.END, DEFAULT_DIRECTORY)  # Set default directory
directory_entry.pack()
directory_entry.config(highlightbackground="red", highlightthickness=1)  # Add red rectangle

# Create the browse button
browse_button = tk.Button(input_frame, text="Browse", command=browse_directory)
browse_button.pack()

# Create the prefix label and entry
prefix_label = tk.Label(input_frame, text="File Prefix:")
prefix_label.pack()
prefix_entry = tk.Entry(input_frame)
prefix_entry.insert(tk.END, DEFAULT_FILE_PREFIX)  # Set default file prefix
prefix_entry.pack()
prefix_entry.config(highlightbackground="red", highlightthickness=1)  # Add red rectangle

# Create the index filename label and entry
index_label = tk.Label(input_frame, text="Index Filename:")
index_label.pack()
index_entry = tk.Entry(input_frame)
index_entry.insert(tk.END, DEFAULT_INDEX_FILENAME)  # Set default index filename
index_entry.pack()
index_entry.config(highlightbackground="red", highlightthickness=1)  # Add red rectangle

# Create a label for the frame
button_label = tk.Label(button_frame, text="Read and save trajectories or load existing data", fg="blue")
button_label.pack()

# Create the load button
load_button = tk.Button(button_frame, text="Read and save XYZ Files and lattices", command=read_and_save_data)
load_button.pack()

# Create the load saved data button
load_saved_data_button = tk.Button(button_frame, text="Load Saved Data", command=load_saved_data, state=tk.DISABLED)
load_saved_data_button.pack()

# Create the saved data label
saved_data_label = tk.Label(button_frame, text="", fg="black")
saved_data_label.pack()

# Create the stepstoplot entry
stepsToPlot_label = tk.Label(button_frame, text="Steps to plot:")
stepsToPlot_label.pack()
stepsToPlot_entry = tk.Entry(button_frame)
stepsToPlot_entry.pack()

# Create the plot button
plot_button = tk.Button(button_frame, text="Plot from loaded data", state=tk.DISABLED, command=plot_distances)
plot_button.pack()

# Create a label for the frame
calc_label = tk.Label(calc_frame, text="Calculate n steps and save results to file for later plotting", fg="green")
calc_label.pack()


# Create the stepstocalc entry
stepsToCalc_label = tk.Label(calc_frame, text="Steps to calc:")
stepsToCalc_label.pack()
stepsToCalc_entry = tk.Entry(calc_frame)
stepsToCalc_entry.pack()

# Create the calculate distances for steps button
calculate_distances_button = tk.Button(calc_frame, text="Calculate Distances for Steps", command=lambda: ask_directory())
calculate_distances_button.pack()

# Create the saved distance data label
saved_distance_data_label = tk.Label(calc_frame, text="", fg="black")
saved_distance_data_label.pack()


# Create a label for the frame
load_distances_label = tk.Label(load_distances_frame, text="Load precalculated distance file and plot chosen centroids", fg="orange")
load_distances_label.pack()

# Create second browse button for distances file
browse_distances_button = tk.Button(load_distances_frame, text="Browse file containing distances", command=browse_distances_file)
browse_distances_button.pack()

# Create the distances file label and entry
distances_label = tk.Label(load_distances_frame, text="Distances File:")
distances_label.pack()
distances_entry = tk.Entry(load_distances_frame)
distances_entry.pack()

# Create the saved distance data label
saved_distance_data_label = tk.Label(load_distances_frame, text="", fg="black")
saved_distance_data_label.pack()

# Create the plotfromdistancefile button
plotfromdistance_button = tk.Button(load_distances_frame, text="Plot from loaded distances", state=tk.DISABLED, command=plot_from_distances_file)
plotfromdistance_button.pack()

# Create a label for the frame
properties_label = tk.Label(properties_frame, text="Properties of simulation", fg="purple")
properties_label.pack()

# Create the labels for numAtoms, number of centroids and lattice parameters
num_atoms_label = tk.Label(properties_frame, text="Number of Atoms: ")
num_atoms_label.pack()

num_centroids_label = tk.Label(properties_frame, text="Number of Centroids: ")
num_centroids_label.pack()

# Check if the default directory already contains the "sourcedata" file
check_saved_data()

# Create the main loop
window.mainloop()

window.mainloop()
