import tkinter as tk
from tkinter import filedialog
import os.path
from readerAndLoader import ReaderAndLoader
from calculateAndPlot import DistancePlotter

DEFAULT_FILE_PREFIX = "MIL68Ga"
DEFAULT_DIRECTORY = "/home/mfi/Desktop/mfi/MIL-68Ga-guest"
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
    plot_button.config(state=tk.NORMAL)  # Enable the plot button

def browse_directory():
    directory = filedialog.askdirectory()
    directory_entry.delete(0, tk.END)
    directory_entry.insert(tk.END, directory)
    check_saved_data()

def check_saved_data():
    file_path = os.path.join(directory_entry.get(), "sourcedata.npy")
    if os.path.isfile(file_path):
        load_saved_data_button.config(state=tk.NORMAL)
        saved_data_label.config(text="Saved data found!", fg="green")
        plot_button.config(state=tk.DISABLED)  # Disable the plot button
    else:
        load_saved_data_button.config(state=tk.DISABLED)
        saved_data_label.config(text="No saved data found.", fg="black")
        plot_button.config(state=tk.DISABLED)  # Disable the plot button

def plot_distances():
    global distancePlotter
    global centroid_vars
    #data, numAtoms, centroid_indices_flat, guest_indices, numOfCentroids, lattice_values = ReaderAndLoader.load_data(directory_entry.get(), index_entry.get(), os.path.join(directory_entry.get(), "sourcedata.npy"))
    #distancePlotter = DistancePlotter(merged_df_values=data, numAtoms=numAtoms, centroid_indices_j=centroid_indices_flat, guest_indices=guest_indices, numOfCentroids=numOfCentroids, lattice_values=lattice_values)
    stepsToPlot = int(stepsToPlot_entry.get())
    print("centroid vards=" ,centroid_vars)
    distancePlotter.plot_distances(stepsToPlot, centroid_vars)
    

# Create the main window
window = tk.Tk()
window.title("XYZ File Loader")

# Create the directory label and entry
directory_label = tk.Label(window, text="Directory:")
directory_label.pack()
directory_entry = tk.Entry(window)
directory_entry.insert(tk.END, DEFAULT_DIRECTORY)  # Set default directory
directory_entry.pack()

# Create the browse button
browse_button = tk.Button(window, text="Browse", command=browse_directory)
browse_button.pack()

# Create the prefix label and entry
prefix_label = tk.Label(window, text="File Prefix:")
prefix_label.pack()
prefix_entry = tk.Entry(window)
prefix_entry.insert(tk.END, DEFAULT_FILE_PREFIX)  # Set default file prefix
prefix_entry.pack()

# Create the index filename label and entry
index_label = tk.Label(window, text="Index Filename:")
index_label.pack()
index_entry = tk.Entry(window)
index_entry.insert(tk.END, DEFAULT_INDEX_FILENAME)  # Set default index filename
index_entry.pack()

# Create the load button
load_button = tk.Button(window, text="Read and save XYZ Files", command=read_and_save_data)
load_button.pack()

# Start the main event loop
def load_saved_data():
    global distancePlotter
    directory = directory_entry.get()
    index_filename = index_entry.get()
    file_path = os.path.join(directory_entry.get(), "sourcedata.npy")
    
    data, numAtoms, centroid_indices_flat, guest_indices, numOfCentroids, lattice_values, centroid_indices = ReaderAndLoader.load_data(directory, index_filename, file_path)  # Call the function to load the saved data
    distancePlotter = DistancePlotter(merged_df_values=data, numAtoms=numAtoms, centroid_indices_j=centroid_indices, guest_indices=guest_indices, numOfCentroids=numOfCentroids, lattice_values=lattice_values)
    
    num_atoms_label.config(text=f"Number of Atoms: {numAtoms}")
    num_centroids_label.config(text=f"Number of Centroids: {numOfCentroids}")
    lattice_params_label.config(text=f"Lattice Parameters: {lattice_values}")
    plot_button.config(state=tk.NORMAL)  # Enable the plot button

    # Create tick boxes for each centroid
    global centroid_vars
    centroid_vars = []
    for i in range(numOfCentroids):
        var = tk.BooleanVar()
        centroid_vars.append(var)
        centroid_checkbox = tk.Checkbutton(window, text=f"Centroid {i+1}", variable=var)
        centroid_checkbox.pack()

# Create the load saved data button
load_saved_data_button = tk.Button(window, text="Load Saved Data", command=load_saved_data, state=tk.DISABLED)
load_saved_data_button.pack()

# Create the stepstoplot entry
stepsToPlot_label = tk.Label(window, text="Steps to plot:")
stepsToPlot_label.pack()
stepsToPlot_entry = tk.Entry(window)
stepsToPlot_entry.pack()

# Create the saved data label
saved_data_label = tk.Label(window, text="", fg="black")
saved_data_label.pack()

# Create the labels for numAtoms, number of centroids and lattice parameters
num_atoms_label = tk.Label(window, text="Number of Atoms: ")
num_atoms_label.pack()

num_centroids_label = tk.Label(window, text="Number of Centroids: ")
num_centroids_label.pack()

lattice_params_label = tk.Label(window, text="Lattice Parameters: ")
lattice_params_label.pack()

# Create the plot button
plot_button = tk.Button(window, text="Plot", state=tk.DISABLED, command = plot_distances)
plot_button.pack()

# Check if the default directory already contains the "sourcedata" file
check_saved_data()

# Create the main loop
window.mainloop()

window.mainloop()
