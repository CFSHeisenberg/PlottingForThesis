import tkinter as tk
from tkinter import filedialog
import os.path
from readerAndLoader import ReaderAndLoader
from calculateAndPlot import DistancePlotter
from calculateAndPlot import plotFromSavedDistances
import tkinter as tk
from tkinter import *
import os

DEFAULT_FILE_PREFIX = "MIL68Ga"
DEFAULT_DIRECTORY = "/home/mfi/Desktop/mfi/MIL-68Ga-guest/2ndguest/3rdguest/3GuestData/heidi"
DEFAULT_INDEX_FILENAME = "indicessorted.dat"  # Set default index filename
#DEFAULT_GUEST_INDICES = [685, 686, 687, 688, 689, 690]
#DEFAULT_GUEST_INDICES = [723, 724, 722, 721, 725, 726]
DEFAULT_GUEST_INDICES = [757, 762, 758, 759, 760, 761]

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
        latticevalues = ReaderAndLoader.getLatticeValues(directory)
        numOfAtoms = ReaderAndLoader.getNumAtoms(directory)
        num_atoms_label.config(text=f"Number of Atoms: {numOfAtoms}")
        num_centroids_label.config(text=f"Number of Centroids: {numOfCentroids}")
        lattice_params_label.config(text=f"Lattice Parameters: {latticevalues}")
        
        # Create a scrollable frame for the tickboxes
        canvas = tk.Canvas(window)
        scrollbar = tk.Scrollbar(window, orient="vertical", command=canvas.yview)
        scrollable_frame = tk.Frame(canvas)

        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        # Create tick boxes for each centroid
        def check_all_centroids():
            for var in centroid_vars:
                var.set(True)

        def uncheck_all_centroids():
            for var in centroid_vars:
                var.set(False)

        check_all_button = tk.Button(scrollable_frame, text="Check All", command=check_all_centroids)
        check_all_button.pack(side="top")

        uncheck_all_button = tk.Button(scrollable_frame, text="Uncheck All", command=uncheck_all_centroids)
        uncheck_all_button.pack(side="top")

        global centroid_vars
        centroid_vars = []
        for i in range(numOfCentroids):
            var = tk.BooleanVar()
            centroid_vars.append(var)
            centroid_checkbox = tk.Checkbutton(scrollable_frame, text=f"Centroid {i+1}", variable=var)
            centroid_checkbox.pack(side="top")

                    
    else:
        saved_distance_data_label.config(text="No saved distances found.", fg="black")
        plotfromdistance_button.config(state=tk.DISABLED)  # Disable the plot button

def check_saved_data():
    file_path = os.path.join(directory_entry.get(), "sourcedata.npy")
    if os.path.isfile(file_path):
        #load_saved_data_button.config(state=tk.NORMAL)
        saved_data_label.config(text="Source data found in chosen directory!", fg="green")

    else:
        #load_saved_data_button.config(state=tk.DISABLED)
        saved_data_label.config(text="No source data found in chosen directory.", fg="black")
        
# Function to ask for directory
def ask_directory():
    directory = filedialog.askdirectory()
    if directory:
        calculate_distances(directory)
        
def plot_from_distances_file():
    file_path_distance = distances_entry.get()  
    plotFromSavedDistances(centroid_vars, file_path_distance)

# Function to calculate distances and save data
def calculate_distances(save_directory):
    
    steps = int(stepsToCalc_entry.get())
    index_filename = index_entry.get()
    directory = directory_entry.get()
    guest_indices_list = guest_entry.get().split()
    guest_indices = [int(i) for i in guest_indices_list]
    guest_identifier = guest_identifier_entry.get()
    file_path = os.path.join(directory_entry.get(), "sourcedata.npy")
    
    data, numAtoms, centroid_indices_flat, numOfCentroids, lattice_values, centroid_indices, numRelevantAtoms = ReaderAndLoader.load_data(directory, index_filename, file_path)  # Call the function to load the saved data
    distancePlotter = DistancePlotter(merged_df_values=data, numAtoms=numAtoms, centroid_indices_j=centroid_indices, numOfCentroids=numOfCentroids, lattice_values=lattice_values, numRelevantAtoms=numRelevantAtoms)
    file_path_distance_guest = distancePlotter.calcAndSaveAllDistancesForNSteps(steps, save_directory, guest_indices, guest_identifier)
    saved_data_label.config(text=f"Distances saved to: {file_path_distance_guest}")

# Create the main window
window = tk.Tk()
window.title("Trajectory Analyzer for centroid distances")

# Create a frame for the directory, prefix, and index filename widgets
input_frame = tk.Frame(window, highlightbackground="red", highlightthickness=4)
input_frame.pack(side="top", padx=20)

# Create a frame for the load and load saved data buttons
button_frame = tk.Frame(window, highlightbackground="blue", highlightthickness=4)
button_frame.pack( side = "top",padx=20)

# Create a new frame for the stepstocalc entry and calculate distances button
calc_frame = tk.Frame(window, highlightbackground="green", highlightthickness=4)
calc_frame.pack( side = "top",padx=20)

# Create a frame for the load precalculated distance file and plot chosen centroids section
load_distances_frame = tk.Frame(window, highlightbackground="orange", highlightthickness=4)
load_distances_frame.pack( side = "top",padx=20)

# Create a frame for the properties of simulation
properties_frame = tk.Frame(window, highlightbackground="purple", highlightthickness=4)
properties_frame.pack( side = "top", padx=20)

# Create a guide label with bold text
guide_label = tk.Label(window, text="""Guide: \n1. Make sure your chosen directory contains xyz files with the given prefix and some sort of integer for sorting. Example: MIL68Ga_1.xyz, MIL68Ga_2.xyz, etc.
2. Ensure that the index file contains n rows of m integers, where n is the number of centroids and m is the number of atoms in the centroid. Sorting not necessary.
3. Enter the indices of the guest atoms in the simulation, separated by spaces.
4. Supply an integer to identify the guest. This integer will identify the guest in the filename of the saved distances.
5. Click the "Read and save XYZ Files" button to save the data to a file. Source data must be present in the chosen directory for calculations to work.
6. Enter the number of steps to calculate and click the "Calculate Distances for Steps" button to calculate and save the distances for the chosen guest and number of steps. The location of the saved file will be displayed.
7. Click the "Browse file containing distances" button to load the distances file. You can then choose which centroids to plot. This will also show the number of atoms, centroids and lattice parameters.
8. ???
9: Profit """, fg="black", font=("Arial", 12, "bold"), justify="left")
guide_label.pack( padx=20)


# Create a label for the frame with bold text
label = tk.Label(input_frame, text="MANDATORY: File locations and guest parameters", fg="red", font=("Arial", 12, "bold"))
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

# Create the indices of guest label and entry
guest_label = tk.Label(input_frame, text="Indices of Guest:")
guest_label.pack()
guest_entry = tk.Entry(input_frame)
guest_entry.insert(tk.END, DEFAULT_GUEST_INDICES)  # Set default index filename
guest_entry.pack()
guest_entry.config(highlightbackground="red", highlightthickness=1)  # Add red rectangle

# Create the guest identifier label and entry
guest_identifier_label = tk.Label(input_frame, text="Guest Identifier:")
guest_identifier_label.pack()
guest_identifier_entry = tk.Entry(input_frame)
guest_identifier_entry.insert(tk.END, "1")  # Set default index filename
guest_identifier_entry.pack()
guest_identifier_entry.config(highlightbackground="red", highlightthickness=1)  # Add red rectangle


# Create a label for the frame
button_label = tk.Label(button_frame, text="Read and save trajectories to sourcefile.", fg="blue")
button_label.pack()

# Create the load button
load_button = tk.Button(button_frame, text="Read and save XYZ Files", command=read_and_save_data)
load_button.pack()

# Insert a red warning
warning_label = tk.Label(button_frame, text="WARNING: This will overwrite any existing source data in the chosen directory!", fg="red")
warning_label.pack()

# Create the saved data label
saved_data_label = tk.Label(button_frame, text="", fg="black")
saved_data_label.pack()

# Create the saved data label
saved_data_label = tk.Label(button_frame, text="", fg="black")
saved_data_label.pack()

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

lattice_params_label = tk.Label(properties_frame, text="Lattice Parameters: ")
lattice_params_label.pack()

# Check if the default directory already contains the "sourcedata" file
check_saved_data()

# Create the main loop
window.mainloop()

window.mainloop()
