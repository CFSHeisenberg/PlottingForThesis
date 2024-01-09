import time
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
import warnings
import numpy as np
from concurrent import futures
import sys
import PQAnalysis as pq
from PQAnalysis.io.trajectoryReader import TrajectoryReader


warnings.filterwarnings("ignore")

# Directory containing the .xyz files to analyze 
directory = '/home/mfi/Desktop/mfi/MIL-68Ga-guest'

# List of dataframes
dfs = []

# Save lattice parameters in separate dataframe
latticefilename = 'MIL68Ga-guest-02.xyz'
lattice = pd.read_csv(os.path.join(directory, latticefilename), delimiter='\s+', header=None, skiprows=range(2, 400000), nrows=1)
# Read the number of atoms from the lattice parameters
numAtoms = lattice[0][0] - 1

# Save indices of carbons of which to form centroids
indexfilename = "indicessorted.dat"
centroid_indices = pd.read_csv(os.path.join(directory, indexfilename), delimiter=" ", header=None)

# Get number of centroids from the number of rows in the dataframe
numOfCentroids = centroid_indices.index.size

# Loop through xyz files, create dataframes from them, clean them up, add indices and add them to the list of dataframes
for i, filename in enumerate(sorted(os.listdir(directory))):
    if filename.endswith('1.xyz') and filename.startswith('MIL68Ga-guest'):
        file_path = os.path.join(directory, filename)
        df = pd.read_csv(file_path, delimiter='\s+', header=None, skiprows=1)
        df['run'] = filename
        # drop all rows that are NaN
        df = df[~df.iloc[:, 1].isna()]
        # drop all rows that are an X element
        df = df[df[0] != 'X']
        # reset the Index of the dataframe
        df.reset_index(drop=True, inplace=True)
        # add cleanedIndex to each row representing its index in the simulation, going from 1 to numAtoms then start at 1 again
        df['cleanedIndex'] = (df.index) % (numAtoms) + 1
        # append to the list of dataframes
        dfs.append(df)

# Merge dataframes into one big dataframe, containing all the coordinates of all the simulations. Keep cleanedIndex to identify atoms
merged_df = pd.concat(dfs, axis=0, ignore_index=True)

# Add new index assigning each row to its respective step in the simulation. First 720 should be 1, next 720 rows should be 2, etc.
merged_df['runningIndex'] = ((merged_df.index) // (numAtoms)).astype(int) + 1

# Find max value of runningIndex column, i.e. number of steps performed in total
numSteps = merged_df['runningIndex'].max()

# Start timer
start_time = time.time()

# Hard coded indices of guest atom carbon ring
guest_indices = [685, 686, 687, 688, 689, 690]

# Get the indices of the atoms making up centroids
centroid_indices_j = centroid_indices.values.flatten()

# Convert data to NumPy arrays outside the loop
merged_df_values = merged_df.values
#print(merged_df_values.shape)

# Number of steps to plot
stepstoplot = 150

# Get lattice values
lattice_values = np.array(lattice.values[0][1:4])

    
def imaging(coords, previous_coords, i):
    difference = np.array(coords - previous_coords)
    difference_new = difference - np.round(np.divide(difference, lattice_values).astype(float)) * lattice_values 
    coords_changed = previous_coords + difference_new
    return coords_changed

# Function to calculate distances from the guest centroid to each host centroid for a given step
def calculate_distance(i):
    if i == 1:  
        
        #Get partial numpy array of current step. Takes only the rows of the current step
        stepframe_values = merged_df_values[(i-1)*numAtoms:i*numAtoms, :]

        # Get partial numpy array of shape (108,7) only containing those entries of stepframe_values whose fifth element, i.e. the cleanedIndex are in centroid_indices_j
        stepframe_j_values = stepframe_values[np.isin(stepframe_values[:, 5], centroid_indices_j)]

        # Get coordinates of the atoms making up centroids in current step
        centroid_coords = stepframe_j_values[:, 1:4]

        # Get partial numpy array of shape (6, 7) only containing those entries of stepframe_values whose fifth element, i.e. the cleanedIndex are in guest_indices
        stepframe_g_values = stepframe_values[np.isin(stepframe_values[:, 5], guest_indices)]
        
        # Get coordinates of the relevant guest atoms in current step
        guest_coords = stepframe_g_values[:, 1:4]

        # Get coordinates of cartesian centroids for both host and guest via the mean of their coordinates
        center_guest_coords = np.mean(guest_coords, axis=0)
        center_host_coords = np.mean(np.split(centroid_coords, numOfCentroids), axis=1)


        # Calculate distance between each centroid to the guest centroid
        distance_to_center = np.linalg.norm((center_guest_coords - center_host_coords).astype(float), axis=1)
        
    else:
    
        #Get partial numpy array of current step. Takes only the rows of the current step
        stepframe_values = merged_df_values[(i-1)*numAtoms:i*numAtoms, :]
        
        #Get partial numpy array of previous step. Takes only the rows of the previous step
        stepframe_values_previous = merged_df_values[(i-2)*numAtoms:(i-1)*numAtoms, :]
        
        # Get partial numpy array of shape (108,7) of the current step only containing those entries of stepframe_values whose fifth element, i.e. the cleanedIndex are in centroid_indices_j
        stepframe_j_values = stepframe_values[np.isin(stepframe_values[:, 5], centroid_indices_j)]
        
        # Get partial numpy array of shape (108,7) of the previous step only containing those entries of stepframe_values_next whose fifth element, i.e. the cleanedIndex are in centroid_indices_j
        stepframe_j_values_previous = stepframe_values_previous[np.isin(stepframe_values_previous[:, 5], centroid_indices_j)]       

        # Get coordinates of the atoms making up centroids in current step
        centroid_coords = stepframe_j_values[:, 1:4]
        # Get coordinates of the atoms making up centroids in previous step
        centroid_coords_previous = stepframe_j_values_previous[:, 1:4]

        # Get partial numpy array of shape (6, 7) only containing those entries of stepframe_values whose fifth element, i.e. the cleanedIndex are in guest_indices
        stepframe_g_values = stepframe_values[np.isin(stepframe_values[:, 5], guest_indices)]
        
        # Get coordinates of the relevant guest atoms in current step
        guest_coords = stepframe_g_values[:, 1:4]
        
        #Compare centroid coordinates in next step to centroid coordinates in current step, image them if necessary
        centroid_coords_imaged = imaging(centroid_coords, centroid_coords_previous, i)
        
        #Overwrite the values in the dataframe accessed by stepframe_j_values with centroid_coords_imaged for the next step
        merged_df_values[(i-1)*numAtoms:i*numAtoms, 1:4][np.isin(stepframe_values[:,5], centroid_indices_j)] = centroid_coords_imaged

        # Get coordinates of cartesian centroids for both host and guest centroids via the mean of their coordinates
        center_guest_coords = np.mean(guest_coords, axis=0)
        center_host_coords_imaged = np.mean(np.split(centroid_coords_imaged, numOfCentroids), axis=1)
        

        # Calculate distance between each centroid to the guest centroid
        distance_to_center = np.linalg.norm((center_guest_coords - center_host_coords_imaged).astype(float), axis=1)

        # Return the calculated distances
    return distance_to_center

# Use ThreadPoolExecutor to parallelize the distance calculation
with futures.ThreadPoolExecutor(1) as executor:
    
    # List to store the future objectst
    futures_list = []

    # Loop over every step of the simulation
    for i in range(1, stepstoplot):
        # Print the step number
        print(i)
        # Submit the distance calculation task to the executor
        future = executor.submit(calculate_distance, i)
        futures_list.append(future)

    # Wait for all the tasks to complete and retrieve the results
    results = [future.result() for future in futures_list]

print("DONE")

# Plot the distances for each centroid to the guest centroid
for j in range(1, numOfCentroids + 1):
    plt.plot(range(1, stepstoplot), [result[j-1] for result in results], label=f'Centroid {j}')

# End timer and print total calculation time
end_time = time.time()
total_time = end_time - start_time
print(f"Total calculation time: {total_time} seconds")

# Edit and show plot
plt.xlabel('Step Number')
plt.ylabel('Distance')
plt.legend()
plt.show()
