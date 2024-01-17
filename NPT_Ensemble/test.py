import os
import numpy as np
import pandas as pd

directory = "/home/mfi/Desktop/mfi/filesManuelOtt"

# List of dataframes
dfs = []
dfs_lattices = []

# Save lattice parameters in separate dataframe
latticefilename = 'MOF5-md-04.xyz'
lattice = pd.read_csv(os.path.join(directory, latticefilename), delimiter='\s+', header=None, skiprows=range(2, 400000), nrows=1)

# Read the number of atoms from the lattice parameters
numAtoms = lattice[0][0] - 1

# Save indices of carbons of which to form centroids
centroid_indices = pd.read_csv(os.path.join(directory, 'indicessorted.dat'), delimiter=" ", header=None)

#print(centroid_indices.dropna(axis=1))

# Get number of centroids from the number of rows in the dataframe
numOfCentroids = centroid_indices.index.size

for i, filename in enumerate(sorted(os.listdir(directory))):
    #print("filename before loop", filename)
    if filename.endswith('.xyz') and filename.startswith("MOF5"):
        file_path = os.path.join(directory, filename)
        df = pd.read_csv(file_path, delimiter='\s+', header=None, skiprows=1, on_bad_lines='skip')
        df_lattices = pd.read_csv(file_path, header=None, on_bad_lines='skip')
        df_lattices = df_lattices[0].str.split('\s+', expand=True)
        df['run'] = filename
        df_lattices['run'] = filename
        # drop all rows that are NaN
        df = df.dropna()
        df_lattices = df_lattices.dropna()
        # drop all rows that are an X element
        df = df[df[0] != 'X']
        # reset the Index of the dataframe
        df.reset_index(drop=True, inplace=True)
        df_lattices.reset_index(drop=True, inplace=True)
        # add cleanedIndex to each row representing its index in the simulation, going from 1 to numAtoms then start at 1 again
        df['cleanedIndex'] = (df.index) % (numAtoms) + 1
        df_lattices['runningIndex'] = (df_lattices.index) + 1
        # append to the list of dataframes
        dfs.append(df)
        dfs_lattices.append(df_lattices)
                

# Concatenate the dataframes into a single dataframe
merged_df = pd.concat(dfs)

merged_df_lattices = pd.concat(dfs_lattices)

# Add new index assigning each row to its respective step in the simulation. First 720 should be 1, next 720 rows should be 2, etc.
merged_df['runningIndex'] = ((merged_df.index) // (numAtoms)).astype(int) + 1

# Find max value of runningIndex column, i.e. number of steps performed in total
numSteps = merged_df['runningIndex'].max()

# Get the indices of the atoms making up centroids
centroid_indices_flat = centroid_indices.values.flatten()

# Convert data to NumPy arrays outside the loop
merged_df_values = merged_df.values

# Convert data to NumPy arrays outside the loop
merged_df_lattices_values = merged_df_lattices.values

#print(merged_df_lattices_values)
print(merged_df_lattices_values[1,1:4])

# Get lattice values
lattice_values = np.array(lattice.values[0][1:4])

# Save the merged_df_values numpy array to a file
saved_path = os.path.join(directory, 'sourcedata.npy')
saved_path_lattices = os.path.join(directory, 'sourcedata_lattices.npy')
np.save(saved_path, merged_df_values)
np.save(saved_path_lattices, merged_df_lattices_values)