import os
import numpy as np
import pandas as pd


class ReaderAndLoader:
    def __init__(self):
        self.numOfCentroids = None
        self.numAtoms = None
        self.numSteps = None
        self.centroid_indices = None
        self.centroid_indices_flat = None
        self.lattice_values = None
        self.guest_indices = [425, 426, 427, 428, 429, 430]

    def read_and_save_data(self, directory, file_prefix, index_filename):
        
        print("Reading and saving data...")

        # List of dataframes
        dfs = []
        dfs_lattices = []

        # Save indices of carbons of which to form centroids
        self.centroid_indices = pd.read_csv(os.path.join(directory, index_filename), header=None)
        
        # Drop columns in centroid_indices that are NaN
        self.centroid_indices = self.centroid_indices.dropna(axis=1)

        print("indices when saving" , self.centroid_indices)

        # Get number of centroids from the number of rows in the dataframe
        self.numOfCentroids = self.centroid_indices.index.size
        
        print(self.centroid_indices)
        
        # Loop through xyz files, create dataframes from them, clean them up, add indices and add them to the list of dataframes
        for i, filename in enumerate(sorted(os.listdir(directory))):
            print("filename before loop" ,filename)
            if filename.endswith('.xyz') and filename.startswith(file_prefix):
                file_path = os.path.join(directory, filename)
                df_lattices = pd.read_csv(file_path, header=None, on_bad_lines='skip')
                df_lattices = df_lattices[0].str.split('\s+', expand=True)
                df_lattices['run'] = filename
                # drop all rows that are NaN
                df_lattices = df_lattices.dropna()
                # drop all rows that are an X element
                # reset the Index of the dataframe
                df_lattices.reset_index(drop=True, inplace=True)
                # add cleanedIndex to each row representing its index in the simulation, going from 1 to numAtoms then start at 1 again
                df_lattices['runningIndex'] = (df_lattices.index) + 1
                # append to the list of dataframes
                dfs_lattices.append(df_lattices)
                
        # Read the number of atoms from the lattice parameters
        self.numAtoms = int(dfs_lattices[0][0][0]) - 1

        # Loop through xyz files, create dataframes from them, clean them up, add indices and add them to the list of dataframes
        for i, filename in enumerate(sorted(os.listdir(directory))):
            print("filename before loop" ,filename)
            if filename.endswith('.xyz') and filename.startswith(file_prefix):
                file_path = os.path.join(directory, filename)
                df = pd.read_csv(file_path, delimiter='\s+', header=None, skiprows=1, on_bad_lines='skip')
                df['run'] = filename
                # drop all rows that are NaN
                df = df[~df.iloc[:, 1].isna()]
                # drop all rows that are an X element
                df = df[df[0] != 'X']
                # reset the Index of the dataframe
                df.reset_index(drop=True, inplace=True)
                # add cleanedIndex to each row representing its index in the simulation, going from 1 to numAtoms then start at 1 again
                df['cleanedIndex'] = (df.index) % (self.numAtoms) + 1
                # append to the list of dataframes
                dfs.append(df)

        # Concatenate the dataframes into a single dataframe
        merged_df = pd.concat(dfs)
        merged_df_lattices = pd.concat(dfs_lattices)

        # Add new index assigning each row to its respective step in the simulation. First 720 should be 1, next 720 rows should be 2, etc.
        merged_df['runningIndex'] = ((merged_df.index) // (self.numAtoms)).astype(int) + 1

        # Find max value of runningIndex column, i.e. number of steps performed in total
        self.numSteps = merged_df['runningIndex'].max()

        # Get the indices of the atoms making up centroids
        self.centroid_indices_flat = self.centroid_indices.values.flatten()

        # Convert data to NumPy arrays outside the loop
        merged_df_values = merged_df.values
        
        # Convert lattice data to NumPy arrays outside the loop
        merged_df_lattices_values = merged_df_lattices.values

        # Save the merged_df_values numpy array to a file
        saved_path = os.path.join(directory, 'sourcedata.npy')
        saved_path_lattices = os.path.join(directory, 'sourcedata_lattices.npy')
        np.save(saved_path, merged_df_values)
        np.save(saved_path_lattices, merged_df_lattices_values)

        return saved_path, saved_path_lattices


    def getNumCentroids (self, directory, index_filename):
        self.numOfCentroids = pd.read_csv(os.path.join(directory, index_filename), delimiter=" ", header=None).index.size
        return self.numOfCentroids

    def load_data(self, directory, index_filename, saved_path, saved_path_lattices):
        # Load the data from the saved file
        data = np.load(saved_path, allow_pickle=True)
        data_lat = np.load(saved_path_lattices, allow_pickle=True)
        
        
        # Save indices of carbons of which to form centroids
        self.centroid_indices = pd.read_csv(os.path.join(directory, index_filename), delimiter=" ", header=None)
        self.centroid_indices = self.centroid_indices.dropna(axis=1)
        
        self.centroid_indices_flat = self.centroid_indices.values.flatten()

        # Get number of centroids from the number of rows in the dataframe
        self.numOfCentroids = self.centroid_indices.index.size
        
        #print the highest element in the sixth column
        
        self.numAtoms = np.max(data[:, 5])
        
        self.numSteps = (np.max(data[:, 6]))
        return data, self.numAtoms, self.centroid_indices_flat, self.guest_indices, self.numOfCentroids, self.lattice_values, self.centroid_indices, data_lat
    
    