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
        self.latticefilename = 'MIL68Ga-3rdguest-06.xyz'
        self.numRelevantAtoms = None

    def read_and_save_data(self, directory, file_prefix, index_filename):
        
        print("Reading and saving data...")

        # List of dataframes
        dfs = []

        # Save lattice parameters in separate dataframe
        lattice = pd.read_csv(os.path.join(directory, self.latticefilename), delimiter='\s+', header=None, skiprows=range(2, 400000), nrows=1)
        
        # Read the number of atoms from the lattice parameters
        self.numAtoms = lattice[0][0] - 1

        # Save indices of carbons of which to form centroids
        self.centroid_indices = pd.read_csv(os.path.join(directory, index_filename), delimiter=" ", header=None)

        # Get number of centroids from the number of rows in the dataframe
        self.numOfCentroids = self.centroid_indices.index.size
        

    

        # Loop through xyz files, create dataframes from them, clean them up, add indices and add them to the list of dataframes
        for i, filename in enumerate(sorted(os.listdir(directory))):
            if filename.endswith('.xyz') and filename.startswith(file_prefix):
                file_path = os.path.join(directory, filename)
                print("file_path: ", file_path)
                df = pd.read_csv(file_path, delimiter='\s+', header=None, skiprows=1)
                # drop all rows that are NaN
                df = df[~df.iloc[:, 1].isna()]
                # drop all rows that are an X element
                df = df[df[0] != 'X']
                # reset the Index of the dataframe
                df.reset_index(drop=True, inplace=True)
                # add cleanedIndex to each row representing its index in the simulation, going from 1 to numAtoms then start at 1 again
                df['cleanedIndex'] = (df.index) % (self.numAtoms) + 1
                # If the first element is Ga, H or O, replace the first three entries with 0
                # drop all rows that are an Ga or H or O element
                #df = df[df[0] != 'Ga']
                df = df[df[0] != 'H']
                df = df[df[0] != 'O']
                # drop the first column
                df = df.drop(columns=[0])
                # append to the list of dataframes
                print(df)
                dfs.append(df)

        # Concatenate the dataframes into a single dataframe
        merged_df = pd.concat(dfs)

        # Reset the index of the merged dataframe
        merged_df.reset_index(drop=True, inplace=True)
        
        # Add new index assigning each row to its respective step in the simulation. First 720 should be 1, next 720 rows should be 2, etc.
        #merged_df['runningIndex'] = ((merged_df.index) // (self.numAtoms)).astype(int) + 1
        #merged_df = merged_df[merged_df[0] != 'Ga']
        #merged_df = merged_df[merged_df[0] != 'O']
        #merged_df = merged_df[merged_df[0] != 'H']
        #merged_df  = merged_df.drop(columns=[0])




        # Find max value of runningIndex column, i.e. number of steps performed in total
        #self.numSteps = merged_df['runningIndex'].max()

        # Get the indices of the atoms making up centroids
        self.centroid_indices_flat = self.centroid_indices.values.flatten()

        # Convert data to NumPy arrays outside the loop
        merged_df_values = merged_df.values
        print(merged_df_values)

        # Get lattice values
        self.lattice_values = np.array(lattice.values[0][1:4])

        # Save the merged_df_values numpy array to a file
        saved_path = os.path.join(directory, 'sourcedata.npy')
        np.save(saved_path, merged_df_values)

        return saved_path
    
    def getNumCentroids (self, directory, index_filename):
        self.numOfCentroids = pd.read_csv(os.path.join(directory, index_filename), delimiter=" ", header=None).index.size
        return self.numOfCentroids
    
    def getLatticeValues (self, directory):
        lattice = pd.read_csv(os.path.join(directory, self.latticefilename), delimiter='\s+', header=None, skiprows=range(2, 400000), nrows=1)
        self.lattice_values = np.array(lattice.values[0][1:4])
        return self.lattice_values
    
    def getNumAtoms (self, directory):
        lattice = pd.read_csv(os.path.join(directory, self.latticefilename), delimiter='\s+', header=None, skiprows=range(2, 400000), nrows=1)
        self.numAtoms = lattice[0][0] - 1
        return self.numAtoms

    def load_data(self, directory, index_filename, saved_path):
        # Load the data from the saved file
        data = np.load(saved_path, allow_pickle=True)
        
        # Read the first numAtoms lines of the lattice file to get the number of relevant atoms. Count only C and N atoms
        self.numAtoms = self.getNumAtoms(directory)
        temp = pd.read_csv(os.path.join(directory, self.latticefilename), delimiter='\s+', header=None, skiprows=1, nrows=self.numAtoms)
        # Count the number of C and N atoms
        self.numRelevantAtoms = temp[temp[0].isin(['C', 'N', 'Ga'])].index.size
        print("relevant atoms: ", self.numRelevantAtoms)
        #print("Number of relevant atoms: ", self.numRelevantAtoms)
        # Save lattice parameters in separate dataframe
        lattice = pd.read_csv(os.path.join(directory, self.latticefilename), delimiter='\s+', header=None, skiprows=range(2, 400000), nrows=1)
        # Get lattice values
        self.lattice_values = np.array(lattice.values[0][1:4])
        
        # Save indices of carbons of which to form centroids
        self.centroid_indices = pd.read_csv(os.path.join(directory, index_filename), delimiter=" ", header=None)

        # Get number of centroids from the number of rows in the dataframe
        self.numOfCentroids = self.centroid_indices.index.size
        
        #print the highest element in the sixth column
        
        self.numAtoms = self.getNumAtoms(directory)
        
        #self.numSteps = (np.max(data[:, 4]))
        
        return data, self.numAtoms, self.centroid_indices_flat, self.numOfCentroids, self.lattice_values, self.centroid_indices, self.numRelevantAtoms
    
    