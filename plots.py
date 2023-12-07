#Defining function to plot centroid distances over the course of the simulation

import matplotlib.pyplot as plt
import pandas as pd
import os
import math
import warnings
warnings.filterwarnings("ignore")

# Directory containing the .xyz files to analyze 
directory = '/home/mfi/Desktop/mfi/MIL-68Ga-guest'

# List of dataframes
dfs = []

# Save lattice parameters in seperate dataframe
latticefilename = 'MIL68Ga-guest-02.xyz'
lattice = pd.read_csv(os.path.join(directory,latticefilename), delimiter='\s+', header=None, skiprows=range(2,400000), nrows=1)
# Read the number of atoms from the lattice parameters
numAtoms = lattice[0][0]-1

# Save indices of carbons of which to form centroids
indexfilename = "indicessorted.dat"
centroid_indices = pd.read_csv(os.path.join(directory, indexfilename), delimiter = " ", header = None)
#Get number of centroids from the number of rows in the dataframe
numOfCentroides  = centroid_indices.index.size

#Loop through xyz files, create dataframes from them, clean them up, add incides and add them to the list of dataframes
for filename in sorted(os.listdir(directory)):
    if filename.endswith('.xyz') and filename.startswith('MIL68Ga-guest'):
        file_path = os.path.join(directory, filename)
        df = pd.read_csv(file_path, delimiter='\s+', header=None, skiprows=1)
        df['run'] = filename
        #drop all rows that are NaN
        df = df[~df.iloc[:,1].isna()]
        #drop all rows that are an X element
        df = df[df[0] != 'X']
        #reset the Index of the dataframe
        df.reset_index(drop=True, inplace=True)
        #add cleanedIndex to each row representing its index in the simulation, going from 1 to numAtoms then start at 1 again
        df['cleanedIndex'] = (df.index)%(numAtoms)+1
        #append to the list of dataframes
        dfs.append(df)
        

#Concat dataframes into one big dataframe, containing all the coordinates of all the simulations. Keep cleanedIndex to identify atoms
merged_df = pd.concat(dfs, axis=0, ignore_index=True)

#Add new index assigning each row to its respective step in the simulation. First 720 should be 1, next 720 rows should be 2, etc.
merged_df['runningIndex'] = ((merged_df.index)//(numAtoms)).astype(int)+1

#Print the dataframe to check if everything worked
print(merged_df)

#Find max value of runningIndex column, i.e. number of steps performed in total
numSteps = merged_df['runningIndex'].max()

#Function that, when given a centroid index and the dataframe of a specific step (i.e. coordinate of all atoms in this step), calculates the coordinates of centroid specified by the index.
#column_of_interest specifies which column of the dataframe to use for the calculation. Usually column 1, 2 or 3, corresponding to x, y or z coordinates
def calculateCentroidCoordinates(centroidNumber, stepframe, column_of_interest):
    sum = 0
    for index in centroid_indices.loc[centroidNumber,:]:
        sum = sum + ((stepframe.loc[stepframe['cleanedIndex'] == index, column_of_interest]).iloc[0]**2)
        continue
    return math.sqrt(sum)


#Loop over every step of the simulation. i is the step number given by the runningIndex column, j is the centroid number that is looped over in each step 
# and calculate coordinates of centroid j in step i. Add these coordinates to the merged_df dataframe
for i in range (1,2):
    stepframe = merged_df[merged_df['runningIndex']==i]
    for j in range (1, numOfCentroides+1):
        print(j)
        merged_df.loc[merged_df['runningIndex']==i, 'x_coords_of centroid ' + str(j)] = calculateCentroidCoordinates(j-1, stepframe, 1)
        merged_df.loc[merged_df['runningIndex']==i, 'y_coords_of centroid ' + str(j)] = calculateCentroidCoordinates(j-1, stepframe, 2)
        merged_df.loc[merged_df['runningIndex']==i, 'z_coords_of centroid ' + str(j)] = calculateCentroidCoordinates(j-1, stepframe, 3)
        continue
    continue

#Define centroid of the guest molecule by the indices of it's carbon ring
guest_indices = [686, 687, 688, 689, 690, 685]
 
#Function that, when given a dataframe of a specific step (i.e. coordinate of all atoms in this step), calculates the coordinates of the guest molecule
def calculateGuestCoordinates(stepframe):
    sum_x = 0
    sum_y = 0
    sum_z = 0
    for index in guest_indices:
        sum_x = sum_x + ((stepframe.loc[stepframe['cleanedIndex'] == index, 1]).iloc[0]**2)
        sum_y = sum_y + ((stepframe.loc[stepframe['cleanedIndex'] == index, 2]).iloc[0]**2)
        sum_z = sum_z + ((stepframe.loc[stepframe['cleanedIndex'] == index, 3]).iloc[0]**2)
        continue
    return math.sqrt(sum_x), math.sqrt(sum_y), math.sqrt(sum_z)

# Loop over every step of the simulation
for i in range(1, 2):
    # Get the stepframe for the current step
    stepframe = merged_df[merged_df['runningIndex'] == i]
        
    # Calculate the Guest coordinates for the current step
    x_coords, y_coords, z_coords = calculateGuestCoordinates(stepframe)
        
    # Write the Guest coordinates into new columns
    merged_df.loc[merged_df['runningIndex'] == i, 'Guest_x_coords'] = x_coords
    merged_df.loc[merged_df['runningIndex'] == i, 'Guest_y_coords'] = y_coords
    merged_df.loc[merged_df['runningIndex'] == i, 'Guest_z_coords'] = z_coords
    
def calculateDistanceToCentroid(stepframe, centroidNumber):
    x_distance = stepframe['Guest_x_coords'] - stepframe['x_coords_of centroid ' + str(centroidNumber)]
    y_distance = stepframe['Guest_y_coords'] - stepframe['y_coords_of centroid ' + str(centroidNumber)]
    z_distance = stepframe['Guest_z_coords'] - stepframe['z_coords_of centroid ' + str(centroidNumber)]
    distance = ((x_distance**2 + y_distance**2 + z_distance**2)**0.5).astype(float)
    return distance

#Finally, in each step, calculate the distance between each centroid and the guest molecule
for i in range(1, 2):
    # Get the stepframe for the current step
    stepframe = merged_df[merged_df['runningIndex'] == i]
    for j in range (1, numOfCentroides+1):
        merged_df.loc[merged_df['runningIndex']==i, 'distance_to_centroid ' + str(j)] = calculateDistanceToCentroid(stepframe, j)
        #stepframe['x_coords_of centroid ' + str(j)] = calculateCentroidCoordinates(j-1, stepframe, 1)
        #stepframe['y_coords_of centroid ' + str(j)] = calculateCentroidCoordinates(j-1, stepframe, 2)
        #stepframe['z_coords_of centroid ' + str(j)] = calculateCentroidCoordinates(j-1, stepframe, 3)
        #print(stepframe)
        continue
    
    