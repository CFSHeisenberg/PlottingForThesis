from concurrent import futures
import os
from tkinter import Entry
from matplotlib import pyplot as plt
import numpy as np
from tqdm import tqdm
import matplotlib.ticker as ticker


class DistancePlotter:
    def __init__(self, merged_df_values, numAtoms, centroid_indices_j, numOfCentroids, lattice_values, numRelevantAtoms):
        self.merged_df_values = merged_df_values
        self.numAtoms = numAtoms
        self.centroid_indices_j = centroid_indices_j
        self.numOfCentroids = numOfCentroids
        self.lattice_values = lattice_values
        self.numRelevantAtoms = numRelevantAtoms

    def imaging(self, coords, previous_coords, i):
        difference = np.array(coords - previous_coords)
        difference_new = difference - np.round(np.divide(difference, self.lattice_values).astype(float)) * self.lattice_values 
        coords_changed = previous_coords + difference_new
        return coords_changed
    
    def imageDistance(self, guest_coords_imaged, host_coords_imaged):
        # Calculate cartesian distance between guest and host centroids
        distance = host_coords_imaged - guest_coords_imaged
        # Image the distance
        distance_to_image = distance - np.round(np.divide(distance, self.lattice_values).astype(float)) * self.lattice_values
        # Add the imaged distance to the guest coordinates
        host_coords_imaged_temp = guest_coords_imaged + distance_to_image
        # Get cartesian distance between guest and imaged host centroids
        distance_to_center = np.linalg.norm((distance_to_image).astype(float), axis=1)
        return distance_to_center
        

    def imageGuestCentroid(self, guest_centroid_coords):
        # Save the first centroid coordinates as a starting point
        centroid_coords_origin = guest_centroid_coords[0]
        
        # Get the difference between the other centroid coordinates and the first centroid coordinates
        centroid_coords_difference = guest_centroid_coords - centroid_coords_origin
        
        # Image the centroid coordinates
        centroid_coords_imaged = centroid_coords_origin + centroid_coords_difference - np.round(np.divide(centroid_coords_difference, self.lattice_values).astype(float)) * self.lattice_values
        
        return centroid_coords_imaged
    
    def imageHostCentroids(self, host_centroid_coords):
        splitted_host_centroid_coords = np.split(host_centroid_coords, self.numOfCentroids)
        imaged_centroids = []
        # Save the first centroid coordinates as a starting point
        for centroid in splitted_host_centroid_coords:
            centroid_coords_origin = centroid[0]

            # Get the difference between the other centroid coordinates and the first centroid coordinates
            centroid_coords_difference = centroid - centroid_coords_origin
        
            # Image the centroid coordinates
            centroid_coords_imaged = centroid_coords_origin + centroid_coords_difference - np.round(np.divide(centroid_coords_difference, self.lattice_values).astype(float)) * self.lattice_values
            
            imaged_centroids.append(centroid_coords_imaged)
        
        return imaged_centroids
        

    def calculate_distance(self, i, selected_centroids):
        # Get partial numpy array of current step. Takes only the rows of the current step
        stepframe_values = self.merged_df_values[(i-1)*self.numRelevantAtoms:i*self.numRelevantAtoms, :]
        #print(i, stepframe_values)

        # Get partial numpy array of previous step. Takes only the rows of the previous step
        #stepframe_values_previous = self.merged_df_values[(i-2)*self.numAtoms:(i-1)*self.numAtoms, :]

        # Get partial numpy array of shape (108,7) of the current step only containing those entries of stepframe_values whose fifth element, i.e. the cleanedIndex are in centroid_indices_j
        stepframe_j_values = np.array([stepframe_values[stepframe_values[:, 3] == index] for index in self.centroid_indices_j.values.ravel()])
        #print("centroid_indies_j: ", self.centroid_indices_j)
        # Get partial numpy array of shape (108,7) of the previous step only containing those entries of stepframe_values_next whose fifth element, i.e. the cleanedIndex are in centroid_indices_j
        #stepframe_j_values_previous = stepframe_values_previous[np.isin(stepframe_values_previous[:, 5], self.centroid_indices_j)]       
        stepframe_j_values = stepframe_j_values.reshape(-1, stepframe_j_values.shape[-1])
        # Get coordinates of the atoms making up centroids in current step
        centroid_coords = stepframe_j_values[:, 0:3]
        #print("stepframe_j_values", stepframe_j_values)
        centroid_coords = self.imageHostCentroids(centroid_coords)
        # Get coordinates of the atoms making up centroids in previous step
        #centroid_coords_previous = stepframe_j_values_previous[:, 1:4]

        # Get partial numpy array of shape (6, 7) only containing those entries of stepframe_values whose fifth element, i.e. the cleanedIndex are in guest_indices
        stepframe_g_values = stepframe_values[np.isin(stepframe_values[:, 3], self.guest_indices)]

        # Get coordinates of the relevant guest atoms in current step
        guest_coords = stepframe_g_values[:, 0:3]
        guest_coords = self.imageGuestCentroid(guest_coords)

        # Compare centroid coordinates in next step to centroid coordinates in current step, image them if necessary
        #centroid_coords_imaged = self.imaging(centroid_coords, centroid_coords_previous, i)

        # Overwrite the values in the dataframe accessed by stepframe_j_values with centroid_coords_imaged for the next step
        #self.merged_df_values[(i-1)*self.numAtoms:i*self.numAtoms, 1:4][np.isin(stepframe_values[:,5], self.centroid_indices_j)] = centroid_coords_imaged

        # Get coordinates of cartesian centroids for both host and guest centroids via the mean of their coordinates
        center_guest_coords = np.mean(guest_coords, axis=0)
        center_host_coords = np.mean(centroid_coords, axis=1)

        # Calculate distance between each centroid to the guest centroid
        #distance_to_center = np.linalg.norm((center_guest_coords - center_host_coords).astype(float), axis=1)
        distance_to_center = self.imageDistance(center_guest_coords, center_host_coords)

        # Return the calculated distances for selected centroids
        return [distance_to_center[j-1] for j in selected_centroids]
    

    def calcAndSaveAllDistancesForNSteps(self, stepsToCalculate, directory, guest_indices_to_calc, guest_identifier):
        filename = 'distances_for_guest_' + str(guest_identifier) +'_and_'+ str(stepsToCalculate) + '_steps_NVT.npy'
        self.guest_indices = guest_indices_to_calc
        # Set the batch size for distance calculations
        batch_size = 1000

        # Use ThreadPoolExecutor to parallelize the distance calculation
        with futures.ThreadPoolExecutor(1) as executor:
            # List to store the future objects
            futures_list = []

            # List to store all the results
            all_results = []

            # Loop over the specified steps in batches
            for start_step in tqdm(range(1, stepsToCalculate + 1, batch_size), desc="Calculating distances"):
                end_step = min(start_step + batch_size - 1, stepsToCalculate)

                # Submit the distance calculation tasks to the executor
                for i in range(start_step, end_step + 1):
                    future = executor.submit(self.calculate_distance, i, np.arange(1, self.numOfCentroids + 1))
                    futures_list.append(future)

                # Wait for all the tasks to complete and retrieve the results
                results = [future.result() for future in futures_list]

                # Add the results to the list of all results
                all_results.extend(results)

                # Clear the futures list to free up memory
                futures_list.clear()

            # Save all the results to file
            savepath = os.path.join(directory, filename)
            np.save(savepath, all_results)

        return savepath
            

    def plot_distances(self, stepsToPlot, centroid_vars, guest_vars):
        # Get the total steps to plot from the entry field
        total_steps = stepsToPlot
        if guest_vars.get() == 0:
            self.guest_indices = self.guest_indices_1
        else:
            self.guest_indices = self.guest_indices_2
        try:
            # Convert the total steps to an integer
            total_steps = int(total_steps)

            # Get the selected centroids to plot
            selected_centroids = [var.get() for var in centroid_vars]
            selected_centroids = [i+1 for i, val in enumerate(selected_centroids) if val == 1]

            # Use ThreadPoolExecutor to parallelize the distance calculation
            with futures.ThreadPoolExecutor(1) as executor:
                # List to store the future objects
                futures_list = []

                # Loop over the specified steps
                for i in range(1, total_steps + 1):
                    # Submit the distance calculation task to the executor
                    future = executor.submit(self.calculate_distance, i, selected_centroids)
                    futures_list.append(future)

                # Wait for all the tasks to complete and retrieve the results
                results = [future.result() for future in futures_list]

            # Plot the distances for each selected centroid to the guest centroid
            for j, centroid in enumerate(selected_centroids):
                plt.plot(range(1, total_steps + 1), [result[j] for result in results], label=f'Centroid {centroid}')

            # Edit and show plot
            plt.xlabel('Step Number')
            plt.ylabel('Distance to Guest Centroid [Å]')
            plt.legend()
            plt.show()
        except ValueError:
            print("Invalid input. Please enter a valid integer for the total steps.")
            
def plotFromSavedDistances(centroid_vars, file_path):
        
    window_size = 1000    
        
    # Open numpy array from file
    results = np.load(file_path)
        
    # Get the steps to plot from the name of the file
    stepsToPlot = int(file_path.split('_')[5])
    
    # Get the selected guest
    guestNumber = int(file_path.split('_')[3])
        
    # Get the selected centroids to plot
    selected_centroids = [var.get() for var in centroid_vars]
    selected_centroids = [i+1 for i, val in enumerate(selected_centroids) if val == 1]
    print(selected_centroids)

    # Define a color palette for the centroids
    colors = plt.get_cmap('Set1').colors

    # Plot the distances for each selected centroid to the guest centroid
    for j, centroid in enumerate(selected_centroids):
        
        centroid_index = centroid - 1

        # Plot the actual data with a different color for each centroid
        plt.plot([i * 4 / 1_000_000 for i in range(1, stepsToPlot + 1)], [result[centroid_index] for result in results], color=colors[j % len(colors)], alpha=0.3)

        # Calculate the running average
        running_avg = np.convolve([result[centroid_index] for result in results], np.ones(window_size)/window_size, mode='valid')

        # Plot the running average with the same color as the source data
        plt.plot([i * 4 / 1_000_000 for i in range(window_size//2, stepsToPlot - window_size//2 + 1)], running_avg, color=colors[j % len(colors)])

    # Edit and show plot
    
    min_x = min([i * 4 / 1_000_000 for i in range(1, stepsToPlot + 1)])
    max_x = max([i * 4 / 1_000_000 for i in range(1, stepsToPlot + 1)])
    
    #plt.axvline(x = 0.032, color = 'k', linestyle = '--')
    #plt.axvline(x = 0.112, color = 'k', linestyle = '--')
    #plt.axvline(x = 0.168, color = 'k', linestyle = '--')
    #plt.axvline(x = 0.24, color = 'k', linestyle = '--')
    #plt.axvline(x = 0.25, color = 'k', linestyle = '--')

    #labels for vertical lines
    #labels = ['S1', 'S2', 'S3', 'S4']
    #labels = ['S1']

    #label positions
    #label_positions = [0.032, 0.112, 0.168, 0.24]
    #label_positions = [0.25]

    #Place labels above figure
    #fig = plt.gcf()
    #ax = plt.gca()

    # Get the current y-axis limits
    #y_lim = ax.get_ylim()

    # Loop through each label and use ax.text()
    #for x_pos, label in zip(label_positions, labels):
        #ax.text(x_pos, y_lim[1], label, horizontalalignment='center', verticalalignment='bottom', fontsize=20)


    
    # Set y-axis ticks and labels
    ticks = np.arange(4, 8.5, 1)
    plt.yticks(ticks)
    plt.gca().yaxis.set_major_locator(ticker.FixedLocator(ticks))
    
    # Format y-axis labels to have no decimal point
    formatter = ticker.FuncFormatter(lambda x, pos: '{:,.0f}'.format(x))
    plt.gca().yaxis.set_major_formatter(formatter)

    plt.yticks(ticks = plt.yticks()[0], labels = plt.yticks()[0].astype(int))

    
    plt.xlim(min_x, max_x)
    plt.xlabel('Simulation Time / ns', fontsize=30)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    plt.ylabel('Distance / Å', fontsize=30)

    #plt.axvline(x = 0.25, color = 'k', linestyle = '--')
    plt.legend()
    plt.show()