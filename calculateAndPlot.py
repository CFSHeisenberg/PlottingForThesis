from concurrent import futures
from tkinter import Entry
from matplotlib import pyplot as plt
import numpy as np


class DistancePlotter:
    def __init__(self, merged_df_values, numAtoms, centroid_indices_j, guest_indices, numOfCentroids, lattice_values):
        self.merged_df_values = merged_df_values
        self.numAtoms = numAtoms
        self.centroid_indices_j = centroid_indices_j
        self.guest_indices = guest_indices
        self.numOfCentroids = numOfCentroids
        self.lattice_values = lattice_values

    def imaging(self, coords, previous_coords, i):
        difference = np.array(coords - previous_coords)
        difference_new = difference - np.round(np.divide(difference, self.lattice_values).astype(float)) * self.lattice_values 
        coords_changed = previous_coords + difference_new
        return coords_changed

    def calculate_distance(self, i, selected_centroids):
        if i == 1:  
            # Get partial numpy array of current step. Takes only the rows of the current step
            stepframe_values = self.merged_df_values[(i-1)*self.numAtoms:i*self.numAtoms, :]

            # Get partial numpy array of shape (108,7) only containing those entries of stepframe_values whose fifth element, i.e. the cleanedIndex are in centroid_indices_j
            stepframe_j_values = stepframe_values[np.isin(stepframe_values[:, 5], self.centroid_indices_j)]
            print(self.centroid_indices_j)

            # Get coordinates of the atoms making up centroids in current step
            centroid_coords = stepframe_j_values[:, 1:4]
            print(centroid_coords)

            # Get partial numpy array of shape (6, 7) only containing those entries of stepframe_values whose fifth element, i.e. the cleanedIndex are in guest_indices
            stepframe_g_values = stepframe_values[np.isin(stepframe_values[:, 5], self.guest_indices)]

            # Get coordinates of the relevant guest atoms in current step
            guest_coords = stepframe_g_values[:, 1:4]
            print(guest_coords)

            # Get coordinates of cartesian centroids for both host and guest via the mean of their coordinates
            center_guest_coords = np.mean(guest_coords, axis=0)
            center_host_coords = np.mean(np.split(centroid_coords, self.numOfCentroids), axis=1)

            # Calculate distance between each centroid to the guest centroid
            distance_to_center = np.linalg.norm((center_guest_coords - center_host_coords).astype(float), axis=1)

        else:
            # Get partial numpy array of current step. Takes only the rows of the current step
            stepframe_values = self.merged_df_values[(i-1)*self.numAtoms:i*self.numAtoms, :]

            # Get partial numpy array of previous step. Takes only the rows of the previous step
            stepframe_values_previous = self.merged_df_values[(i-2)*self.numAtoms:(i-1)*self.numAtoms, :]

            # Get partial numpy array of shape (108,7) of the current step only containing those entries of stepframe_values whose fifth element, i.e. the cleanedIndex are in centroid_indices_j
            stepframe_j_values = stepframe_values[np.isin(stepframe_values[:, 5], self.centroid_indices_j)]

            # Get partial numpy array of shape (108,7) of the previous step only containing those entries of stepframe_values_next whose fifth element, i.e. the cleanedIndex are in centroid_indices_j
            stepframe_j_values_previous = stepframe_values_previous[np.isin(stepframe_values_previous[:, 5], self.centroid_indices_j)]       

            # Get coordinates of the atoms making up centroids in current step
            centroid_coords = stepframe_j_values[:, 1:4]
            # Get coordinates of the atoms making up centroids in previous step
            centroid_coords_previous = stepframe_j_values_previous[:, 1:4]

            # Get partial numpy array of shape (6, 7) only containing those entries of stepframe_values whose fifth element, i.e. the cleanedIndex are in guest_indices
            stepframe_g_values = stepframe_values[np.isin(stepframe_values[:, 5], self.guest_indices)]

            # Get coordinates of the relevant guest atoms in current step
            guest_coords = stepframe_g_values[:, 1:4]

            # Compare centroid coordinates in next step to centroid coordinates in current step, image them if necessary
            centroid_coords_imaged = self.imaging(centroid_coords, centroid_coords_previous, i)

            # Overwrite the values in the dataframe accessed by stepframe_j_values with centroid_coords_imaged for the next step
            self.merged_df_values[(i-1)*self.numAtoms:i*self.numAtoms, 1:4][np.isin(stepframe_values[:,5], self.centroid_indices_j)] = centroid_coords_imaged

            # Get coordinates of cartesian centroids for both host and guest centroids via the mean of their coordinates
            center_guest_coords = np.mean(guest_coords, axis=0)
            center_host_coords_imaged = np.mean(np.split(centroid_coords_imaged, self.numOfCentroids), axis=1)

            # Calculate distance between each centroid to the guest centroid
            distance_to_center = np.linalg.norm((center_guest_coords - center_host_coords_imaged).astype(float), axis=1)

        # Return the calculated distances for selected centroids
        return [distance_to_center[j-1] for j in selected_centroids]

    def plot_distances(self, stepsToPlot, centroid_vars):
        # Get the total steps to plot from the entry field
        total_steps = stepsToPlot
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
