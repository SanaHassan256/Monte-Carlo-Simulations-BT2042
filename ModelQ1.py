import numpy as np
import matplotlib.pyplot as plt
import random
import matplotlib.animation as animation
import copy
import math

#random_seed = 12  # You can choose any seed value
#np.random.seed(random_seed)
#random.seed(random_seed)

lattice_size = 30

temperature = 1.0 # taking KbT =1


def generate_polymer():  # Function to generate a polymer
    x = np.random.randint(1, lattice_size-3)
    y = np.random.randint(1, lattice_size-3)
    polymer_coordinates = [] # We are Storing the coordinates of the polymer in a list

    polymer_coordinates.append([x, y])
    polymer_coordinates.append([x + 1, y])
    polymer_coordinates.append([x + 2, y])
    polymer_coordinates.append([x + 3, y])

    #  we are generating a  linear polymer of length 4 and returning the coordinates
    # we are storing evreything in a nested list list causs it is easier to access the coordinates of the beads and mutable
    return polymer_coordinates


def calculate_energy_difference(old_energy, new_energy): # Simple function to calculate the energy difference
    return new_energy - old_energy



def calculate_interaction_energy(polymers): # Function to calculate the interaction energy
    total_interaction_energy = 0 # intialing the total interaction energy to 0
    for polymer_name, polymer in polymers.items(): # accessing the polymer name and the coordinates of the polymer
        interaction_energy = 0 # intializing the interaction energy to 0
        for bead_index, bead in enumerate(polymer): # the enumerate function returns the index and the value of the list so the bead index and the bead coordinates are returned
            x, y = bead # the x and y coordinates of the bead are stored in x and y
            for other_polymer_name, other_polymer in polymers.items(): # we are accessing the other polymers
                if polymer_name != other_polymer_name:  #  we are Avoiding  comparing the same polymer so that we can compare 2 polymers
                    other_bead = other_polymer[bead_index]  # we are accessing the bead of the other polymer at the same index for comaparing with orginal bead
                    distance = abs(x - other_bead[0]) + abs(y - other_bead[1]) # calculating the distance between the beads of same type
                    if distance == 1:  # condition for adjuscent beads
                        interaction_energy -=  10# Increment energy for adjacent beads
        total_interaction_energy += interaction_energy # adding the interaction energy to the total interaction energy
    return total_interaction_energy


def metropolis_criterion(delta_energy): # the metroplis criterion function self explanatory
    if delta_energy < 0:
        return True
    if delta_energy == 0:
        return random.choice([True, False])
    else:
        probability = math.exp(-delta_energy / temperature)
        return random.random() < probability


def end_movement1(polymers):   # Function for the end movement of the polymer U can see the writeup for the  intial algo we came up with
    for polymer_name, polymer in polymers.items():
        new_polymers = copy.deepcopy(polymers)
        newcoord = [0, 0]
        for i in range(0, 1):
            if polymer[i][0] != polymer[i + 1][0] and polymer[i][1] == polymer[i + 1][1]:
                newcoord[0] = polymer[i + 1][0]
                newcoord[1] = polymer[i + 1][1] + random.choice([-1, 1])
                if 0 <= newcoord[0] < lattice_size and 0 <= newcoord[1] < lattice_size \
                        and newcoord not in polymer \
                        and not any(newcoord in other_polymer for other_polymer in polymers.values() if
                                    other_polymer != polymer):
                    new_polymers[polymer_name][i][0] = newcoord[0]
                    new_polymers[polymer_name][i][1] = newcoord[1]
                    new_energy = calculate_interaction_energy(new_polymers)
                    old_energy = calculate_interaction_energy(polymers)
                    delta_energy = calculate_energy_difference(old_energy, new_energy)
                    if metropolis_criterion(delta_energy):
                        polymers[polymer_name][i][0] = newcoord[0]
                        polymers[polymer_name][i][1] = newcoord[1]

            if polymer[i][1] != polymer[i + 1][1] and polymer[i][0] == polymer[i + 1][0]:
                newcoord[0] = polymer[i + 1][0] + random.choice([-1, 1])
                newcoord[1] = polymer[i + 1][1]
                if 0 <= newcoord[0] < lattice_size and 0 <= newcoord[1] < lattice_size \
                        and newcoord not in polymer \
                        and not any(newcoord in other_polymer for other_polymer in polymers.values() if
                                    other_polymer != polymer):
                    new_polymers[polymer_name][i][0] = newcoord[0]
                    new_polymers[polymer_name][i][1] = newcoord[1]
                    new_energy = calculate_interaction_energy(new_polymers)
                    old_energy = calculate_interaction_energy(polymers)
                    delta_energy = calculate_energy_difference(old_energy, new_energy)
                    if metropolis_criterion(delta_energy):
                        polymers[polymer_name][i][0] = newcoord[0]
                        polymers[polymer_name][i][1] = newcoord[1]

        last_index = len(polymer) - 1
        if polymer[last_index][0] != polymer[last_index - 1][0] and polymer[last_index][1] == polymer[last_index - 1][
            1]:
            newcoord[0] = polymer[last_index - 1][0]
            newcoord[1] = polymer[last_index - 1][1] + random.choice([-1, 1])
            if 0 <= newcoord[0] < lattice_size and 0 <= newcoord[1] < lattice_size \
                    and newcoord not in polymer \
                    and not any(newcoord in other_polymer for other_polymer in polymers.values() if
                                other_polymer != polymer):
                new_polymers[polymer_name][last_index][0] = newcoord[0]
                new_polymers[polymer_name][last_index][1] = newcoord[1]
                new_energy = calculate_interaction_energy(new_polymers)
                old_energy = calculate_interaction_energy(polymers)
                delta_energy = calculate_energy_difference(old_energy, new_energy)
                if metropolis_criterion(delta_energy):
                    polymers[polymer_name][last_index][0] = newcoord[0]
                    polymers[polymer_name][last_index][1] = newcoord[1]
        if polymer[last_index][1] != polymer[last_index - 1][1] and polymer[last_index][0] == polymer[last_index - 1][
            0]:
            newcoord[0] = polymer[last_index - 1][0] + random.choice([-1, 1])
            newcoord[1] = polymer[last_index - 1][1]
            if 0 <= newcoord[0] < lattice_size and 0 <= newcoord[1] < lattice_size \
                    and newcoord not in polymer \
                    and not any(newcoord in other_polymer for other_polymer in polymers.values() if
                                other_polymer != polymer):
                new_polymers[polymer_name][last_index][0] = newcoord[0]
                new_polymers[polymer_name][last_index][1] = newcoord[1]
                new_energy = calculate_interaction_energy(new_polymers)
                old_energy = calculate_interaction_energy(polymers)
                delta_energy = calculate_energy_difference(old_energy, new_energy)
                if metropolis_criterion(delta_energy):
                    polymers[polymer_name][last_index][0] = newcoord[0]
                    polymers[polymer_name][last_index][1] = newcoord[1]


def cornermovement(polymers): # similarly for the corner moves algo can be found in the writeup
    for polymer_name, polymer in polymers.items():
        new_polymers = copy.deepcopy(polymers) # we create a copy to act as a dumy variable for comparison of interaction energy
        newcoord = [0, 0]
        for i in range(1, 3):
            if polymer[i][0] == polymer[i - 1][0] and polymer[i][0] != polymer[i + 1][0]:
                newcoord[0] = polymer[i - 1][0]
                newcoord[1] = polymer[i + 1][1]
                if 0 <= newcoord[0] < lattice_size and 0 <= newcoord[1] < lattice_size \
                        and newcoord not in polymer \
                        and not any(newcoord in other_polymer for other_polymer in polymers.values() if
                                    other_polymer != polymer):
                    new_polymers[polymer_name][i][0] = newcoord[0]
                    new_polymers[polymer_name][i][1] = newcoord[1]
                    new_energy = calculate_interaction_energy(new_polymers)
                    old_energy = calculate_interaction_energy(polymers)
                    delta_energy = calculate_energy_difference(old_energy, new_energy)
                    if metropolis_criterion(delta_energy):
                        polymers[polymer_name][i][0] = newcoord[0]
                        polymers[polymer_name][i][1] = newcoord[1]

            if polymer[i][1] == polymer[i - 1][1] and polymer[i][1] != polymer[i + 1][1]:
                newcoord[0] = polymer[i - 1][0]
                newcoord[1] = polymer[i + 1][1]
                if 0 <= newcoord[0] < lattice_size and 0 <= newcoord[1] < lattice_size \
                        and newcoord not in polymer \
                        and not any(newcoord in other_polymer for other_polymer in polymers.values() if
                                    other_polymer != polymer):
                    new_polymers[polymer_name][i][0] = newcoord[0]
                    new_polymers[polymer_name][i][1] = newcoord[1]
                    new_energy = calculate_interaction_energy(new_polymers)
                    old_energy = calculate_interaction_energy(polymers)
                    delta_energy = calculate_energy_difference(old_energy, new_energy)
                    if metropolis_criterion(delta_energy):
                        polymers[polymer_name][i][0] = newcoord[0]
                        polymers[polymer_name][i][1] = newcoord[1]


def plot_polymers(polymers):
    fig, ax = plt.subplots()
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_title(f'Polymers (Total Interaction Energy: {calculate_interaction_energy(polymers)})')
    ax.grid(True)
    ax.set_xlim(0, lattice_size)
    ax.set_ylim(0, lattice_size)

    #the first part sets the basic condition for plotting reffererd from the documentaion page from geeks for geeks

    lines = []
    scatters = []
    color_map = ['blue', 'green', 'red', 'purple']  # List of colors for different positions in the polymer corresponding  to the beads

    for polymer_name, polymer in polymers.items(): # accessing the polymer name and the coordinates of the polymer
        x_values = [coord[0] for coord in polymer] # storing the x coordinates of the polymer
        y_values = [coord[1] for coord in polymer] # storing the y coordinates of the polymer
        line, = ax.plot(x_values, y_values, color='black', linewidth=2) #storing the line plot of the polymer Documentaion from official matplotlib page
        lines.append(line)
        scatter_collection = [] # intializing the scatter collection
        for bead_index, (x, y) in enumerate(polymer): # the enumerate function returns the index and the value of the list so the bead index and the bead coordinates are returned
            scatter = ax.scatter(x, y, color=color_map[bead_index]) # Imp Step we are mapping the colors with the beads and the index values
            scatter_collection.append(scatter) # appending the scatter plot to the scatter collection
        scatters.append(scatter_collection) # appending the scatter collection to the scatters list so the list stores the scatter collection of all the polymers

    def update(frame): # the update function for the animation so that we continusly update the plot
        end_movement1(polymers) # the end movement function is called remember the metroplis criterion is implemented in the end movement function
        cornermovement(polymers) # the corner movement function is called remember the metroplis criterion is implemented in the corner movement function
        total_interaction_energy = calculate_interaction_energy(polymers) # the total interaction energy is calculated
        ax.set_title(f'Polymers (Total Interaction Energy: {total_interaction_energy})')
        # the title is updated with the new interaction energy

        # This part is taken from the documenataion page of matplotlib for animation

        for i, (polymer, scatter_collection) in enumerate(zip(polymers.values(), scatters)):
            x_values = [coord[0] for coord in polymer]
            y_values = [coord[1] for coord in polymer]
            lines[i].set_data(x_values, y_values)
            for scatter, coord in zip(scatter_collection, polymer):
                scatter.set_offsets([coord[0], coord[1]])  # Update scatter plot positions

        return lines + [scatter for scatter_collection in scatters for scatter in scatter_collection]

    ani = animation.FuncAnimation(fig, update, frames=1000, interval=1, repeat=False) # the animation function is called
    plt.show()

    # initialising list


ini_list = ['a', 'b', 'c', 'd'] # the list of the beads for the keys for the dictionary final`_dict_beads


def bead_coord(polymers):
    bead_coordinates = {}

    for polymer, coordinates in polymers.items(): # accessing the polymer name and the coordinates of the polymer
        for bead_index, coord in enumerate(coordinates): # the enumerate function returns the index and the value of the list so the bead index and the bead coordinates are returned
            if bead_index not in bead_coordinates: # if the bead index is not in the bead coordinates so this way every coordinate from the polymer gets assigned to the bead index
                bead_coordinates[bead_index] = [] # intializing the bead index to an empty list if it doesnt exist
            bead_coordinates[bead_index].append(coord) # appending the coordinates to the bead index
# the Zip function is used to combine the two lists and create a dictionary
    final_dict_beads = dict(zip(ini_list, list(bead_coordinates.values()))) # the final dictionary is created with the bead index as the key and the coordinates as the values
    return final_dict_beads


polymers = {}
for i in range(6): # creating 6 polymers
    polymers[f'polymer_{i + 1}'] = generate_polymer()
plot_polymers(polymers)

