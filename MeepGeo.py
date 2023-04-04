import meep as mp
import h5py
import numpy as np
import matplotlib.pyplot as plt
import Genetic
import GenerateArrays


resolution = 5
cell = mp.Vector3(20, 20, 0)  # just work in 2D for this
#TODO calculate the desired freq in meep units
#TODO add proper materials
#TODO add border around the grid to stop radiation from escaping
freq = 0.66713
length = 5


algo = Genetic.GeneticAlgorithm()

sources = [mp.Source(mp.ContinuousSource(frequency=freq),
                             center=mp.Vector3(x=-0.5 * length, y=0, z=0),
                             component=mp.Ez)]  # 1mA source amplitude=1.0

pml_layers = [mp.PML(1.0)]

coord_mat_array=[]

for i in range(2):   # how many gens we want

    for j in range(algo.POPULATION_SIZE):

        coord_mat_array= GenerateArrays.create_coordinates_mat(1, algo.population)


        geometry = [mp.Block(mp.Vector3(length, length, 0),  # define an infinite block
                             center=mp.Vector3(xi, yi, 0),  # centered at the origin
                             material=mp.Medium(epsilon=index))
                    for xi, yi, index in coord_mat_array[j]]

        sim = mp.Simulation(cell_size=cell,
                            boundary_layers=pml_layers,
                            geometry=geometry,
                            sources=sources,
                            resolution=resolution)

        sim.run(until=5)

        # eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
        # plt.figure()
        # plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
        # plt.axis('off')
        # plt.show()



    algo.advance_generation()



