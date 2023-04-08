import meep as mp
import h5py
import numpy as np
import matplotlib.pyplot as plt
import Genetic
import GenerateArrays


resolution = 2
cell = mp.Vector3(100, 100, 0)  # just work in 2D for this
#TODO calculate the desired freq in meep units
#TODO add proper materials
#TODO add border around the grid to stop radiation from escaping

length = 32
coord_array=[(-16,-16),(16,-16),(-16,16),(16,16)]

algo = Genetic.GeneticAlgorithm()


central_freq = 0.3 # pulse center frequency
pulse_width = 70   # pulse width (in frequency)
sources = [mp.Source(mp.GaussianSource(central_freq,fwidth=pulse_width),
                     component=mp.Ez,
                     center=mp.Vector3(-37,0,0),
                     size=mp.Vector3(0,4,0))]

pml_layers = [mp.PML(1.0)]

coord_mat_array=[]

for i in range(4):   # how many gens we want

    coord_mat_array,coord_array = GenerateArrays.create_coordinates_mat(coord_array,algo.population,i)
    print(coord_array)
    for j in range(algo.POPULATION_SIZE):

        geometry = [mp.Block(mp.Vector3(length, length, 0),
                             center=mp.Vector3(xi, yi, 0),
                             material=mp.Medium(epsilon=index))
                    for xi, yi, index in coord_mat_array[j]]

        geometry.append(mp.Block(mp.Vector3(10, 4, 0),
                             center=mp.Vector3(-37, 0, 0),
                             material=mp.Medium(epsilon=3.9)))

        sim = mp.Simulation(cell_size=cell,
                            boundary_layers=pml_layers,
                            geometry=geometry,
                            sources=sources,
                            resolution=resolution)

        sim.run(until=180)

        eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)

        # plt.figure()
        # plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
        # plt.axis('on')
        # plt.show()

        ez_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
        plt.figure()
        plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
        plt.imshow(ez_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9)
        plt.axis('off')
        plt.show()

    length = length/2
    algo.advance_generation()



