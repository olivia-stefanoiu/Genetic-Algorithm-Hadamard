import os
import time
import meep as mp
import h5py
import numpy as np
import matplotlib.pyplot as plt
import GenerateArrays
import multiprocessing


# TODO send code
def run_sim_paralel(square_length, coord_mat_array, j, q):
    cells = mp.Vector3(100, 100, 0)
    pml_layers = [mp.PML(1.0)]

    central_freq = 0.153  # pulse center frequency
    pulse_width = 0.11  # pulse width (in frequency)
    sources = [mp.Source(mp.GaussianSource(central_freq, fwidth=pulse_width),
                         component=mp.Ex,
                         center=mp.Vector3(-37, 0, 0),
                         size=mp.Vector3(0, 4, 0))]
    resolution = 1

    geometry = [mp.Block(mp.Vector3(square_length, square_length, 0),
                         center=mp.Vector3(xi, yi, 0),
                         material=mp.Medium(epsilon=index))
                for xi, yi, index in coord_mat_array[j]]

    sim = mp.Simulation(cell_size=cells,
                        boundary_layers=pml_layers,
                        geometry=geometry,
                        sources=sources,
                        resolution=resolution)

    sim.run(until=170)

    nonpml_vol = mp.Vector3(1, 64)

    ex_data = sim.get_array(center=mp.Vector3(32.5, 0), size=nonpml_vol, component=mp.Ex)
    ey_data = sim.get_array(center=mp.Vector3(32.5, 0), size=nonpml_vol, component=mp.Ey)

    q.put((ex_data, ey_data, j))


class SimulationStrategy:

    def __init__(self):
        # Simulation
        self.coord_array = [(-16, -16), (16, -16), (-16, 16), (16, 16)]

    def run_simulation(self, generation, population, square_length):
        coord_mat_array, self.coord_array = \
            GenerateArrays.create_coordinates_mat(self.coord_array, population, generation)

        ex_data_list = [0 for i in range(len(population))]
        ey_data_list = [0 for i in range(len(population))]

        process_list = []
        ctx = multiprocessing.get_context('fork')
        q = ctx.Queue()

        for j in range(len(population)):

            p = ctx.Process(target=run_sim_paralel,
                                        args=(square_length, coord_mat_array, j, q))
            p.start()
            process_list.append(p)

        for i in range(len(population)):
            ex_data, ey_data, j = q.get()
            ex_data_list[j] = ex_data
            ey_data_list[j] = ey_data

        for p in process_list:
            p.join()

        return ex_data_list, ey_data_list
