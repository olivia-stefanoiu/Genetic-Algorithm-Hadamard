import meep as mp
import h5py
import numpy as np
import matplotlib.pyplot as plt
import GenerateArrays
#TODO send code

class SimulationStrategy:

    def __init__(self):
        # Simulation
        self.coord_array = [(-16, -16), (16, -16), (-16, 16), (16, 16)]
        self.cells = mp.Vector3(100, 100, 0)  # just work in 2D for this
        self.pml_layers = [mp.PML(1.0)]

        central_freq = 0.3  # pulse center frequency
        pulse_width = 30  # pulse width (in frequency)
        self.sources = [mp.Source(mp.GaussianSource(central_freq, fwidth=pulse_width),
                                  component=mp.Ey,
                                  center=mp.Vector3(-37, 0, 0),
                                  size=mp.Vector3(0, 4, 0))]
        self.resolution = 1

    # TODO calculate the desired freq in meep units
    # TODO add proper materials
    # TODO add border around the grid to stop radiation from escaping

    # sources = [mp.EigenModeSource(src=mp.GaussianSource(central_freq,fwidth=pulse_width),
    #                               center=mp.Vector3(-37,0),
    #                               size=mp.Vector3(0,3),
    #                               eig_parity=mp.EVEN_Y+mp.ODD_Z,
    #                               eig_band=1,
    #                               eig_match_freq=True)]

    def run_simulation(self, generation, population, square_length):
        coord_mat_array, self.coord_array = \
            GenerateArrays.create_coordinates_mat(self.coord_array, population, generation)
        print(self.coord_array)

        ex_data_list = []
        ey_data_list = []

        for j in range(len(population)):
            geometry = [mp.Block(mp.Vector3(square_length, square_length, 0),
                                 center=mp.Vector3(xi, yi, 0),
                                 material=mp.Medium(epsilon=index))
                        for xi, yi, index in coord_mat_array[j]]

            geometry.append(mp.Block(mp.Vector3(10, 4, 0),
                                     center=mp.Vector3(-37, 0, 0),
                                     material=mp.Medium(epsilon=3.9)))

            geometry.append(mp.Block(mp.Vector3(64, 6, 0),
                                     center=mp.Vector3(0, -35, 0),
                                     material=mp.Medium(epsilon=1.4)))
            geometry.append(mp.Block(mp.Vector3(64, 6, 0),
                                     center=mp.Vector3(0, 35, 0),
                                     material=mp.Medium(epsilon=1.4)))
            # geometry.append(mp.Block(mp.Vector3(6,76, 0),
            #                          center=mp.Vector3(35,0, 0),
            #                          material=mp.Medium(epsilon=3.9)))

            sim = mp.Simulation(cell_size=self.cells,
                                boundary_layers=self.pml_layers,
                                geometry=geometry,
                                sources=self.sources,
                                resolution=self.resolution)

            sim.run(until=170)

            eps_data = sim.get_array(center=mp.Vector3(), size=self.cells, component=mp.Ey)
            # print (eps_data)
            #
            # plt.figure()
            # plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
            # plt.axis('on')
            # plt.show()

            nonpml_vol = mp.Vector3(1, 64)
            # ex_data = sim.get_efield_x(cell_list=nonpml_vol)
            # ey_data = sim.get_efield_x(cell_list=nonpml_vol)

            # print(len(ex_data))
            # print(len(ey_data))

            ex_data = sim.get_array(center=mp.Vector3(-32.5, 0), size=nonpml_vol, component=mp.Ex)
            ey_data = sim.get_array(center=mp.Vector3(-32.5, 0), size=nonpml_vol, component=mp.Ey)
            # TODO export arrays
            # print(len(ex_data))
            # print(len(ey_data))
            # print(ex_data)
            # print(
            #      "--------------------------------------------------------------------------------------------------------------------"
            #  )
            # print(ey_data)

            ex_data_list.append(ex_data)
            ey_data_list.append(ey_data)

            # plt.figure()
            # plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9)
            # plt.axis('off')
            # plt.show()

            # plt.figure()
            # plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
            # plt.imshow(ex_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9)
            # plt.axis('off')
            # plt.show()
        return ex_data_list, ey_data_list
        # coord_mat_array, coord_array = GenerateArrays.create_coordinates_mat(coord_array, algo.population, i)
        # print(coord_array)
        # for j in range(algo.POPULATION_SIZE):
        #     geometry = [mp.Block(mp.Vector3(length, length
        #     , 0),
        #                          center=mp.Vector3(xi, yi, 0),
        #                          material=mp.Medium(epsilon=index))
        #                 for xi, yi, index in coord_mat_array[j]]
        #
        #     geometry.append(mp.Block(mp.Vector3(10, 4, 0),
        #                              center=mp.Vector3(-37, 0, 0),
        #                              material=mp.Medium(epsilon=3.9)))
        #
        #     geometry.append(mp.Block(mp.Vector3(64, 6, 0),
        #                              center=mp.Vector3(0, -35, 0),
        #                              material=mp.Medium(epsilon=1.4)))
        #     geometry.append(mp.Block(mp.Vector3(64, 6, 0),
        #                              center=mp.Vector3(0, 35, 0),
        #                              material=mp.Medium(epsilon=1.4)))
        #     # geometry.append(mp.Block(mp.Vector3(6,76, 0),
        #     #                          center=mp.Vector3(35,0, 0),
        #     #                          material=mp.Medium(epsilon=3.9)))
        #
        #     sim = mp.Simulation(cell_size=cell,
        #                         boundary_layers=pml_layers,
        #                         geometry=geometry,
        #                         sources=sources,
        #                         resolution=resolution)
        #
        #     sim.run(until=180)
        #
        #     eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ey)
        #     # print (eps_data)
        #     #
        #     # plt.figure()
        #     # plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
        #     # plt.axis('on')
        #     # plt.show()
        #
        #     nonpml_vol = mp.Vector3(1, 64)
        #     # ex_data = sim.get_efield_x(cell_list=nonpml_vol)
        #     # ey_data = sim.get_efield_x(cell_list=nonpml_vol)
        #
        #     # print(len(ex_data))
        #     # print(len(ey_data))
        #
        #     ex_data = sim.get_array(center=mp.Vector3(-32.5, 0), size=nonpml_vol, component=mp.Ex)
        #     ey_data = sim.get_array(center=mp.Vector3(-32.5, 0), size=nonpml_vol, component=mp.Ey)
        #     #TODO export arrays
        #     # print(len(ex_data))
        #     # print(len(ey_data))
        #     # print(ex_data)
        #     # print(
        #     #      "--------------------------------------------------------------------------------------------------------------------"
        #     #  )
        #     # print(ey_data)
        #
        #
        #
        # algo.population[j].calculate_fitness( ex_data, ey_data)




            # plt.figure()
            # plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9)
            # plt.axis('off')
            # plt.show()

            # plt.figure()
            # plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
            # plt.imshow(ex_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9)
            # plt.axis('off')
            # plt.show()




