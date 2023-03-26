import meep as mp
import h5py
import numpy as np
import matplotlib.pyplot as plt


resolution = 10  # 50 pixels per unit a.
cell = mp.Vector3(20, 20, 0)  # just work in 2D for this

freq = 0.66713
length = 1
index = 1
for i in range(3):

    for j in range(populationSize):
    sources = [mp.Source(mp.ContinuousSource(frequency=freq),
                         center=mp.Vector3(x=-0.5 * length, y=0, z=0),
                         component=mp.Ez)]  # 1mA source amplitude=1.0

    geometry = [mp.Block(mp.Vector3(length, length, 0),  # define an infinite block
                         center=mp.Vector3(xi, yi, 0),  # centered at the origin
                         material=mp.Medium(epsilon=index))
                for xi, yi, index in
                [(0, 0, 1), (0, length, 1), (0, 2 * length, 1), (0, 3 * length, 1)]]  # that is a vacuum

    pml_layers = [mp.PML(1.0)]

    sim = mp.Simulation(cell_size=cell,
                        boundary_layers=pml_layers,
                        geometry=geometry,
                        sources=sources,
                        resolution=resolution)
    #
    # sim.run(mp.to_appended("ez", mp.at_every(0.5, mp.output_efield_z)),
    #         until=10)

    sim.run(until=20)

    eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
    plt.figure()
    plt.imshow(eps_data.transpose(), interpolation="spline36", cmap="binary")
    plt.axis("off")
    plt.show()


# f = h5py.File('MeepGeo-ez.h5', 'r')
#
# efield_z = np.array(f.get('ez')) # a float tensor of shape (600, 600, 350)
#
# plt.figure("ez")
# plt.imshow(efield_z.transpose())
# plt.show()
#
# f.close()
