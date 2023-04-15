import numpy as np

COORD_FACTORS = [8, 4, 2, 1]
GEN = 2

def create_coordinates_mat(coord_array, population, generation):
    square_side = int(np.sqrt(len(population[0].genes)))
    global GEN
    coord = []

    if GEN != generation:
        GEN = generation
        for i in range(int((len(population[0].genes)/4))):
            coord.append((coord_array[i][0] - COORD_FACTORS[generation-3],
                          coord_array[i][1] - COORD_FACTORS[generation-3]))

            coord.append((coord_array[i][0] + COORD_FACTORS[generation-3],
                          coord_array[i][1] - COORD_FACTORS[generation-3]))

            coord.append((coord_array[i][0] - COORD_FACTORS[generation-3],
                          coord_array[i][1] + COORD_FACTORS[generation-3]))

            coord.append((coord_array[i][0] + COORD_FACTORS[generation-3],
                          coord_array[i][1] + COORD_FACTORS[generation-3]))


    else:
        coord = coord_array

    coord_mat_array = [
        [(coord[i][0], coord[i][1], chromo.genes[i])
         for i in range(len(chromo.genes))]
        for chromo in population]
    # for i in range(square_side):
    #     x = x + size
    #     y = 0
    #     for j in range(square_side):
    #         coord.append((x, y))
    #         y = y + size

    return coord_mat_array, coord
