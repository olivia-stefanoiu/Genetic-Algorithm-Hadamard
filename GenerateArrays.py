import numpy as np

COORD_FACTORS = [8, 4, 2, 1]


def create_coordinates_mat(coord_array, population, generation):
    square_side = int(np.sqrt(len(population[0].genes)))

    coord = []
    x = 0
    coord_mat_array = [
        [(coord_array[i][0], coord_array[i][1], chromo.genes[i])
         for i in range(len(chromo.genes))]
        for chromo in population]

    for i in range(len(population[0].genes)):
        coord.append((coord_array[i][0] - COORD_FACTORS[generation],
                      coord_array[i][1] - COORD_FACTORS[generation]))

        coord.append((coord_array[i][0] + COORD_FACTORS[generation],
                      coord_array[i][1] - COORD_FACTORS[generation]))

        coord.append((coord_array[i][0] - COORD_FACTORS[generation],
                      coord_array[i][1] + COORD_FACTORS[generation]))

        coord.append((coord_array[i][0] + COORD_FACTORS[generation],
                      coord_array[i][1] + COORD_FACTORS[generation]))

    # for i in range(square_side):
    #     x = x + size
    #     y = 0
    #     for j in range(square_side):
    #         coord.append((x, y))
    #         y = y + size

    return coord_mat_array, coord
