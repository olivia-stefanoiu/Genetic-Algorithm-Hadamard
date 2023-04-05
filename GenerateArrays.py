import numpy as np


def create_coordinates_mat(size, population):
    square_side = int(np.sqrt(len(population[0].genes)))

    coord = []
    x = 0

    for i in range(square_side):
        x = x + size
        y = 0
        for j in range(square_side):
            coord.append((x, y))
            y = y + size

    return [  #How to format nestedd list comprehensions?
        [(coord[i][0], coord[i][1], chromo.genes[i])
         for i in range(len(chromo.genes))]
        for chromo in population
    ]
