def create_coordinates_mat(iterator, size, population):
    array = []
    coord = []
    aux = []
    x = 0
    for i in range(pow(iterator, 2)):
        x = x + size
        y = 0
        for j in range(iterator):
            coord.append((x, y))
            y = y + size

    for i in range(len(population)):
        aux.append(zip(coord, population[i]))

    return aux

