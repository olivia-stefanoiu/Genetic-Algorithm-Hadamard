import copy
import secrets
import random
import numpy as np

Genes = [0, 1]
Target_x = 0.5
Target_y = 0.5
populationSize = 10
generation = 2
chromosome_size = 4
INITIAL_CROMO_SIZE = 4
MUTATION_PROBABILITY = 0.001


class Chromosome:
    def __init__(self, initial_genes = None):
        self.genes = initial_genes if initial_genes is not None else np.array([secrets.choice(Genes) for _ in range(INITIAL_CROMO_SIZE)])
        self.fitness = self.get_fitness()

    def get_fitness(self):
        # fitness = 0.5 - (Target_x - polarisation_x) + 0.5 - (Target_y - polarization_y)
        # TODO: Calculate fitness
        return random.uniform(0.01, 0.99)

    def newGeneration(self, iterator):
        global chromosome_size
        new_genes = np.array([5]*chromosome_size)
        nr = 0
        i = 0
        while i < chromosome_size:
            new_genes[i] = self.genes[nr]
            new_genes[i + 1] = self.genes[nr]
            new_genes[i + pow(2, generation-1)] = self.genes[nr]
            new_genes[i + pow(2, generation-1) + 1] = self.genes[nr]
            if (nr+1) % iterator == 0:
                i = i + pow(2, generation)-2  # !!!
            else:
                i = i + 2
            nr = nr + 1

        self.genes = new_genes

    def split_cells(self):
        global chromosome_size
        new_genes = np.array([None] * chromosome_size)
        nr = 0
        i = 0
        while i < chromosome_size:
            new_genes[i] = self.genes[nr]
            new_genes[i + 1] = self.genes[nr]
            new_genes[i + pow(2, generation - 1)] = self.genes[nr]
            new_genes[i + pow(2, generation - 1) + 1] = self.genes[nr]
            i = i + 2
            while i < len(new_genes) and new_genes[i] is not None:
                i = i+1
            nr = nr + 1

        self.genes = new_genes

    def mutation(self):
        global chromosome_size
        global MUTATION_PROBABILITY
        for i in range(chromosome_size):
            prob = random.uniform(0, 1)
            if prob <= MUTATION_PROBABILITY:
                # we have only 2 values at the moment, so I just flip the value, else I use random choice
                self.genes[i] = not self.genes[i]


def get_percentage(population):  # calculating the fitness percentages
    fitness = [chromo.fitness for chromo in population]
    sum = np.sum(fitness)

    return [fit / sum for fit in fitness]


def crossover_candidates(population):
    relative_fitnesses = get_percentage(population)

    probability = []
    result = []

    for i in range(len(relative_fitnesses)):
        probability.extend([i] * int(relative_fitnesses[i] * 100.0))  # !!! does not reach 100 sadly

    for _ in range(populationSize):
        secret = secrets.choice(range(len(probability)))
        rel_fit = probability[secret]
        chromosomes = population[rel_fit]
        result.append(copy.deepcopy(chromosomes))

    return result


def crossover(chromosome_one, chromosome_two):
    # TODO: no crossover with self
    global chromosome_size
    secret = secrets.choice(range(chromosome_size - 1))
    for i in range(secret, chromosome_size):
        aux = chromosome_one.genes[i]
        chromosome_one.genes[i] = chromosome_two.genes[i]
        chromosome_two.genes[i] = aux


def main():
    global generation
    global chromosome_size
    # generate the first gen of chromosomes, sort, crossover, mutate

    population = [Chromosome() for _ in range(populationSize)]
    population = sorted(population, key=lambda x: x.fitness)
    crossover_cand = crossover_candidates(population)
    for i in range(0, populationSize, 2):
        crossover(crossover_cand[i], crossover_cand[i + 1])
    for i in range(populationSize):
        crossover_cand[i].mutation()
    print([c.genes for c in crossover_cand])
    population = crossover_cand
    aux=population
    generation = generation + 1

    for _ in range(2):
        iterator = np.sqrt(chromosome_size)
        exponent = 2 * (generation - 1)
        chromosome_size = pow(2, exponent)
        for i in range(populationSize):
            population[i].split_cells()
            print(population[i].genes)
            # aux[i].newGeneration(iterator)
            # print(aux[i].genes)
        generation = generation + 1


    # print([c.genes for c in population])

# def split_cell_test():
#     global chromosome_size
#     chromosome_size = 16
#     c = Chromosome(np.array([1, 0, 1, 0]))
#     c.split_cells()
#     print(c.genes)
#

main()
# split_cell_test()