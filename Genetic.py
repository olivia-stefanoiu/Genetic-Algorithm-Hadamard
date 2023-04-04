import copy
import secrets
import random
import numpy as np
from GenerateArrays import create_coordinates_mat

CHROMOSOME_SIZE = 4
INITIAL_CROMO_SIZE = 4
NUMBER_GENERATIONS = 2


class Chromosome:
    MUTATION_PROBABILITY = 0.001
    GENES = [5, 12]

    def __init__(self, initial_genes=None):
        self.genes = initial_genes if initial_genes is not None else np.array(
            [secrets.choice(Chromosome.GENES) for _ in range(INITIAL_CROMO_SIZE)])
        self.fitness = self.get_fitness()

    def get_fitness(self):
        # fitness = 0.5 - (Target_x - polarisation_x) + 0.5 - (Target_y - polarization_y)
        # TODO: Calculate fitness
        #return np.random.uniform(0.01, 0.98)
        return random.uniform(0.01, 0.98)

    # def newGeneration(self, iterator):
    #     global chromosome_size
    #     new_genes = np.array([5] * chromosome_size)
    #     nr = 0
    #     i = 0
    #     while i < chromosome_size:
    #         new_genes[i] = self.genes[nr]
    #         new_genes[i + 1] = self.genes[nr]
    #         new_genes[i + pow(2, self.generation - 1)] = self.genes[nr]
    #         new_genes[i + pow(2, self.generation - 1) + 1] = self.genes[nr]
    #         if (nr + 1) % iterator == 0:
    #             i = i + pow(2, self.generation) - 2  # !!!
    #         else:
    #             i = i + 2
    #         nr = nr + 1
    #
    #     self.genes = new_genes

    def split_cells(self, new_chromosome_size, new_generation):
        new_genes = np.array([None] * new_chromosome_size)
        nr = 0
        i = 0
        while i < new_chromosome_size:
            new_genes[i] = self.genes[nr]
            new_genes[i + 1] = self.genes[nr]
            new_genes[i + pow(2, new_generation - 1)] = self.genes[nr]
            new_genes[i + pow(2, new_generation - 1) + 1] = self.genes[nr]
            i = i + 2
            while i < len(new_genes) and new_genes[i] is not None:
                i = i + 1
            nr = nr + 1

        self.genes = new_genes

    def mutation(self):
        for i in range(len(self.genes)):
            prob = random.uniform(0, 1)
            if prob <= Chromosome.MUTATION_PROBABILITY:
                # we have only 2 values at the moment, so I just flip the value, else I use random choice
                self.genes[i] = not self.genes[i]

    @staticmethod
    def crossover(chromosome_one, chromosome_two):
        # TODO: no crossover with self, should I though?
        chromosome_size = len(chromosome_one.genes)
        secret = secrets.choice(range(chromosome_size - 1))
        for i in range(secret, chromosome_size):
            aux = chromosome_one.genes[i]
            chromosome_one.genes[i] = chromosome_two.genes[i]
            chromosome_two.genes[i] = aux


class GeneticAlgorithm:
    Target_x = 0.5
    Target_y = 0.5
    POPULATION_SIZE = 100
    @staticmethod
    def _get_initial_generation(population_size):
        return [Chromosome() for _ in range(population_size)]

    def __init__(self, population_size=POPULATION_SIZE, chromo_size=CHROMOSOME_SIZE,
                 population: list[Chromosome] = None,
                 generation=2):
        self.population_size = population_size
        self.chromo_size = chromo_size
        self.population = GeneticAlgorithm._get_initial_generation(
            self.population_size) if population is None else population
        self.generation = generation

    def _get_precision_factor(self):  # calculate the amount to multiply to get
        populationSize = self.population_size
        aux = populationSize
        nr = 0
        while aux > 9:
            aux = aux / 10
            nr = nr + 1
        return pow(10, nr)

    def _get_percentage(self):  # calculating the fitness percentages
        fitness = [chromo.fitness for chromo in self.population]
        sum = np.sum(fitness)

        return [fit / sum for fit in fitness]

    def _crossover_candidates(self, precision_factor: int):
        relative_fitnesses = self._get_percentage()
        probability = []
        result = []
        for i in range(len(relative_fitnesses)):
            probability.extend([i] * int(relative_fitnesses[i] * precision_factor))  # !!! does not reach 100 sadly
        for i in range(self.population_size):
            rand = random.randrange(len(probability))
            chromosomes = self.population[probability[rand]]
            result.append(copy.deepcopy(chromosomes))
        return result

    def advance_generation(self):

        precision_factor = self._get_precision_factor()

        self.population = sorted(self.population, key=lambda x: x.fitness)
        crossover_cand = self._crossover_candidates(precision_factor)

        for i in range(0, self.population_size, 2):
            Chromosome.crossover(crossover_cand[i], crossover_cand[i + 1])
        for i in range(self.population_size):
            crossover_cand[i].mutation()

        # TODO: Recalculate fitness
        # print([c.fitness for c in crossover_cand])
        # print([c.genes for c in crossover_cand])
        self.population = crossover_cand

        self.generation = self.generation+1
        exponent = 2 * (self.generation - 1)
        self.chromo_size = pow(2, exponent)

        for i in range(self.population_size):
            self.population[i].split_cells(self.chromo_size, self.generation)
            #print(self.population[i].genes)

            # aux[i].newGeneration(iterator) experiencing mathematical issues
            # print(aux[i].genes)

        # calculate the number of chromosomes on a row


# print([c.genes for c in population])


# def split_cell_test():
#     global chromosome_size
#     chromosome_size = 16
#     c = Chromosome(np.array([1, 0, 1, 0]))
#     c.split_cells()
#     print(c.genes)


# split_cell_test()
#TODO: add sort, and fitness for each generation
algo = GeneticAlgorithm()
for i in range(3):
    algo.advance_generation()
