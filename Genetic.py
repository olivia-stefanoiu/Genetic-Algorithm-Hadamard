import copy
import math
import secrets
import random
import numpy as np
from MeepGeo import SimulationStrategy

CHROMOSOME_SIZE = 4
INITIAL_CROMO_SIZE = 4


class Chromosome:
    MUTATION_PROBABILITY = 0.001
    GENES = [2.25, 2.62]

    def __init__(self, initial_genes=None):
        self.genes = initial_genes if initial_genes is not None else np.array(
            [secrets.choice(Chromosome.GENES) for _ in range(INITIAL_CROMO_SIZE)])
        self.fitness = -math.inf

    def calculate_fitness(self, ex_data, ey_data):

        sum_x = 0
        sum_y = 0
        for i in range(len(ex_data)):
            for j in range(len(ex_data[i])):
                for k in range(len(ex_data[i][j])):
                    sum_x += pow(ex_data[i][j][k], 2)
                    sum_y += pow(ey_data[i][j][k], 2)

        # Filter elements lower than 0.5
        # ex_data = np.array([x for x in ex_data if x > 0.5])
        if sum_x > sum_y:
            sum_x, sum_y = sum_y, sum_x
        self.fitness = sum_x / sum_y

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
                if self.genes[i] == 2.25:
                    self.genes[i] = 2.62
                else:
                    self.genes[i] = 2.25

    @staticmethod
    def crossover(chromosome_one, chromosome_two):
#TODO add crossover on intervals
        chromosome_size = len(chromosome_one.genes)
        secret = secrets.choice(range(chromosome_size - 1))
        for i in range(secret, chromosome_size):
            aux = chromosome_one.genes[i]
            chromosome_one.genes[i] = chromosome_two.genes[i]
            chromosome_two.genes[i] = aux


class GeneticAlgorithm:
    POPULATION_SIZE = 10 #TODO why on earth is this attribute alone here?

    @staticmethod
    def _get_initial_generation(population_size):
        return [Chromosome() for _ in range(population_size)]

    def __init__(self, population_size=POPULATION_SIZE, chromo_size=CHROMOSOME_SIZE,
                 population: list[Chromosome] = None,
                 generation=2, fitness_strategy=None):

        self.population_size = population_size
        self.chromo_size = chromo_size
        self.population = GeneticAlgorithm._get_initial_generation(
            self.population_size) if population is None else population
        self.generation = generation
        self.fitness_strategy = SimulationStrategy() if (fitness_strategy is None) else fitness_strategy
        self.square_length = 32
        self.CROSSOVERS = [10, 15, 20,40]
        self.sourceFile = open("aux3.txt", "w")

    def __del__(self):
        self.sourceFile.close()

    def _get_percentage(self):  # calculating the fitness percentages
        fitness = [chromo.fitness for chromo in self.population]
        sum_f = np.sum(fitness)

        return [fit / sum_f for fit in fitness]

    def _crossover_candidates(self, precision_factor: int):
        relative_fitnesses = self._get_percentage()
        probability = []
        result = []
        for i in range(len(relative_fitnesses)):
            probability.extend([i] * int(relative_fitnesses[
                                             i] * precision_factor))

        for i in range(self.population_size):
            rand = random.randrange(len(probability))
            chromosomes = self.population[probability[rand]]
            result.append(copy.deepcopy(chromosomes))
        return result

    def save_fitness(self):
        fitness = [chromo.fitness for chromo in self.population]
        self.sourceFile.write(str(fitness))
        self.sourceFile.write("\n")
        self.sourceFile.flush()

    def _update_population_fitness(self):
        ex_data_list, ey_data_list = self.fitness_strategy.run_simulation(self.generation, self.population,
                                                                          self.square_length)
        for i in range(self.population_size):
            self.population[i].calculate_fitness(ex_data_list[i], ey_data_list[i])

    def advance_generation(self):

        if self.population[0].fitness == -math.inf:
            # First time we run the algorithm, need to init the fitness
            self._update_population_fitness()

        for i in range(self.CROSSOVERS[self.generation - 2]):

            # self.population = sorted(self.population, key=lambda x: x.fitness)
            crossover_cand = self._crossover_candidates(self.population_size * 10)

            for j in range(0, self.population_size, 2):
                Chromosome.crossover(crossover_cand[j], crossover_cand[j + 1])
            for j in range(self.population_size):
                crossover_cand[j].mutation()
            # eliminate the worst chromosomes
            self.population = crossover_cand
            # run simulation
            self._update_population_fitness()

        self.generation = self.generation + 1
        exponent = 2 * (self.generation - 1)
        self.chromo_size = pow(2, exponent)

        for i in range(self.population_size):
            self.population[i].split_cells(self.chromo_size, self.generation)

        self.square_length = self.square_length / 2
