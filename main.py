

from Genetic import GeneticAlgorithm

# multiprocessing.set_start_method('spawn')



def main():
    algo = GeneticAlgorithm()
    for i in range(4):  # how many gens we want
        algo.advance_generation()
        algo.save_fitness()


main()
