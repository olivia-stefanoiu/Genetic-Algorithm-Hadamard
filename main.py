# Pana aici

from Genetic import GeneticAlgorithm


# Asta ar foi main-uul tau. L-am lasat aici ca am presupus ca asta era fisierul pe care-l rulai
def main():
    algo = GeneticAlgorithm()

    for i in range(4):  # how many gens we want
        algo.advance_generation()
        algo.save_fitness()


main()
