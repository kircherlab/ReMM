import random


def getSeedsForTraining(training_run, n):
    """
    Getting n number of integers in the range of 0,100000 using the seed of the given training run.
    """
    seed = config["training"][training_run]["config"]["seed"]
    random.seed(seed)

    seeds = [random.randint(0, 100000) for x in range(n)]

    return seeds
