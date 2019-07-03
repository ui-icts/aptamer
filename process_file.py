import itertools
import RNA

DISTANCE_FILE = "data/14903-Hoinka/distance.txt"


def run_loop():
    count_a = 0
    count_b = 0
    with open(DISTANCE_FILE, "r") as in_file:
        pairs = itertools.combinations(in_file, 2)
        for a, b in pairs:
            count_a = RNA.make_tree(RNA.expand_Full(a))
            count_b = RNA.make_tree(RNA.expand_Full(b))
            RNA.free_tree(count_a)
            RNA.free_tree(count_b)

    print("DONE")


run_loop()
