#!/usr/bin/env python3

from random import randint, sample, choice
from itertools import product
import sys

#number of terminals range
t_lim = (2000, 2001)
#number of routers range
r_lim = (1500, 2000)
#capacity range
cap_lim = (1, 40)
#cost range
cost_lim = (1, 50)

def gen(write_to_file=False):
    #number of terminals, routers, edges
    n_t = randint(*t_lim)
    n_r = randint(*r_lim)
    n_e = randint(min(n_t, n_r), n_t*n_r)

    if write_to_file:
        f = open("{}_{}.in".format(n_t, n_r), "w")
    else:
        f = sys.stdout

    #generating random capacities
    t_caps = [randint(*cap_lim) for __ in range(n_t)]
    r_caps = [randint(*cap_lim) for __ in range(n_r)]
    diff = sum(t_caps) - sum(r_caps)
    smaller = t_caps if diff < 0 else r_caps
    for __ in range(abs(diff)):
        smaller[randint(0, len(smaller)-1)] += 1

    #terminals/routers names and capacities
    t = {"t{}".format(i+1): cap for i, cap in enumerate(t_caps)}
    r = {"r{}".format(i+1): cap for i, cap in enumerate(r_caps)}

    #printing example
    print("{} {} {}".format(n_t, n_r, n_e), file=f)
    for k, v in t.items():
        print("{} {}".format(k, v), file=f)
    for k, v in r.items():
        print("{} {}".format(k, v), file=f)
    for u, v in sample(list(product(t.keys(), r.keys())), n_e):
        print("{} {} {}".format(u, v, randint(*cost_lim)), file=f)

def main():
    gen(write_to_file=True)

if __name__ == "__main__":
    main()
