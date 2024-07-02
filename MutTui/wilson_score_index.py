import argparse
from math import sqrt

def wilsonSI(s, pr, n):
    z = 1.96

    if n == 0:
        return(0)
    
    if s:
        if int(s) == 0:
            return(0)
        p = float(s)/n
    else:
        if float(pr) == float(0):
            return(0)
        p = float(pr)

    denominator = 1 + z**2/n
    centre_adjusted_probability = p + z*z / (2*n)
    adjusted_standard_deviation = sqrt((p*(1 - p) + z*z / (4*n)) / n)
    lower_bound = (centre_adjusted_probability - z*adjusted_standard_deviation) / denominator
    upper_bound = (centre_adjusted_probability + z*adjusted_standard_deviation) / denominator

    if s:
        print((lower_bound * n), (upper_bound * n))
    else:
        print(lower_bound, upper_bound)

parser = argparse.ArgumentParser()

ns = parser.add_mutually_exclusive_group(required = True)

ns.add_argument("-s", help = "Number of sucesses")
ns.add_argument("-p", help = "Proportion of sucesses")

parser.add_argument("-n", help = "Number of trials")
args = parser.parse_args()

wilsonSI(args.s, args.p, float(args.n))