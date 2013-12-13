import sys
import argparse

def main(args):
    # The input is assumed to be 3 columns and of the form:
    # name    length    expression
    # the last column should be convertable to a float.
    exps = None
    try:
        exps = [ (l.split()[0], float(l.split()[1]), float(l.split()[2])) for l in args.infile]
    except ValueError:
        print "Could not convert entry in column 2 or 3 to a float"

    # The total expression is just the sum of all expressions
    tot = sum([e[2] for e in exps])
    norm = 1.0 / tot

    outfile = args.outfile

    # Write the resulting expressed fractions to the output file
    for e in exps:
        l = float(e[1])
        rc = float(e[2])
        frac = rc * norm
        outfile.write("{}\t{}\t{}\n".format(e[0], l, frac))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Normalize last column of input file to sum to 1")
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),\
                         default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),\
                         default=sys.stdout)
    args = parser.parse_args()
    main(args)
