import sys
import argparse

def main(args):
    # The input is assumed to be 3 columns and of the form:
    # name    length    estimated_read_count
    # the last column should be convertable to a float.
    exps = None
    try:
        exps = [ (l.split()[0], float(l.split()[1]), float(l.split()[2])) for l in args.infile]
    except ValueError:
        print "Could not convert entry in column 2 or 3 to a float"

    outfile = args.outfile
    tot = sum([e[2] for e in exps])
    denom = tot / 10**6
    
    for e in exps:
        l = e[1]
        rc = e[2]
        rpkm = (rc / (l / 1000.0)) / denom if l > 0.0 else 0.0
        outfile.write("{}\t{}\t{}\n".format(e[0], l, rpkm))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute RPKM from transcript length and read count.")
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),\
                         default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),\
                         default=sys.stdout)
    args = parser.parse_args()
    main(args)
