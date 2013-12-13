"""PlotCorrelation

Usage: PlotCorrelation.py EXP1 --name=<name> [--norm=<norm>] [--m1=<min1>] EXP2 --name=<name> [--norm=<norm>] [--m2=<min2>] [--scatter] [--out=<out>] [--noylabel] [--noplot] [--strip_accession]

Arguments:
  EXP1     The first expression file
  EXP2     The second expression file

Options:
  -h                     Display help
  --name=<name>          The name of the expression estimation method
  --norm=<norm>          The normalization to perform
  --m1=<min1>            The minimum expression to consider from the first dataset [default: -inf]
  --m2=<min2>            The minimum expression to consider from the second dataset [default: -inf]
  --noplot               If this option is present, no plot is drawn or written to file, only the stdout text is produced
  --out=<out>            Save the plot to the specified file rather than draw it.
  --scatter              If set, draw a scatter plot, otherwise draw a hexbin
  --noylabel             Don't draw the y-axis label
  --strip_accession      Strip variant information from transcript identifiers
"""
from docopt import docopt

from matplotlib import pyplot as plt
from matplotlib import rcParams
from mpltools import style
from scipy import stats
import scipy as sp
import numpy as np
import ExpressionTools
import itertools

def correlationInRange(x, y, srange):
    '''
    Find the correlation of all scores
    where x is in the range of [score[0], score[1])

    Args:
        x: set of scores
        y: set of scores in 1-1 correspondence with x
        srange: tuple consiting of (minS, maxS), which are the minimum
                and maximum scores in x which should be considered
    '''

    import bisect
    minS, maxS = srange

    minI = bisect.bisect_left(x, minS)
    maxI = bisect.bisect_left(x, maxS)

    xs = x[minI:maxI]
    ys = y[minI:maxI]

    pr = sp.stats.pearsonr(xs, ys)[0]
    sr = sp.stats.spearmanr(xs, ys)[0]

    return pr, sr


def differencesInRange(x, y, srange, relative=True):
    '''
    Find the relative differences of all scores
    where x is in the range of [score[0], score[1])

    Args:
        x: set of scores
        y: set of scores in 1-1 correspondence with x
        srange: tuple consiting of (minS, maxS), which are the minimum
                and maximum scores in x which should be considered
    '''

    import bisect
    minS, maxS = srange

    minI = bisect.bisect_left(x, minS)
    maxI = bisect.bisect_left(x, maxS)

    xs = np.array(x[minI:maxI])
    ys = np.array(y[minI:maxI])

    return ((xs-ys) / xs) if relative else (xs - ys)


def main(args):
    stripAccession = args['--strip_accession']

    ##
    # Parameters for reading in the expression files
    ##
    exp1 = ExpressionTools.parseExpressionFile(args['EXP1'], stripAccession=stripAccession)
    name1 = args['--name'][0]
    norm1 = args['--norm'][0]

    exp2 = ExpressionTools.parseExpressionFile(args['EXP2'], stripAccession=stripAccession)
    name2 = args['--name'][1]
    norm2 = args['--norm'][1]

    ExpressionTools.normalizeSets(exp1, exp2)

    ##
    # Keep track of number of points who's values we truncate
    ##
    trunc1, trunc2 = 0, 0
    trunc1NZ, trunc2NZ = 0, 0
    e1 = float(args['--m1'])
    e2 = float(args['--m2'])
    newExp = []

    # Read in the values from the first dataset, truncating the expression
    # level if necessary
    for e in exp1.exps_:
        ev = e.expression
        if e.expression < e1:
            ev = 0.0
            trunc1 += 1
            trunc1NZ += 1 if e.expression > 0.0 else 0
        
        newExp.append(ExpressionTools.ExpressionDatum(e.name, e.length, ev))

    exp1.exps_ = newExp

    # Read in the values from the second dataset, truncating the expression
    # level if necessary
    newExp = []
    for e in exp2.exps_:
        ev = e.expression
        if e.expression < e2:
            ev = 0.0
            trunc2 += 1
            trunc2NZ += 1 if e.expression > 0.0 else 0

        newExp.append(ExpressionTools.ExpressionDatum(e.name, e.length, ev))

    exp2.exps_ = newExp

    ##
    ## normalization and pairing
    ##
    exp1.normalize(norm1)
    exp2.normalize(norm2)

    # zip the results back up into matches
    matches = exp1.zipWithMatching(exp2)

    # get just the expression values as x and y
    x, y = zip(*[(e[0].expression, e[1].expression) for e in matches])

    sx = sorted(exp1.exps_, key = lambda x: x.expression)
    sy = sorted(exp2.exps_, key = lambda x: x.expression)

    ## Some summary statistics about the data
    print("{}: min = {}({}), max = {}({})".format(name1, sx[0].expression, sx[0].name, sx[-1].expression, sx[-1].name))
    print("{}: min = {}({}), max = {}({})".format(name2, sy[0].expression, sy[0].name, sy[-1].expression, sy[-1].name))

    ##
    # If we need to 'align' the datasets
    ##

    #yscale = sum([e.expression for e in exp1.exps_]) / sum([e.expression for e in exp2.exps_])
    #ys = [yscale * e[1].expression for e in matches]

    #ys = [e[1].expression for e in matches]
    #meanX = sum(x) / float(len(x))
    #meanY = sum(ys) / float(len(ys))
    #yshift = meanX - meanY

    ##
    # Don't do any 'alignment' for now
    ##
    yscale = 1.0
    yshift = 0.0

    print("yscale = {}".format(yscale))
    print("yshift = {}".format(yshift))

    matches = [(ExpressionTools.ExpressionDatum(e[0].name, e[0].length, e[0].expression),\
                ExpressionTools.ExpressionDatum(e[1].name, e[1].length, yscale * e[1].expression + yshift ) )\
                for e in matches]
    x, y = zip(*[(e[0].expression, e[1].expression) for e in matches])

    ##
    # Compute the RMSE and median percentage error
    ##
    rmse = ExpressionTools.RMSE(matches)
    medPE = ExpressionTools.medPE(matches)

    print("RMSE {} vs {} is {:0.2f}".format(name1, name2, rmse))
    #print("Trimmed RMSE {} vs {} is {}".format(name1, name2, ExpressionTools.TrimmedRMSE(matches, 0.005)))
    print("MedPE {} vs {} is {:0.2f}".format(name1, name2, medPE))

    percentError, _ = ExpressionTools.PE(matches)
    print("Percentage error computed on {} points".format(len(percentError)))

    ##
    # Plot the histogram of percentage errors
    ##
    #plt.hist(percentError, bins=100, range=(0,100))
    #plt.ylim(0, 3100)
    #plt.show()

    ##
    # We don't use these numbers right now
    ##
    #relativeErrors = ExpressionTools.relativeErrors(matches)
    #print("EF15 = {}".format(ExpressionTools.errorFraction(relativeErrors, 0.15)))
    #print("IsoEM MedPE {}".format(ExpressionTools.medPEIsoEM(relativeErrors)))

    print("Cutoff removed {} points from dataset 1 ({} were nonzero)".format(trunc1, trunc1NZ))
    print("Cutoff removed {} points from dataset 2 ({} were nonzero)".format(trunc2, trunc2NZ))

    ##
    # Compute and print the correlation statistics
    ##
    pr = sp.stats.pearsonr(x, y)[0]
    sr = sp.stats.spearmanr(x, y)[0]

    print("Pearson r = {0}".format(pr))
    print("Spearman r = {0}".format(sr))

    if args['--noplot']:
        return

    font = {'family': 'normal',
            'weight': 'normal',
            'size': 20}

    ## R style!
    style.use('ggplot')

    plt.rc('font', **font)

    nstr = {"id": "", "log": "($\\log_2$)", "rpkm": "(RPKM)", "lrpkm": "($\\log_2$ RPKM)", "frac": "(tfrac)", "fracLengthNorm": "(tfrac)", "rpkm2tpm": "(TPM)", "tpm" : "(TPM)", "tpmlog" : "($\\log_2$ TPM)"}

    xlabel = "{0}{1}".format(name1, nstr[norm1])
    ylabel = "RPKM{0}".format(nstr[norm2])
    #ylabel = "{0}{1}".format(name2, nstr[norm2]) if name2 != "" else ""

    # with open("datapoints.txt", "wb") as ofile:
    #     ofile.write('"{}","{}"\n'.format(xlabel, ylabel))
    #     for xi, yi in itertools.izip(x, y):
    #         ofile.write('{},{}\n'.format(xi, yi))


    #minVal = 0
    #maxVal = 20.0

    minVal = min(min(x), min(y))
    maxVal = max(max(x), max(y))

    #plt.axis([minVal, maxVal, minVal, maxVal])
    plt.axis([min(x), max(x), min(y), max(y)])

    # Don't draw the top and right axes
    for loc, spine in plt.gca().spines.items():
        if loc in ['left','bottom']:
            spine.set_position(('outward',10)) # outward by 10 points
            spine.set_color('black')
        elif loc in ['right','top']:
            spine.set_color('none') # don't draw spine
        else:
            raise ValueError('unknown spine location: %s'%loc)

    # ticks point outward
    plt.gca().tick_params(axis='both', direction='out')

    # remove unneeded ticks
    plt.gca().get_xaxis().tick_bottom()
    plt.gca().get_yaxis().tick_left()


    # plt.text(0.30, 0.90, r"$\sigma = {0:.3},\, \rho = {1:.3}$".format(pr, sr),
    #     fontsize=24,
    #     horizontalalignment='center',
    #     verticalalignment='center',
    #     transform = plt.axes().transAxes)

    plt.xlabel(xlabel)
    if not args['--noylabel']:
        plt.ylabel(ylabel)

    ## Draw the correlation values inside the plot
    #plt.gca().text(min(x), max(y)-2, r"$\sigma$={0:0.2f}, $\rho$={1:0.2f}".format(pr, sr))
    #plt.gca().text(-4, 20.5, r"$\sigma$={0:0.2f}, $\rho$={1:0.2f}, RMSE={2:0.2f}, mPE={3:0.2f}".format(pr, sr, rmse, medPE))
    #plt.gca().text(-2, 19, r"$\sigma$={0:0.2f}, $\rho$={1:0.2f}".format(pr, sr, rmse, medPE))

    #plt.ylabel("{0}{1}".format(name2, nstr[norm2]))

    plt.gca().set_aspect('equal')

    #cb = plt.colorbar()
    if args['--scatter']:

        # Why hexbin here you say?  Well, it's because I don't know how else
        # to get dotted lines (like those of the grid) around the plot.  Hexbin
        # seems to do this automatically. So for the sake of consistency with the
        # hexbin plots below, I want these plots to have the dotted boundary grid
        # lines as well.
        plt.hexbin([0], [0], extent=(minVal, maxVal, minVal, maxVal), mincnt=1.0, alpha=0.0)
        plt.scatter(x, y, alpha=0.7)

    else:
        plt.hexbin(x, y, extent=(minVal, maxVal, minVal, maxVal), mincnt=1.0,
                   gridsize=100, bins='log', cmap=plt.cm.YlOrRd, alpha=0.8)


    plt.gca().xaxis.grid(color='gray', linestyle='dashed')
    plt.gca().yaxis.grid(color='gray', linestyle='dashed')

    # import random
    # for d in matches:
    #     if d[0].expression > 8.0 and d[0].expression < 8.5 and d[1].expression > 2.5 and d[1].expression < 3.0:
    #         plt.annotate(
    #             d[0].name,
    #             xy = (d[0].expression + random.random(), d[1].expression + random.random()), xytext = (-20, 20),
    #             textcoords = 'offset points', ha = 'right', va = 'bottom',
    #             bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
    #             arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

    of = open('strange_blob.txt', 'wb')
    for xi, yi in matches:
        if xi.expression >= 4 and xi.expression <= 9:
            if yi.expression >= 11 and yi.expression <= 15:
                of.write("{}\t{}\t{}\n".format(xi.name, xi.expression, yi.expression))
    of.close()

    ## uncomment this to print the title
    #plt.title(args['--name'][1])

    plt.tight_layout()
    if args['--out']:
        print("saving to {}".format(args['--out']))
        plt.savefig(args['--out'])
    else:
        plt.show()

    # smallestScore = 0.0

    # pcs, prs, sps = [], [], []
    # cpr, csr = [], []
    # relDiffs = []
    # step = 5

    # for percentile in xrange(step, 100 + step, step):
    #     percentile = min(100, percentile)
    #     score = stats.scoreatpercentile(x, percentile)
    #     srange = (smallestScore, score)
    #     smallestScore = score
    #     pr, sr = correlationInRange(x, y, srange)

    #     relDiffs.append(differencesInRange(x, y, srange))

    #     print("Percentile range [{}, {})".format(percentile - step, percentile))
    #     print("Score range [{}, {})".format(srange[0], srange[1]))
    #     print("pearson = {}, spearman = {}".format(pr, sr))
    #     print("=" * 80)

    #     pcs.append(percentile)
    #     prs.append(pr)
    #     sps.append(sr)

    #     pr, sr = correlationInRange(x, y, (0.0, srange[1]))
    #     cpr.append(pr)
    #     csr.append(sr)

    # plt.scatter(pcs, prs)
    # plt.title('Pearson correlations at percentile ({})'.format(name2))
    # plt.xlim(0, 101)
    # plt.show()

    # relErrors = [(0, 0.0)] * len(x)
    # for i in xrange(len(x)):
    #     xi, yi = x[i], y[i]
    #     relError = abs(x[i]-y[i]) / abs(x[i])
    #     relErrors[i] = (i, relError)

    # relErrors = sorted(relErrors, key=lambda x: x[1], reverse=True)
    # idxs = [i for i,e in relErrors]
    # xp = np.array([x[i] for i in idxs])
    # yp = np.array([y[i] for i in idxs])
    # truncPrs = [pr]
    # truncSrs = [sr]
    # for i in xrange(0, len(xp)):
    #     prp = sp.stats.pearsonr(xp[i:], yp[i:])[0]
    #     srp = sp.stats.spearmanr(xp[i:], yp[i:])[0]

    #     truncPrs.append(prp)
    #     truncSrs.append(srp)

    # plt.plot(truncPrs, label='Pearson Correlations')
    # plt.plot(truncSrs, label='Spearman Correlations')
    # plt.axhline(y=0.931)
    # plt.axhline(y=0.941)
    # plt.legend(loc=2)
    # plt.title("Variation in correlation")
    # plt.xlabel("# of datapoints removed (most violated first)")
    # plt.ylabel("correlation")
    # plt.show()

if __name__ == "__main__":
    args = docopt(__doc__, version="PlotCorrelation 1.0")
    main(args)
