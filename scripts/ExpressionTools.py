from collections import namedtuple
from collections import defaultdict
import numpy as np

ExpressionDatum = namedtuple('ExpressionDatum', 'name length expression')

class ExpresionDataset(object):

  def __init__(self):
    self.exps_ = []

  def addElement(self, datum):
    """
    Add the element datum to this ExpresionDataset
    """
    self.exps_.append(datum)

  def aggregateWith(self, edict):
    """
    If we have a set of different observations for a datapoins (e.g. isoforms of a gene), we aggregate
    them together using this function
    """
    aexp = ExpresionDataset()
    adict = defaultdict(list)
    missing = 0

    haveNone = self.exps_[0].length is None
    for d in self.exps_:
        if d.name in edict:
            k = edict[d.name]
            adict[k].append(d)
        else: 
            missing += 1

    if haveNone:
        aggregatedExps = [ExpressionDatum(k, None, sum([v.expression for v in e])) for k, e in adict.iteritems()]
    else:
        aggregatedExps = [ExpressionDatum(k, sum([v.length for v in e]), sum([v.expression for v in e])) for k, e in adict.iteritems()]

    aexp.exps_ = aggregatedExps

    print("The expression dictionary was missing a total of {} keys".format(missing))
    return aexp

  def normalize(self, method="id"):
    """
    How the datapoints should be "normalized".  The relevant options (currently) are 
    id:   this is a passthrough, the expression values remain unchanged
    log:  those datapoints with a value > 0 will have their expression logged (the others are discarded)
    frac: the expression value of each datum is set to the fraction of total expression in the dataset it accounts for
    tpm:  the expressions are normalized to represent transcript fractions, and then multiplied by 10^6 to represent TPM
    """
    billion = np.power(10,9)
    million = np.power(10,6)        

    if method == "id":
        pass
    elif method == "log":
        slen = len(self.exps_)
        self.exps_ = filter(lambda x: x.expression > 0, self.exps_)
        print("Removed {0} of {1} original expression values for log transform".format(slen-len(self.exps_), slen))
        print("Minimum expression value is {0}".format(min([e.expression for e in self.exps_])))
        self.exps_ = [ExpressionDatum(e.name, e.length, np.log2(e.expression)) for e in self.exps_]
        print("Minimum logged expression value is {0}".format(min([e.expression for e in self.exps_])))
    elif method == "frac":
        tot = np.sum(d.expression for d in self.exps_)
        print("tot = {}".format(tot))
        self.exps_ = [ExpressionDatum(e.name, e.length, e.expression / float(tot)) for e in self.exps_]
    elif method == "tpm":
        self.normalize("frac")
        self.exps_ = [ExpressionDatum(e.name, e.length, e.expression * million) for e in self.exps_]
    # elif method == "tpmlog":
    #     self.normalize("tpm")
    #     self.exps_ = filter(lambda x: x.expression > 0, self.exps_)
    #     self.exps_ = [ExpressionDatum(e.name, e.length, np.log2(e.expression)) for e in self.exps_ if e.expression > 0.0]
    # elif method == "rpkm2tpm":
    #     rpkms = np.array([e.expression for e in self.exps_])
    #     tpms = million * (rpkms / rpkms.sum())
    #     self.exps_ = [ExpressionDatum(e.name, e.length, tpms[i]) for i,e in enumerate(self.exps_)]
    # elif method == "rpkmFrac":
    #     self.normalize("rpkm")
    #     self.normalize("frac")
    # elif method == "fracDivLength":
    #     tot = np.sum(d.expression / d.length for d in self.exps_)
    #     self.exps_ = [ExpressionDatum(e.name, e.length, (e.expression / e.length) / float(tot)) for e in self.exps_]
    # elif method == "fracLengthNorm":
    #     for d in self.exps_:
    #         if d.length == 0.0:
    #             if(d.expression != 0.0):
    #                 print(d)
    #     tot = np.sum(d.expression * d.length for d in self.exps_)
    #     print("tot = {}".format(tot))
    #     self.exps_ = [ExpressionDatum(e.name, e.length, e.expression * e.length / float(tot)) for e in self.exps_]
    else:
        raise ValueError("The requested normalization method ({}) was not an acceptable type".format(method))

  def zipWithMatching(self, other):
    """
    This function returns a list of tuples, where each tuple contains
    an ExpressionDatum from this datset and an ExpressionDatum from other, where the elements 
    of each tuple share the same name.  Any data present exclusively in either dataset (i.e. those
    whose names can't be matched) will not be present in the returned list.
    """
    self.exps_ = sorted(self.exps_, key=lambda x: x.name)
    other.exps_ = sorted(other.exps_, key=lambda x: x.name)

    sind, oind = 0, 0
    slen, olen = len(self.exps_), len(other.exps_)
    matches = []
    
    droppedSelf, droppedOther = [], []

    while sind < slen and oind < olen:
        sdat = self.exps_[sind]
        odat = other.exps_[oind]
        if sdat.name == odat.name:
            matches.append((sdat, odat))
            oind += 1
            sind += 1
        elif sdat.name < odat.name:
            droppedSelf.append(ExpressionDatum(sdat.name, sdat.length, sdat.expression))
            sind += 1
        else:
            droppedOther.append(ExpressionDatum(odat.name, odat.length, odat.expression))
            oind += 1

    #print("discarded {} samples from the ground truth".format(len(droppedSelf)))
    #print("discarded {} samples from the estimates truth".format(len(droppedOther)))
    doPlot = False
    if doPlot:
      from matplotlib import pyplot as plt
      if len(droppedOther) > 0:
        plt.hist([e.expression for e in droppedOther], range=(-20, 20.0), bins=50, color=[(85/255.0, 212/255.0, 210/255.0, 50/255.0)], label=['dropped from estimated'])
      if len(droppedSelf) > 0:
        plt.hist([e.expression for e in droppedSelf], range=(-20, 20.0), bins=50, color=[(212/255.0, 85/255.0, 87/255.0, 50/255.0)], label=['dropped from ground truth'])

      if len(droppedOther) > 0 or len(droppedSelf) > 0:
        plt.xlabel("$\log_2$(RPKM)")
        plt.legend()
        plt.show()
 
    return matches


def normalizeSets(s1, s2):
    """
    This function takes two expression datasets, s1 and s2, and
    "normalizes" them to be on the same set of transcripts.  Any transcript
    present in s1 but not s2 will be inserted into s2 with an expression of 0 and
    vice versa.
    """

    expDictSelf = {x.name : x for x in s1.exps_}
    expDictOther = {x.name : x for x in s2.exps_}

    addToSelf = set(expDictOther.keys()) - set(expDictSelf.keys())
    addToOther = set(expDictSelf.keys()) - set(expDictOther.keys())

    for k in addToSelf:
      e = expDictOther[k]
      s1.addElement(ExpressionDatum(e.name, e.length, 0.0))

    for k in addToOther:
      e = expDictSelf[k]
      s2.addElement(ExpressionDatum(e.name, e.length, 0.0))


import numpy as np
def RMSE(zippedList):
    """
    Given a list of paired datapoints (e.g. as produced by zipWithMatching), compute the Root Mean Squared Error
    (RMSE) of the expression values in this list.
    """
    return np.sqrt(np.array([(x.expression-y.expression)**2 for x,y in zippedList ]).sum() / len(zippedList))

import scipy as sp
from scipy import stats
def TrimmedRMSE(zippedList, proportionCut):
    """
    Compute the RMSE where the top "propotionCut" percent of the errors are discarded (to account for outliers)
    """
    squaredErrors = sorted(np.array([(x.expression-y.expression)**2 for x,y in zippedList ]))
    squaredErrors = sp.stats.trim1(squaredErrors, proportionCut, 'right')
    return np.sqrt(squaredErrors.sum() / len(squaredErrors))

def PE(zippedList, truncVal=0.0):
    """
    This function returns two lists; percentageErrors and aug.
    The percentageErrors list is a list of the percentage error for each pair of points in zipped list
    where the percentage error is defined as: 
        PE(x,y) = 100.0 * |x-y| / x
    When x is 0, PE(x,y) is undefined and such points are not part of the returned list

    The aug list is simply the input zippedList, where each 2-tuple is replaced by a 3-tuple where the
    3rd element is the computed percentage error.  As mentioned above, any pair of points with an undefined
    percentage error are omitted from the list.
    """
    relativeErrors = np.ones(len(zippedList)) * np.nan
    for i, (x, y) in enumerate(zippedList):
        if x.expression > 0.0:
            relativeErrors[i] = abs(x.expression - y.expression) / x.expression

    percentageErrors = 100.0 * relativeErrors
    
    aug = [ (x, y, pe) for (x, y), pe in zip(zippedList, percentageErrors) if not np.isnan(pe) ]
    percentageErrors = percentageErrors[ np.isfinite(percentageErrors) ]
    percentageErrors = np.sort(percentageErrors)
    return percentageErrors, aug

def relativeErrors(zippedList):

    relativeErrors = []
    # million = np.power(10,6)
    # billion = np.power(10,6)    
    # xexp, _ = zip(*zippedList)
    # xs = np.array([x.expression for x in xexp if x.expression > 0.0])
    # c = sp.stats.scoreatpercentile(xs, 10)
    # npmDivLen = [ (x.expression) / billion for x,y in zippedList]
    # norm = sum(npmDivLen)
    # tpm = [million * npm / norm for npm in npmDivLen]

    for i, (x, y) in enumerate(zippedList):
        if abs(x.expression) > 0.0:
            relativeErrors.append( abs(x.expression - y.expression) / x.expression )
        elif x.expression == y.expression == 0.0:
            relativeErrors.append(0.0)
        elif x.expression == 0.0 and y.expression > x.expression:
            relativeErrors.append(np.inf)

    relativeErrors = np.sort(relativeErrors)
    return relativeErrors

def errorFraction(relativeErrors, t):
    return len( relativeErrors[relativeErrors >= t] ) / float(relativeErrors.shape[0])

def medPEIsoEM(relativeErrors):
    t = -1.0
    for ret in relativeErrors:
        if ret <= t:
           continue 
        eft = errorFraction(relativeErrors, ret)
        #print(ret, eft)
        if eft <= 0.5:
            return ret
        t = ret

def efracList(relativeErrors):
    t = -1.0
    thresh, efs = [], []
    for ret in relativeErrors:
        if ret <= t:
           continue 
        if ret >= 1.0:
            break
        eft = errorFraction(relativeErrors, ret)
        #print(ret, eft)
        thresh.append(ret)
        efs.append(eft)
    return thresh, efs

def medPE(zippedList):
    """
    Compute the median percentage error of the paired data in zippedList.
    """
    perror, auglist = PE(zippedList)
    # sl = sorted(auglist, key=lambda x: x[2])
    # sl = sl[int(0.99*len(sl)):-1]
    # print('\n'.join([ "Transcript: {}, ground truth = {}, predicted = {}, relative error = {}".format(
    #                 x.name, x.expression, y.expression, pe) for x, y, pe in sl]))
    return np.median(perror)


def truncMedPE(zippedList, proportionCut):
    predErrors = PE(zippedList, proportionCut)
    return np.median(predErrors)

def accession(tname):
    return tname.split('.')[0]

def version(tname):
    return int(tname.split('.')[-1])

def parseExpressionFile(fname, epsilon=0.0, stripAccession=False):
    """
    Parse the file given by fname into an ExpressionDataset object.  All expression values
    (values in the 3rd column) will be augmented by epsilion (default value of 0).
    
    If stripAccession is set to True, then the accession is stripped of the variant information,
    assumed to be any part of the string after the first '.'
    """
    dset = ExpresionDataset()
    
    ifile = None
    if isinstance(fname, str):
      ifile = open(fname,'rb')
    else:
      ifile = fname

    for l in ifile:
      toks = l.rstrip().split()
      ntoks = len(toks)
      hasLength = ntoks == 3
      
      dname = accession(toks[0]) if stripAccession else toks[0]
      dlength = float(toks[1]) if hasLength else None
      dexp = float(toks[-1]) if toks[-1] != "NA" else 0.0
      dset.addElement(ExpressionDatum(dname, dlength, dexp + epsilon))

    ifile.close()

    return dset
