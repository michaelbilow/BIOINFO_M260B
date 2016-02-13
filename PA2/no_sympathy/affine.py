"""This code was taken from
https://www.cs.umd.edu/class/fall2010/cmsc423/align423.py"""

Infinity = float('inf')


def make_matrix(sizex, sizey):
    """Creates a sizex by sizey matrix filled with zeros."""
    return [[0]*sizey for i in xrange(sizex)]


def print_matrix(x, y, A):
    """Print the matrix with the (0,0) entry in the top left
    corner. Will label the rows by the sequence and add in the
    0-row if appropriate."""

    # decide whether there is a 0th row/column
    if len(x) == len(A):
        print "%5s" % (" "),
    else:
        print "%5s %5s" % (" ","*"),
        y = "*" + y

    # print the top row
    for c in x:
        print "%5s" % (c),
    print

    for j in xrange(len(A[0])):
        print "%5s" % (y[j]),
        for i in xrange(len(A)):
            print "%5.0f" % (A[i][j]),
        print


def is_complement(a, b):
    """Return True if character a is complmentary to character b"""
    assert len(a) == len(b) == 1
    return (a.upper(), b.upper()) in [
        ("A", "T"), ("T", "A"),
        ("C", "G"), ("G", "C"),
        ("A", "U"), ("U", "A")
    ]


#=============================================================
# Alignment Parameters
#=============================================================

class ScoreParam:
    """Stores the parameters for an alignment scoring function"""
    def __init__(self, match, mismatch, gap, gap_start=0):
        self.gap_start = gap_start
        self.gap = gap
        self.match = match
        self.mismatch = mismatch

    def matchchar(self, a,b):
        """Return the score for aligning character a with b"""
        assert len(a) == len(b) == 1
        if a==b:
            return self.match
        else:
            return self.mismatch

    def __str__(self):
        return "match = %d; mismatch = %d; gap_start = %d; gap_extend = %d" % (
                self.match, self.mismatch, self.gap_start, self.gap
        )

#=============================================================
# Sequence Alignment
#=============================================================

def local_align(x, y, score=ScoreParam(10, -5, -7)):
    """Do a local alignment between x and y with the given scoring parameters.
    We assume we are MAXIMIZING."""

    # create a zero-filled matrix
    A = make_matrix(len(x) + 1, len(y) + 1)

    best = 0
    optloc = (0,0)

    # fill in A in the right order
    for i in xrange(1, len(x)+1):
        for j in xrange(1, len(y)+1):

            # the local alignment recurrance rule:
            A[i][j] = max(
               A[i][j-1] + score.gap,
               A[i-1][j] + score.gap,
               A[i-1][j-1] + score.matchchar(x[i-1], y[j-1]),
               0
            )

            # track the cell with the largest score
            if A[i][j] >= best:
                best = A[i][j]
                optloc = (i,j)

    print "Scoring:", str(score)
    print "A matrix ="
    print_matrix(x, y, A)
    print "Optimal Score =", best
    print "Max location in matrix =", optloc
    # return the opt score and the best location
    return best, optloc


def affine_align(x, y, score=ScoreParam(10, -2, -7, -15)):
    """Global alignment with affine penalties. We assume we are maximizing."""
    M = make_matrix(len(x) + 1, len(y) + 1)
    X = make_matrix(len(x) + 1, len(y) + 1)
    Y = make_matrix(len(x) + 1, len(y) + 1)

    for i in xrange(1, len(x)+1):
        M[i][0] = -Infinity
        X[i][0] = -Infinity
        Y[i][0] = score.gap_start + i * score.gap

    for i in xrange(1, len(y)+1):
        M[0][i] = -Infinity
        X[0][i] = score.gap_start + i * score.gap
        Y[0][i] = -Infinity

    for i in xrange(1, len(x)+1):
        for j in xrange(1, len(y)+1):

            M[i][j] = score.matchchar(x[i-1], y[j-1]) + max(
                    M[i-1][j-1],
                    X[i-1][j-1],
                    Y[i-1][j-1]
            )

            X[i][j] = max(
                    score.gap_start + score.gap + M[i][j-1],
                    score.gap + X[i][j-1],
                    score.gap_start + score.gap + Y[i][j-1]
            )

            Y[i][j] = max(
                    score.gap_start + score.gap + M[i-1][j],
                    score.gap_start + score.gap + X[i-1][j],
                    score.gap + Y[i-1][j]
            )

    opt = max(M[len(x)][len(y)], X[len(x)][len(y)], Y[len(x)][len(y)])

    print "x = %s & y = %s" % (x,y)
    print "Scoring:", str(score)
    print "M matrix ="
    print_matrix(x,y,M)
    print "X matrix ="
    print_matrix(x,y,X)
    print "Y matrix ="
    print_matrix(x,y,Y)
    print "Optimal =", opt

    return opt


if __name__ == "__main__":
    affine_align("GTACC", "GCC", ScoreParam(gap_start=-6, gap=-1, match=3, mismatch=-2))
