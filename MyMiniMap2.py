import sys
import getopt
import math
from bio import *
#https://dgenies.toulouse.inra.fr/documentation/formats#:~:text=PAF%20is%20the%20default%20output%20format%20of%20minimap2.&text=If%20PAF%20is%20generated%20from,the%20BLAST%2Dlike%20alignment%20identity.
def collectMinimizers(s, k:Static[int], w:Static[int]):
    minimizers:set[tuple[int, int, int]] = set()
    for i in range(len(s) - w - k + 1):
        min_kmers = [(int(Kmer[k](s[i:i + k])), i, 0)]
        for j, kmer in kmers_with_pos(s[i: i + k + w -1], k=k, step=1):
            if int(kmer) < min_kmers[0][0]:
                min_kmers = [(int(kmer), i + j, 0)]
            elif int(kmer) == min_kmers[0][0]:
                min_kmers.append((int(kmer), i + j, 0))
            if int(~kmer) < min_kmers[0][0]:
                min_kmers = [(int(~kmer), i + j, 1)]
            elif int(~kmer) == min_kmers[0][0]:
                min_kmers.append((int(~kmer), i + j, 1))
        for x in min_kmers:
            minimizers.add(x)
    return minimizers

         
def stripAmbiguousBases(s:seq) -> seq:
    t:str = str(s).upper()
    t = ''.join([x for x in t if x in ['A', 'C', 'G', 'T']])
    return seq(t)

def bwt(s:str):
    alphabet = '$ACTG'
    rc_alphabet = '$TGCA'
    shifts = []
    for i in range(len(s)):
        shifts.append(s[i:] + s[:i])
    shifts.sort()
    Occ = [{c: 0 for c in alphabet} for i in range(len(shifts) + 1)]
    for i, s in enumerate(shifts):
        Occ[i + 1] = {c: j for c, j in Occ[i].items()}
        Occ[i + 1][s[-1]] += 1
    C = {c: sum([1 for s in shifts if s[0] < c]) for c in alphabet}
    return shifts, C, Occ

class Chain:
    QSeqName: str
    QSeqLength: int
    Qstart: int
    Qend: int
    R: str
    RefSeqName: str
    RefSeqLength: int
    RefStart: int
    RefEnd: int
    Matches: int
    Bases: int
    qualityScore:float
    def __init__(self, PrimaryChain:tuple[list[tuple[float, tuple[int, int, int, int]]], list[float]], g, k, seqName, seqLength, qScore, refNames):
        self.QSeqName= seqName
        self.QSeqLength = seqLength
        self.Qstart = PrimaryChain[0][-1][1][3]
        self.Qend = PrimaryChain[0][0][1][3] + k
        self.R = '+' if PrimaryChain[0][0][1][1] == 0 else '-'
        self.RefSeqName = refNames[PrimaryChain[0][0][1][0]]
        self.RefStart = PrimaryChain[0][-1][1][2]
        self.RefEnd = PrimaryChain[0][0][1][2] + k
        self.Bases = self.Qend - self.Qstart
        self.Matches = self.Bases - g
        self.qualityScore = qScore
        return
    
    def __str__(self):
        return f"{self.QSeqName}\t{self.QSeqLength}\t{self.Qstart}\t{self.Qend}\t{self.R}\t{self.RefSeqName}\t{self.RefStart}\t{self.RefEnd}\t{self.Matches}\t{self.Bases}\t{self.qualityScore}"
    
    def __repr__(self):
        return str(self)

class Aligner:
    minimizers: dict[int, set[tuple[int, int, int]]]
    sequences: list[seq]
    names: list[str]
    #k:Static[int]
    #w:Static[int]
    k:int
    w:int
    G:int
    def __init__(self, refFilename, k:Static[int], w:Static[int]):
        self.minimizers = dict()
        self.sequences = []
        self.names = []
        self.k = k
        self.w = w
        self.G = 50
        t:int = 0
        for read in FASTA(refFilename, fai=False):
            sequence = stripAmbiguousBases(read.seq)
            self.sequences.append(sequence)
            self.names.append(read.name)
            seqMinimizers = collectMinimizers(sequence, k, w)
            for kmer, i, r in seqMinimizers:
                if kmer not in self.minimizers:
                    self.minimizers.update({kmer: set[tuple[int, int, int]]()})
                self.minimizers[kmer].add((t, i, r))
            t += 1

    def map(self, query:seq, name, epsilon:int = 500, k:Static[int] = 14, w:Static[int] = 20):
        def gaps(q:tuple[list[tuple[float, tuple[int, int, int, int]]], list[float]]):
            total = 0
            cur = g = q[0][-1][1][3]
            for i in range(len(q[0]) - 1, 0, -1):
                total += max((q[0][i][1][3] - k - cur, 0))
                cur = q[0][i][1][3]
            return total
        def mapQuality(fp:float, fs:float, m:int):
            return 40*(1-(fs/fp)) * min((1.0, m/10))*math.log(fp)
        
        def a(anc, i, j):
            return min((anc[i][3] - anc[j][3], anc[i][2] - anc[j][2], self.k))
        
        def b(anc, i, j):
            if anc[j][3] >= anc[i][3]:
                return math.inf
            if max((anc[i][3] - anc[j][3], anc[i][2] - anc[j][2])) > self.G:
                return math.inf
            
            l = (anc[i][3] - anc[j][3]) - (anc[i][2] - anc[j][2])

            if l == 0:
                return 0
            else:
                return (0.01*self.k* abs(l)) + (0.5 * math.log2(abs(l)))

        def chainingScores(anchors:list):
            scores:dict[int, tuple[float, int, int]] = dict()
            anchors.sort()
            for i in range(len(anchors)):

                maxF:tuple[float, int, int, int] = (float(self.k), i, -1, 0)
                for j in range(i-1, max(i-1-51, -1), -1):
                    if anchors[i][0] != anchors[j][0] or anchors[i][1] != anchors[j][1]:
                        break
                    f = scores[j][0] + a(anchors, i, j) - b(anchors, i, j)
                    if f > maxF[0]:
                        maxF = (float(f), i, j, 0)   

                scores.update({maxF[1]: (maxF[0], maxF[2], maxF[3])})
            return scores
        
        queryMinimizers = collectMinimizers(query, k, w)
        anchors:list[tuple[int, int, int, int]] = []
        for kmer, i, r in queryMinimizers:
            if kmer not in self.minimizers:
                continue
            for targetSeq, targetIndex, targetR in self.minimizers[kmer]:
                if r == targetR:
                    #anchors.append((targetSeq, 0, i - targetIndex, targetIndex))
                    anchors.append((targetSeq, 0, targetIndex, i))
                else:
                    #anchors.append((targetSeq, 1, i + targetIndex, targetIndex))
                    anchors.append((targetSeq, 1, targetIndex, i))
        scores = chainingScores(anchors)
        chains = []
        scoresKeys = list(scores.keys())
        scoresKeys.sort(reverse=True)
        for i in scoresKeys:
            cur = i
            chain = []
            while True:
                if scores[cur][2] == 1:
                    break
                scores[cur] = (scores[cur][0], scores[cur][1], 1)
                chain.append((scores[cur][0], anchors[cur]))
                if scores[cur][1] == -1:
                    break
                cur = scores[cur][1]
            if len(chain) > 0:
                chains.append(chain)
        
        chains.sort(key=lambda x: x[0], reverse=True)
        Q:list[tuple[list[tuple[float, tuple[int, int, int, int]]], list[float]]] = []
        for c in chains:
            primary = True
            cLen = c[0][1][3] - c[-1][1][3] + k #Position on query sequence of last anchor minus position of first anchor

            for q, s in Q:
                qLen = q[0][1][3] - c[-1][1][3] + k
                overlap = min((c[0][1][3], q[0][1][3])) + k - max((c[-1][1][3], q[-1][1][3]))
                if overlap/min((cLen, qLen)) >= 0.5:
                    primary = False
                    s.append(c[0][0])

            if primary:
                Q.append((c, [0.0]))
        
        PrimaryChains = [Chain(q, gaps(q), self.k, name, len(query), mapQuality(q[0][0][0], max(q[1]), len(q[0][0][1])), self.names) for q in Q]
        
        for p in PrimaryChains:
            print(p)

        '''
        anchors.sort()
        colinearSubsets:list[list[tuple[int, int, int, int]]] = []
        b:int = 0
        
        for e in range(len(anchors)):
            if (e == len(anchors) - 1)\
                or (anchors[e + 1][0] != anchors[e][0])\
                    or (anchors[e + 1][1] != anchors[e][1])\
                        or ((anchors[e + 1][2] - anchors[e][2]) >= epsilon):
                colinearSubsets.append(anchors[b:e + 1])
                b = e + 1
        for x in colinearSubsets:
            print(x)
        '''

def main(argv):
    #Get optional arguments from command line, opts are the flags, args are the values
    opts, args = getopt.getopt(argv[1:], "x:n:m:k:w:r:c")
    if len(args) < 2:
        print("Usage: minimap2.py [options] <ref.fa>|<ref.mmi> <query.fq>")
        print("Options:")
        print("  -x STR      preset: sr, map-pb, map-ont, asm5, asm10 or splice")
        print("  -n INT      mininum number of minimizers")
        print("  -m INT      mininum chaining score")
        print("  -k INT      k-mer length")
        print("  -w INT      minimizer window length")
        print("  -r INT      band width")
        print("  -c          output the cs tag")
        sys.exit(1)

    #Set default values for argumnets here
    preset = None
    miin_chain_score = None
    min_cnt = None
    min_sc = None
    bw = None
    out_cs = False
    k:Static[int] = 14
    w:Static[int] = 5

    for opt, arg in opts:
        if opt == '-x': preset = arg
        elif opt == '-n': min_cnt = int(arg)
        elif opt == '-m': min_chain_score = int(arg)
        elif opt == '-r': bw = int(arg)
        elif opt == '-c': out_cs = True
    a = Aligner(args[0], k = k, w = w)
    print("QSeqName\tQSeqLength\tQstart\tQend\tR\tRefSeqName\tRefStart\tRefEnd\tMatches\tBases\tQuality")
    for read in FASTQ(args[1]):
        a.map(read.seq, read.name, k= k, w=w)
    '''
    if not a: raise Exception("ERROR: failed to load/build index file '{}'".format(args[0]))
    for name, seq, qual in mp.fastx_read(args[1]): # read one sequence
          for h in a.map(seq, cs=out_cs): # traverse hits
              print('{}\t{}\t{}'.format(name, len(seq), h))
    '''
    return


if __name__== '__main__':
   main(sys.argv)