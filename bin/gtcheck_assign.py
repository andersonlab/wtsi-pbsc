#!/usr/bin/env python

## pick out the closes matching pair from output by bcftools gtcheck

DEBUG = False

import sys
import statistics
import copy

SMVAL = 1.0E-29

class ScoreMatrix:
    def __init__(self):
        self.is_sum = False
        self.ref_ids = set()
        self.gtchkd = {}

    def parse_gtcheck_scores(self, infh):
        ctr = 0
        for lin in infh:
            if lin[:2] == 'DC' or lin[:2] == 'DCv2':
                fld = lin.strip().split('\t')
                if DEBUG:
                    print(fld)
                if len(fld) > 3 and 'DC' in fld[0].strip():
                    ctr += 1
                    donor_query = fld[1].strip()
                    donor_ref = fld[2].strip()
                    score = float(fld[3])
                    tup = (donor_ref, score)
                    try:
                        self.gtchkd[donor_query].append(tup)
                    except KeyError:
                        self.gtchkd[donor_query] = [tup]
        return ctr

    def sum_scores(self):
        for d in tuple(self.gtchkd.keys()):
            scoresd = {}
            for ref_id, score in self.gtchkd[d]:
                try:
                    scoresd[ref_id] += score
                except KeyError:
                    scoresd[ref_id] = score
            self.gtchkd[d] = list(scoresd.items())
            self.ref_ids = self.ref_ids.union(set(scoresd.keys()))
        self.is_sum = True
        return

    def write_score_matrix(self, oufh):
        # get names of reference panel
        mtx = []
        donor_ids = tuple(self.gtchkd.keys())
        # print(gtchkd.values())
        genotype_ids = tuple(self.ref_ids)
        n_donors = len(donor_ids)
        n_genotypes = len(genotype_ids)
        row_init = [float(0)]*n_genotypes

        for d in tuple(gtchkd.keys()):
            donor_id = donor_ids[i]
            mtx.append([0 for k in range(n_genotypes)])
            for (gtid, score) in gtchkd[donor_id]:
                try:
                    j = genotype_ids.index(gtid)
                except ValueError:
                    sys.exit("ERROR: genotype id '{:s}' not found.".format(gtid))
                mtx[i][j] = score
        for d in donor_ids:
            oufh.write(",{:s}".format(d))
        oufh.write("\n")
        for j in range(n_genotypes):
            oufh.write("{:s}".format(genotype_ids[j]))
            for i in range(n_donors):
                oufh.write(",{:.1f}".format(mtx[i][j]))
            oufh.write("\n")
        return mtx, donor_ids, genotype_ids

    def write_score_array(self, oufh):
        oufh.write("gtcheck_score,donor_id\n")
        for d in self.gtchkd:
            for (gtid, score) in self.gtchkd[d]:
                oufh.write("{:s},{:.1f},{:s}\n".format(d, score, gtid))
        return

    def find_matches(self):
        rsltd = {}
        for k in self.gtchkd:
            print(k)
            self.gtchkd[k].sort(key = lambda a: a[1])
            if DEBUG:
                print(self.gtchkd[k])
            donor0, score0 = self.gtchkd[k][0]
            donor1, score1 = None, 0
            score_n = 0
            n = 0
            d0 = float(0)
            d1 = float(0)
            sd = float(0)
            m = float(0)
            sd = SMVAL
            if len(self.gtchkd[k]) > 1:
                donor1, score1 = self.gtchkd[k][1]
                dd1 = score1 - score0
                scores = []
                for i in range(1, len(self.gtchkd[k])):
                    scores.append(self.gtchkd[k][i][1])
                m = statistics.mean(scores)
                d0 = m - score0
                d1 = m - score1
                score_n = scores[-1]
                n = len(scores)
                if n > 1:
                    sd = statistics.stdev(scores)
                    if sd < SMVAL:
                        sd = SMVAL
                    if d0 < 1.5*sd:
                        print("WARNING: ambiguous assignment", donor0, donor1, score0, score1, m, sd, d0, d1)
                        sys.stderr.write("ERROR: {:s} assignment ambiguous.\n".format(donor0))
                else:
                    if(score1==score0):
                        sd = 0.0000000001
                    else:
                        sd = score1 - score0
            z0 = d0/sd
            z1 = d1/sd
            rsltd[k] = (donor0, donor1, score0, score1, score_n, n, m, sd, z0, z1, dd1/sd)
        return rsltd

if __name__ == '__main__':
    narg = len(sys.argv)
    if narg < 3:
        sys.exit("usage: {:s} <ouput file prefix> <input file> [<input file 2>] ...".format(sys.argv[0]))

    oufnprfx = sys.argv[1]
    fnin_arr = tuple(sys.argv[2:])

    fnmtx = oufnprfx + '_gtcheck_score_table.csv'
    fnout = oufnprfx + '_gtcheck_donor_assignments.csv'

    scm = ScoreMatrix()

    if len(fnin_arr) < 2 and \
        fnin_arr[0] == '-':
        infh = sys.stdin
        scm.parse_gtcheck_scores(infh)
    else:
        for fn in fnin_arr:
            infh = open(fn, 'r')
            nlin = scm.parse_gtcheck_scores(infh)
            sys.stderr.write("# {:d} lines read from file '{:s}'\n".format(nlin, fn))
            infh.close()
        sys.stderr.write("Summing scores over {:d} output files ...\n".format(len(fnin_arr)))
        scm.sum_scores()
    oufh_mtx = open(fnmtx, 'w')
    #write_score_matrix(oufh_mtx, gtchkd)
    scm.write_score_array(oufh_mtx)
    oufh_mtx.close()

    oufh = open(fnout, 'w')
    oufh.write("donor_query,donor_gt,donor_1,score0,score1,score_n,n,mean,sd,z0,z1\n")

    rsltd = scm.find_matches()
    if DEBUG:
        print(rsltd)
    sys.stdout.write("\nAssignments:\n")
    for k in rsltd:
        (donorid, donor1_id, score0, score1, score_n, n, m, sd, z0, z1, zz1) = rsltd[k]
        sys.stdout.write(
            "{:s} -> {:s} (Best), {:s} (Second-best), score0 = {:.1f}, score1 = {:.1f}, score_last = {:.1f}, n_samples = {:d}, "
            "mean = {:.1f}, sd = {:.1f}, z0 = {:.1f}, z1 = {:.1f}, zz1 = {:.1f}\n"
            .format(k, donorid, donor1_id, score0, score1, score_n, n, m, sd, z0, z1, zz1)
        )
        oufh.write("{:s},{:s},{:s},{:.1f},{:.1f},{:.1f},{:d},{:.1f},{:.1f},{:.1f},{:.1f}\n"
            .format(k, donorid, donor1_id, score0, score1, score_n, n, m, sd, z0, z1)
        )
    oufh.close()

    sys.exit(0)
