#!/usr/bin/env python

import glob
import re

def combine_per_panel():
    file_list=glob.glob('pool_*_panel_*_gtcheck_donor_assignments.csv')
    with open ("combined_gt_donor_assignment_per_panel.csv", "w") as fout:
        fout.write('pool,panel,donor_query,donor_gt,donor_1,score0,score1,score_n,n,mean,sd,z0,z1\n')
        pattern = re.compile(r"^pool_(.+?)_panel_(.+?)_gtcheck_donor_assignments\.csv$")
        for file_name in file_list:
            m = pattern.search(file_name)
            pool_id = m.group(1)
            panel_id = m.group(2)
            with open (file_name, "r") as fin:
                next(fin)
                for line in fin:
                    fout.write(pool_id+","+panel_id+","+line)

def combine_overall():
    file_list=glob.glob('stats_*_gt_donor_assignments.csv')
    with open ("combined_stats_gt_donor_assignments.csv", "w") as fout:
        fout.write('pool,donor_query,donor_gt,score0,score1,score_n,n,mean,sd,z0,z1,final_panel\n')
        pattern = re.compile(r"^stats_(.+?)_gt_donor_assignments\.csv$")
        for file_name in file_list:
            print(file_name)
            m = pattern.search(file_name)
            pool_id = m.group(1)
            print(pool_id)
            with open (file_name, "r") as fin:
                next(fin)
                for line in fin:
                    fout.write(pool_id+","+line)
def main ():
    combine_per_panel()
    if glob.glob('stats_*_gt_donor_assignments.csv'):
        combine_overall()

if __name__ == '__main__':
    main()
