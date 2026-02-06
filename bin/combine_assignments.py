#!/usr/bin/env python

import glob
import re
import sys

def combine_per_panel(not_assigned, z_thresh, z_dist_thresh):
    failed_z_dist=""
    file_list=glob.glob('pool_*_panel_*_gtcheck_donor_assignments.csv')
    with open ("combined_gt_donor_assignments_per_panel.csv", "w") as fout:
        fout.write('pool,donor_query,donor_gt,donor_1,score0,score1,score_n,n,mean,sd,z0,z1,panel\n')
        pattern = re.compile(r"^pool_(.+?)_panel_(.+?)_gtcheck_donor_assignments\.csv$")
        for file_name in file_list:
            m = pattern.search(file_name)
            pool_id = m.group(1)
            panel_id = m.group(2)
            with open (file_name, "r") as fin:
                next(fin)
                for line in fin:
                    fout.write(pool_id+","+line.strip()+","+panel_id+"\n")
                    key=pool_id+","+line.split(",")[0]
                    if key in not_assigned:
                        z0=float(line.split(",")[-2])
                        z1=float(line.split(",")[-1])
                        z_dist=z0-z1
                        if z0>=z_thresh and z_dist<z_dist_thresh:
                            failed_z_dist=failed_z_dist+pool_id+","+line.strip()+","+panel_id+"\n"
    if len(failed_z_dist)>0:
        with open ("combined_gt_donor_assignments_failed_z_score_distance.csv", "w") as fout:
            fout.write('pool,donor_query,donor_gt,donor_1,score0,score1,score_n,n,mean,sd,z0,z1,panel\n'+failed_z_dist)


def combine_overall():
    not_assigned=[]
    file_list=glob.glob('stats_*_gt_donor_assignments.csv')
    with open ("combined_gt_donor_assignments_overall.csv", "w") as fout:
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
                    if line.strip().split(",")[-1]=="NONE":
                        not_assigned.append(pool_id+","+line.split(",")[0])
    return not_assigned
def main ():
    z_thresh = float(sys.argv[1])
    z_dist_thresh = float(sys.argv[2])

    not_assigned=combine_overall()
    combine_per_panel(not_assigned, z_thresh, z_dist_thresh)

if __name__ == '__main__':
    main()
