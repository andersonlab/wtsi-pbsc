import numpy as np
import pandas as pd



chrom='chr2'
unmapped_bed_f='/lustre/scratch126/humgen/projects/isogut/output/results/mapped_bam/Isogut14548275.mapped.realcells_only_true_unmapped.bed'
chrom_sizes=pd.read_csv('/lustre/scratch126/humgen/projects/isogut/utils/hg38.chrom.sizes',sep='\t',header=None)

chrom_length=chrom_sizes.loc[chrom_sizes[0]==chrom,1][1]
unmapped_dat=pd.read_csv(unmapped_bed_f,sep='\t',header=None)
unmapped_dat['length']=unmapped_dat[2]-unmapped_dat[1]
unmapped_dat['midpoint']=unmapped_dat[1]+np.floor( ((unmapped_dat[2]-unmapped_dat[1])/2) +1).astype(int)
chrom_midpoint=int(np.floor(chrom_length/2)+1)
unmapped_dat['dist_to_chrom_midpoint']=np.abs(unmapped_dat['midpoint']-chrom_midpoint)
print(unmapped_dat.sort_values(by=['dist_to_chrom_midpoint']))

