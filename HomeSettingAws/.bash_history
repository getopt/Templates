python scripts/get_rpm_read_starts_and_filter_by_fold_change.py m6A-TNFalpha2_S16.crosslink_sites_w_read_starts.bed m6A-TNFalpha2_S16.read_starts.bed RNAseq-TNFalpha2_S25.read_starts.bed 10 m6A-TNFalpha2_S16.crosslink_sites_w_read_starts_cutoff10.bed 
python scripts/get_rpm_read_starts_and_filter_by_fold_change.py m6A-TNFalpha2_S16.crosslink_sites_w_read_starts.bed m6A-TNFalpha2_S16.combined_w_uniquemap.rmDup.sam.parsed RNAseq-TNFalpha2_S25.combined_w_uniquemap.rmDup.sam.parsed 10 m6A-TNFalpha2_S16.crosslink_sites_w_read_starts_cutoff10.bed 
aws s3 cp m6A-TNFalpha2_S16.crosslink_sites_w_read_starts_cutoff10.bed s3://eclipsebio-turkey-m6a040319/single_nucleotide_analysis/pureclip_w_read_start_cutoffs_analysis/
clear
bedtools intersect -a m6A-TNFalpha3_S17.crosslink_sites.bed -b m6A-TNFalpha3_S17.read_starts.bed -s -c > m6A-TNFalpha3_S17.crosslink_sites_w_m6A_read_starts.bed 
bedtools intersect -a m6A-TNFalpha3_S17.crosslink_sites_w_m6A_read_starts.bed -b RNAseq-TNFalpha3_S26.read_starts.bed -s -c > m6A-TNFalpha3_S17.crosslink_sites_w_read_starts.bed 
aws s3 cp m6A-TNFalpha3_S17.crosslink_sites_w_read_starts.bed s3://eclipsebio-turkey-m6a040319/single_nucleotide_analysis/pureclip_w_read_start_cutoffs_analysis/
python scripts/get_rpm_read_starts_and_filter_by_fold_change.py m6A-TNFalpha3_S17.crosslink_sites_w_read_starts.bed m6A-TNFalpha3_S17.combined_w_uniquemap.rmDup.sam.parsed RNAseq-TNFalpha3_S26.combined_w_uniquemap.rmDup.sam.parsed 10 m6A-TNFalpha3_S17.crosslink_sites_w_read_starts_cutoff10.bed 
aws s3 cp m6A-TNFalpha3_S17.crosslink_sites_w_read_starts_cutoff10.bed s3://eclipsebio-turkey-m6a040319/single_nucleotide_analysis/pureclip_w_read_start_cutoffs_analysis/
rm m6A-*
clear
ls
rm RNAseq-*
clear
s
ls
logout
source activate eclipsebio
ls
clear
python scripts/get_DRACH_from_motif_heatmap.py m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_5.heatmap.csv 
clear
head m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_5.heatmap.csv 
clear
head -1 m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_5.heatmap.csv 
python scripts/get_DRACH_from_motif_heatmap.py m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_5.heatmap.csv 
clear
python scripts/get_DRACH_from_motif_heatmap.py m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_5.heatmap.csv m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_DRACH.tsv 
head m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_DRACH.tsv 
aws s3 cp m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_DRACH.tsv s3://eclipsebio-turkey-m6a040319/single_nucleotide_analysis/pureclip_w_read_start_cutoffs_analysis/
python scripts/get_DRACH_from_motif_heatmap.py m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_5.heatmap.csv m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_DRACH.tsv 
aws s3 cp m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_DRACH.tsv s3://eclipsebio-turkey-m6a040319/single_nucleotide_analysis/pureclip_w_read_start_cutoffs_analysis/
ls scripts/perl_scripts/
cp scripts/perl_scripts/motifs_by_read_density_heatmap.pl scripts/perl_scripts/motifs_by_read_density_heatmap_for_pureclip_sites.pl
vim scripts/perl_scripts/motifs_by_read_density_heatmap_for_pureclip_sites.pl 
aws s3 cp s3://eclipsebio-turkey-m6a040319/single_nucleotide_analysis/pureclip_w_read_start_cutoffs_analysis/m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.neg.bg .
aws s3 cp s3://eclipsebio-turkey-m6a040319/single_nucleotide_analysis/pureclip_w_read_start_cutoffs_analysis/m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg .
clear
mkdir eclipsebio-turkey-pureclip
mv *bg eclipsebio-turkey-pureclip/
perl scripts/perl_scripts/motifs_by_read_density_heatmap_for_pureclip_sites.pl /home/ec2-user/eclipsebio-turkey-pureclip/ /home/ec2-user/eclipsebio-turkey-pureclip/m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg /home/ec2-user/eclipsebio-turkey-pureclip/m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.neg.bg 5 hg38
ls eclipsebio-turkey-pureclip/
aws s3 cp s3://eclipsebio-turkey-m6a040319/single_nucleotide_analysis/pureclip_w_read_start_cutoffs_analysis/m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_5.all_output.csv .
aws s3 cp s3://eclipsebio-turkey-m6a040319/single_nucleotide_analysis/pureclip_w_read_start_cutoffs_analysis/m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_5.heatmap.csv .
diff m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_5.all_output.csv eclipsebio-turkey-pureclip/m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_5.all_output.csv 
diff m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_5.heatmap.csv eclipsebio-turkey-pureclip/m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_5.heatmap.csv 
clera
clear
grep 'GGACT' eclipsebio-turkey-pureclip/m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_5.heatmap.csv 
clear
grep 'AAACA' eclipsebio-turkey-pureclip/m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_5.heatmap.csv 
clear
ls
rm m6A-*
clear
ls
rm -r eclipsebio-turkey-pureclip/
clear
ls
logout
source activate eclipsebio
aws s3 cp s3://eclipsebio-turkey-m6a040319/single_nucleotide_analysis/pureclip_w_read_start_cutoffs_analysis/m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.pos.bg.kmer_5.heatmap.csv .
cp scripts/create_summary.py scripts/get_DRACH_from_motif_heatmap.py
vim scripts/get_DRACH_from_motif_heatmap.py 
clear
aws s3 cp s3://eclipsebio-turkey-m6a040319/single_nucleotide_analysis/pureclip_w_read_start_cutoffs_analysis/ . --recursive --exclude "*" --include "*cutoff10.bed"
clear
grep '+' m6A-CP2_S13.crosslink_sites_w_read_starts_cutoff10.bed | awk '{print $1,$2,$3,1}' OFS='\t' > m6A-CP2_S13.crosslink_sites_w_read_starts_cutoff10.pos.bg
grep '-' m6A-CP2_S13.crosslink_sites_w_read_starts_cutoff10.bed | awk '{print $1,$2,$3,1}' OFS='\t' > m6A-CP2_S13.crosslink_sites_w_read_starts_cutoff10.ne.bg
grep '+' m6A-CP3_S14.crosslink_sites_w_read_starts_cutoff10.bed | awk '{print $1,$2,$3,1}' OFS='\t' > m6A-CP3_S14.crosslink_sites_w_read_starts_cutoff10.pos.bg
grep '-' m6A-CP3_S14.crosslink_sites_w_read_starts_cutoff10.bed | awk '{print $1,$2,$3,1}' OFS='\t' > m6A-CP3_S14.crosslink_sites_w_read_starts_cutoff10.neg.bg
grep '+' m6A-DMSO1_S9.crosslink_sites_w_read_starts_cutoff10.bed | awk '{print $1,$2,$3,1}' OFS='\t' > m6A-DMSO1_S9.crosslink_sites_w_read_starts_cutoff10.pos.bg
grep '-' m6A-DMSO1_S9.crosslink_sites_w_read_starts_cutoff10.bed | awk '{print $1,$2,$3,1}' OFS='\t' > m6A-DMSO1_S9.crosslink_sites_w_read_starts_cutoff10.neg.bg
grep '+' m6A-DMSO2_S10.crosslink_sites_w_read_starts_cutoff10.bed | awk '{print $1,$2,$3,1}' OFS='\t' > m6A-DMSO2_S10.crosslink_sites_w_read_starts_cutoff10.pos.bg
grep '-' m6A-DMSO2_S10.crosslink_sites_w_read_starts_cutoff10.bed | awk '{print $1,$2,$3,1}' OFS='\t' > m6A-DMSO2_S10.crosslink_sites_w_read_starts_cutoff10.neg.bg
clear
grep '+' m6A-DMSO3_S11.crosslink_sites_w_read_starts_cutoff10.bed | awk '{print $1,$2,$3,1}' OFS='\t' > m6A-DMSO3_S11.crosslink_sites_w_read_starts_cutoff10.pos.bg
grep '-' m6A-DMSO3_S11.crosslink_sites_w_read_starts_cutoff10.bed | awk '{print $1,$2,$3,1}' OFS='\t' > m6A-DMSO3_S11.crosslink_sites_w_read_starts_cutoff10.neg.bg
grep '+' m6A-TNFalpha1_S15.crosslink_sites_w_read_starts_cutoff10.bed | awk '{print $1,$2,$3,1}' OFS='\t' > m6A-TNFalpha1_S15.crosslink_sites_w_read_starts_cutoff10.pos.bg
grep '-' m6A-TNFalpha1_S15.crosslink_sites_w_read_starts_cutoff10.bed | awk '{print $1,$2,$3,1}' OFS='\t' > m6A-TNFalpha1_S15.crosslink_sites_w_read_starts_cutoff10.neg.bg
grep '+' m6A-TNFalpha2_S16.crosslink_sites_w_read_starts_cutoff10.bed | awk '{print $1,$2,$3,1}' OFS='\t' > m6A-TNFalpha2_S16.crosslink_sites_w_read_starts_cutoff10.pos.bg
grep '-' m6A-TNFalpha2_S16.crosslink_sites_w_read_starts_cutoff10.bed | awk '{print $1,$2,$3,1}' OFS='\t' > m6A-TNFalpha2_S16.crosslink_sites_w_read_starts_cutoff10.neg.bg
grep '+' m6A-TNFalpha3_S17.crosslink_sites_w_read_starts_cutoff10.bed | awk '{print $1,$2,$3,1}' OFS='\t' > m6A-TNFalpha3_S17.crosslink_sites_w_read_starts_cutoff10.pos.bg
grep '-' m6A-TNFalpha3_S17.crosslink_sites_w_read_starts_cutoff10.bed | awk '{print $1,$2,$3,1}' OFS='\t' > m6A-TNFalpha3_S17.crosslink_sites_w_read_starts_cutoff10.neg.bg
clear
ls *bg
aws s3 cp . s3://eclipsebio-turkey-m6a040319/single_nucleotide_analysis/pureclip_w_read_start_cutoffs_analysis/ --recursive --exclude "*" --include "m6A-*crosslink_sites_w_read_starts_cutoff10*bg"
logout
source activate eclipsebio
aws s3 cp s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-newNCmembrane101619/data_analysis/m6A-5IP400-UV_enriched_5_4m6A_newNC_IP1_not_enriched.bed.annotated_proxdist_miRlncRNA .
clear
head m6A-5IP400-UV_enriched_5_4m6A_newNC_IP1_not_enriched.bed.annotated_proxdist_miRlncRNA 
clear
sort -nrk7 m6A-5IP400-UV_enriched_5_4m6A_newNC_IP1_not_enriched.bed.annotated_proxdist_miRlncRNA | head
clear
sort -nk8 m6A-5IP400-UV_enriched_5_4m6A_newNC_IP1_not_enriched.bed.annotated_proxdist_miRlncRNA | head
mkdir files_for_dropbox
aws s3 cp s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-newNCmembrane101619/data_analysis/m6A-7IP400-UV_enriched_7_m6A_newNC_IP2_not_enriched.bed.annotated_proxdist_miRlncRNA files_for_dropbox/
aws s3 cp s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-newNCmembrane101619/data_analysis/m6A-8IP400-UV_enriched_8_m6A_newNC_IP3_not_enriched.bed.annotated_proxdist_miRlncRNA files_for_dropbox/
clear
ls files_for_dropbox/
mv m6A-5IP400-UV_enriched_5_4m6A_newNC_IP1_not_enriched.bed.annotated_proxdist_miRlncRNA files_for_dropbox/
clear
ls files_for_dropbox/
source activate rclone_env
rclone lsd remote:
rclone lsd remote:"Eclipse Bio - m6A data folder"
rclone copy files_for_dropbox/ remote:"Eclipse Bio - m6A data folder"/m6A_newNC_101819/
wc -l files_for_dropbox/*
clear
sort -nrk5 files_for_dropbox/m6A-5IP400-UV_enriched_5_4m6A_newNC_IP1_not_enriched.bed.annotated_proxdist_miRlncRNA | head
clear
sort -nk8 files_for_dropbox/m6A-5IP400-UV_enriched_5_4m6A_newNC_IP1_not_enriched.bed.annotated_proxdist_miRlncRNA | head
source activate eclipsebio
python scripts/create_encode_idr_html_report.py files_for_dropbox/m6A-7IP400-UV_enriched_7_m6A_newNC_IP2_not_enriched.bed.annotated_proxdist_miRlncRNA 
aws s3 cp m6A-7IP400-UV_enriched_7_m6A_newNC_IP2_not_enriched.report.html s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-newNCmembrane101619/data_analysis/
clear
python scripts/create_encode_idr_html_report.py files_for_dropbox/m6A-8IP400-UV_enriched_8_m6A_newNC_IP3_not_enriched.bed.annotated_proxdist_miRlncRNA 
aws s3 cp m6A-8IP400-UV_enriched_8_m6A_newNC_IP3_not_enriched.report.html s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-newNCmembrane101619/data_analysis/
vim scripts/perl_scripts/motifs_by_read_density_heatmap.pl 
clear
ls
rm m6A-*
rm -r files_for_dropbox/
clear
ls
ls scripts/perl_scripts/
clear
ls scripts/perl_scripts/
clear
ls
logout
mkdir files_for_dropbox
source activate rclone_env
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/data_analysis/ . --recursive --exclude "*" --include "*d*best_peaks.bed"
clear
ls
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/normalization_outputs_w_3rep_idr/IDR/Lippi-dIP_P18_vs_input.01v02.IDR.out .
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/normalization_outputs_w_3rep_idr/IDR/Lippi-dIP_P18_vs_input.01v03.IDR.out .
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/normalization_outputs_w_3rep_idr/IDR/Lippi-dIP_P18_vs_input.02v03.IDR.out .
cp scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py files_for_dropbox/
cp scripts/perl_scripts/run_and_parse_IDR_3rep.pl files_for_dropbox/
vim scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py 
ls
ls scripts/
ls scripts/perl_scripts/
vim scripts/perl_scripts/get_reproducing_peaks.pl 
vim scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py 
cp scripts/perl_scripts/get_reproducing_peaks.pl scripts/perl_scripts/get_reproducing_peaks_3rep.pl 
vim scripts/perl_scripts/get_reproducing_peaks_3rep.pl 
clear
ls
ls scripts/
vim scripts/perl_scripts/get_reproducing_peaks_3rep.pl 
logout
source activate eclipsebio
mkdir eclipsebio-lippi100819
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/Lippi-dIP5_S5_L001_R1_001.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam eclipsebio-lippi100819/
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/Lippi-dIP5_S5_L001_R1_001.adapterTrim.round2.rmRep.sorted.rmDup.sorted.peaks.bed eclipsebio-lippi100819/
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/Lippi-dIP6_S6_L001_R1_001.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam eclipsebio-lippi100819/
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/Lippi-dIP6_S6_L001_R1_001.adapterTrim.round2.rmRep.sorted.rmDup.sorted.peaks.bed eclipsebio-lippi100819/
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/Lippi-dIP7_S7_L001_R1_001.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam eclipsebio-lippi100819/
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/Lippi-dIP7_S7_L001_R1_001.adapterTrim.round2.rmRep.sorted.rmDup.sorted.peaks.bed eclipsebio-lippi100819/
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt eclipsebio-lippi100819/
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/Lippi-inp1_S9_L001_R1_001.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam eclipsebio-lippi100819/
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/Lippi-inp2_S10_L001_R1_001.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam eclipsebio-lippi100819/
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/Lippi-inp3_S11_L001_R1_001.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam eclipsebio-lippi100819/
clear
ls eclipsebio-lippi100819/
python scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py /home/ec2-user/eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt /home/ec2-user/normalization_outputs_redo/ mm10
exit
vim scripts/perl_scripts/get_reproducing_peaks.pl 
logout
source activate eclipsebio
python scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py /home/ec2-user/eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt /home/ec2-user/normalization_outputs_redo/ mm10
wc -l normalization_outputs_redo/IDR/Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed
clear
sort normalization_outputs_redo/IDR/Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed
sort normalization_outputs_redo/IDR/Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed | uniq | wc -l
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/normalization_outputs_w_3rep_idr/IDR/Lippi-IP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed .
exit
screen
ls normalization_outputs_redo/
clear
ls normalization_outputs_redo/IDR/
screen -r
vim scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py 
cd normalization_outputs_redo/IDR/
wc -l *out
cd ../../
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/normalization_outputs_w_3rep_idr/IDR/Lippi-IP_P18_vs_input.01v02.IDR.out .
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/normalization_outputs_w_3rep_idr/IDR/Lippi-IP_P18_vs_input.01v03.IDR.out .
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/normalization_outputs_w_3rep_idr/IDR/Lippi-IP_P18_vs_input.02v03.IDR.out .
wc -l *out
ls
vim scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py 
clear
ls normalization_outputs_redo/IDR/
wc -l normalization_outputs_redo/IDR/Lippi-dIP_P18_vs_input.01v02.02v03.IDR.out.overlaps.txt 
wc -l normalization_outputs_redo/IDR/Lippi-dIP_P18_vs_input.IDR.out.overlaps.all.txt 
wc -l Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed.annotated_proxdist_miRlncRNA
wc -l normalization_outputs_redo/IDR/Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed.annotated_proxdist_miRlncRNA
vim scripts/perl_scripts/get_reproducing_peaks_3rep.pl 
clear
ls normalization_outputs_redo/IDR/
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/normalization_outputs_w_3rep_idr/IDR/Lippi-IP_P18_vs_input.IDR.out.overlaps.all.txt .
wc -l Lippi-IP_P18_vs_input.IDR.out.overlaps.all.txt 
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/normalization_outputs_w_3rep_idr/IDR/Lippi-dIP_P18_vs_input.IDR.out.overlaps.all.txt .
wc -l Lippi-dIP_P18_vs_input.IDR.out.overlaps.all.txt 
wc -l normalization_outputs_redo/IDR/Lippi-dIP_P18_vs_input.IDR.out.overlaps.all.txt 
vim scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py 
vim scripts/perl_scripts/get_reproducing_peaks_3rep.pl 
vim scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py 
head Lippi-dIP_P18_vs_input.IDR.out.overlaps.all.txt 
head normalization_outputs_redo/IDR/Lippi-dIP_P18_vs_input.IDR.out.overlaps.all.txt 
vim scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py 
vim scripts/perl_scripts/get_reproducing_peaks_3rep.pl 
clear
ls normalization_outputs_redo/IDR/
wc -l normalization_outputs_redo/IDR/*RNA
vim scripts/perl_scripts/get_reproducing_peaks_3rep.pl 
source activate eclipsebio
rm -r normalization_outputs_redo/*
clea
clear
python scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py /home/ec2-user/eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt /home/ec2-user/normalization_outputs_redo/ mm10
rm -r normalization_outputs_redo/*
rm eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt.*
clear
python scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py /home/ec2-user/eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt /home/ec2-user/normalization_outputs_redo/ mm10
rm -r normalization_outputs_redo/*
rm eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt.*
clear
python scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py /home/ec2-user/eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt /home/ec2-user/normalization_outputs_redo/ mm10
clear
ls normalization_outputs_redo/IDR/
rm -r normalization_outputs_redo/*
rm eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt.*
clear
python scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py /home/ec2-user/eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt /home/ec2-user/normalization_outputs_redo/ mm10
ls normalization_outputs_redo/IDR/
wc -l normalization_outputs_redo/IDR/Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed
rm -r normalization_outputs_redo/*
rm eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt.*
screen
clear
ls
cp normalization_outputs_redo/IDR/Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed new_Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed
new_Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed
clear
vim normalization_outputs_redo/IDR/Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed
vim scripts/perl_scripts/get_reproducing_peaks_3rep.pl 
ls
mkdir lippi_compressed_bed_files
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/normalization_outputs/ lippi_compressed_bed_files/ --recursive --exclude "*" --include "*compressed.bed"
cd lippi_compressed_bed_files/
clear
sort -nrk5 Lippi-dIP5_vs_dIP8_01.basedon.peaks.l2inputnormnew.bed.compressed.bed | head
clear
sort -nrk5 Lippi-dIP5_vs_dIP8_01.basedon.peaks.l2inputnormnew.bed.compressed.bed | head
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/normalization_outputs/Lippi-dIP5_vs_dIP8_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA .
clear
sort -nrk5 Lippi-dIP5_vs_dIP8_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA | head
cleafr
clear
sort -nrk12 Lippi-dIP5_vs_dIP8_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA | head
rm -r normalization_outputs_redo/*
clear
ls
ls eclipsebio-lippi100819/
rm eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt.*
clear
source activate eclipsebio
python scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py /home/ec2-user/eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt /home/ec2-user/normalization_outputs_redo/ mm10
source activate eclipsebio
python scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py /home/ec2-user/eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt /home/ec2-user/normalization_outputs_redo/ mm10
clear
ls
rm -r normalization_outputs_redo/*
rm eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt.*
clear
python scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py /home/ec2-user/eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt /home/ec2-user/normalization_outputs_redo/ mm10
exit
vim scripts/perl_scripts/get_reproducing_peaks_3rep.pl 
logout
rm -r normalization_outputs_redo/*
rm eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt.*
screen
cp normalization_outputs_redo/IDR/Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed rerun_Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed
clear
ls
source activate eclipsebio
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/normalization_outputs_w_3rep_idr/IDR/Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed .
bedtools intersect -a Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed -b rerun_Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed -v -s 
wc -l Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed 
wc -l rerun_Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed 
clear
vim Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed 
rm -r normalization_outputs_redo/*
rm eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt.*
python scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py /home/ec2-user/eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt /home/ec2-user/normalization_outputs_redo/ mm10
clear
rm -r normalization_outputs_redo/*
rm eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt.*
python scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py /home/ec2-user/eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt /home/ec2-user/normalization_outputs_redo/ mm10
rm -r normalization_outputs_redo/*
rm eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt.*
clear
python scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py /home/ec2-user/eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt /home/ec2-user/normalization_outputs_redo/ mm10
rm -r normalization_outputs_redo/*
rm eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt.*
clear
python scripts/perl_scripts/Peak_normalization_wrapper_updated_3rep.py /home/ec2-user/eclipsebio-lippi100819/Lippi-dIP_P18_vs_input.manifest.txt /home/ec2-user/normalization_outputs_redo/ mm10
cd eclipsebio-lippi100819/
clear
ls
aws s3 cp . s3://eclipsebio-scripps/eclipsebio-lippi100819/ --recursive --exclude "*" --include "Lippi-dIP_P18_vs_input.manifest.txt*"
cd ../normalization_outputs_redo/
clear
ls
rm *bam
rm *peaks.bed
clear
ls
aws s3 cp . s3://eclipsebio-scripps/eclipsebio-lippi100819/corrected_normalization_outputs_3rep/ --recursive
clear
ls
cd ../
ls
ls files_for_dropbox/
vim files_for_dropbox/run_and_parse_IDR_3rep.pl 
vim scripts/perl_scripts/merge_idr_entropy_3rep.pl 
cp scripts/perl_scripts/merge_idr_entropy_3rep.pl files_for_dropbox/
cp scripts/perl_scripts/get_reproducing_peaks_3rep.pl files_for_dropbox/
clear
ls
clear
ls scripts/perl_scripts/
source activate rclone_env
clear
ls
rm Lippi-*
rm rerun_Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed 
rm new_Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed 
clear
ls
rm -r eclipsebio-lippi100819/
ls lippi_compressed_bed_files/
ls normalization_outputs_redo/*compressed.bed
cp normalization_outputs_redo/*compressed.bed files_for_dropbox/
cp normalization_outputs_redo/IDR/Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed files_for_dropbox/
rclone lsd remote:
rclone lsd remote:"EclipseBioComputational"
rclone copy files_for_dropbox/ remote:"EclipseBioComputational"/IDR_3rep/
clear
cd normalization_outputs_redo/
wc -l *compressed.bed
cd IDR/
ls
wc -l Lippi-dIP_P18_vs_input.01v02.IDR.out
wc -l Lippi-dIP_P18_vs_input.02v03.IDR.out
wc -l Lippi-dIP_P18_vs_input.01v03.IDR.out
wc -l Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.
wc -l Lippi-dIP_P18_vs_input.01v02.IDR.out.all_reps.merged.bed
cd ../../
clear
ls
rm -r lippi_compressed_bed_files/
rm -r normalization_outputs_redo/
clear
ls
logout
source activate eclipsebio
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/miRNA_outputs/ . --recursive --exclude "*" --include "*full.compressed2.bed.full*RNA"
clear
ls
head mirbase_files/mmu.gff3 
cp mirbase_files/mmu.gff3 .
vim mmu.gff3 
clear
head mmu.gff3 
vim mmu.gff3 
head mmu.gff3 
awk {'print $1,$4,$5,$9,$3,$7}' OFS='\t' mmu.gff3 | head
clear
awk {'print $1,$4,$5,$9,$3,$7}' OFS='\t' mmu.gff3 | head
awk {'print $1,$4,$5,$9,$3,$7}' OFS='\t' mmu.gff3 > mmu.bed
head Lippi-IP1_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
head mmu.bed 
head Lippi-IP1_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
python scripts/parse_mirna_normalized_peaks_to_get_rpms.py mmu.bed output Lippi-IP1_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
python scripts/parse_mirna_normalized_peaks_to_get_rpms.py mmu.bed output Lippi-IP1_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
less mmu.bed 
clear
python scripts/parse_mirna_normalized_peaks_to_get_rpms.py mmu.bed output Lippi-IP1_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
python scripts/parse_mirna_normalized_peaks_to_get_rpms.py mmu.bed output Lippi-IP1_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
python scripts/parse_mirna_normalized_peaks_to_get_rpms.py mmu.bed Lippi_mirna_counts.tsv Lippi-IP1_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP2_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP3_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP4_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP5_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP6_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP7_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP8_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
python scripts/parse_mirna_normalized_peaks_to_get_rpms.py mmu.bed Lippi_mirna_counts.tsv Lippi-IP1_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP2_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP3_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP4_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP5_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP6_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP7_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP8_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
python scripts/parse_mirna_normalized_peaks_to_get_rpms.py mmu.bed Lippi_mirna_counts.tsv Lippi-IP1_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP2_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP3_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP4_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP5_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP6_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP7_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP8_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
head Lippi-IP1_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
python scripts/parse_mirna_normalized_peaks_to_get_rpms.py mmu.bed Lippi_mirna_counts.tsv Lippi-IP1_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP2_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP3_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP4_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP5_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP6_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP7_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP8_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
python scripts/parse_mirna_normalized_peaks_to_get_rpms.py mmu.bed Lippi_mirna_counts.tsv Lippi-IP1_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP2_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP3_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP4_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP5_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP6_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP7_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP8_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
python scripts/parse_mirna_normalized_peaks_to_get_rpms.py mmu.bed Lippi_mirna_counts.tsv Lippi-IP1_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP2_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP3_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP4_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP5_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP6_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP7_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP8_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
python scripts/parse_mirna_normalized_peaks_to_get_rpms.py mmu.bed Lippi_mirna_counts.tsv Lippi-IP1_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP2_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP3_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-IP4_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP5_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP6_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP7_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA Lippi-dIP8_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
ls
head Lippi_mirna_counts.tsv 
vim Lippi_mirna_counts.tsv 
grep 'mmu-mir-30a' Lippi-IP1_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
grep 'mmu-mir-30a' Lippi-IP2_miRNA_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
aws s3 cp Lippi_mirna_counts.tsv s3://eclipsebio-scripps/eclipsebio-lippi100819/miRNA_outputs/
clear
aws s3 cp . s3://eclipsebio-scripps/eclipsebio-lippi100819/ . --recursive --exclude "*" --include "*miRNA_manifest.txt.mapped_read_num"
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/ . --recursive --exclude "*" --include "*miRNA_manifest.txt.mapped_read_num"
aws s3 cp s3://eclipsebio-scripps/eclipsebio-lippi100819/ . --recursive --exclude "*" --include "*miRNA.manifest.txt.mapped_read_num"
clear
head Lippi-IP1_miRNA.manifest.txt.mapped_read_num 
head Lippi-dIP5_miRNA.manifest.txt.mapped_read_num 
clear
cat Lippi-IP2_miRNA.manifest.txt.mapped_read_num 
cat Lippi-IP3_miRNA.manifest.txt.mapped_read_num 
clear
cat Lippi-IP4_miRNA.manifest.txt.mapped_read_num 
clear
cat Lippi-dIP6_miRNA.manifest.txt.mapped_read_num 
cat Lippi-dIP7_miRNA.manifest.txt.mapped_read_num 
clear
cat Lippi-dIP8_miRNA.manifest.txt.mapped_read_num 
logout
cp scripts/create_summary.py scripts/parse_mirna_normalized_peaks_to_get_rpms.py
vim scripts/parse_mirna_normalized_peaks_to_get_rpms.py 
logout
ls
rm Lippi*
rm mmu.gff3 
mv mmu.bed mirbase_files/
clear
ls
logout
sudo yum update
source activate eclipsebio
aws s3 cp s3://eclipsebio-turkey-m6a040319/ . --recursive --exclude "*" --include "*bedgraph"
aws s3 cp s3://eclipsebio-turkey-m6a040319/ . --recursive --exclude "*" --include "*mapped_read_num"
clear
head m6A-CP1_S12.manifest.txt.mapped_read_num 
python scripts/calculate_log2foldchange_read_starts.py m6A-CP1_S12.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph RNAseq-CP1_S21.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph m6A-CP1_S12.manifest.txt.mapped_read_num m6A-CP1_S12.read_starts_log2fc.bed
clear
head m6A-CP1_S12.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph 
python scripts/calculate_log2foldchange_read_starts.py m6A-CP1_S12.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph RNAseq-CP1_S21.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph m6A-CP1_S12.manifest.txt.mapped_read_num m6A-CP1_S12.read_starts_log2fc.bedgraph
python scripts/calculate_foldchange_read_starts.py m6A-CP1_S12.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph RNAseq-CP1_S21.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph m6A-CP1_S12.manifest.txt.mapped_read_num m6A-CP1_S12.read_starts_log2fc.bedgraph
clear
python scripts/calculate_foldchange_read_starts.py m6A-CP1_S12.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph RNAseq-CP1_S21.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph m6A-CP1_S12.manifest.txt.mapped_read_num m6A-CP1_S12.read_starts_log2fc.bedgraph
cat m6A-CP1_S12.read_starts_log2fc.bedgraph 
grep 'chr1' m6A-CP1_S12.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph | grep '605529'
grep 'chr1' RNAseq-CP1_S21.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph | grep '605529'
cat m6A-CP1_S12.manifest.txt.mapped_read_num 
clear
python scripts/calculate_foldchange_read_starts.py m6A-CP1_S12.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph RNAseq-CP1_S21.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph m6A-CP1_S12.manifest.txt.mapped_read_num m6A-CP1_S12.read_starts_log2fc.bedgraph
clear
head m6A-CP1_S12.read_starts_log2fc.bedgraph 
head m6A-CP1_S12.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph 
clear
sort -nrk4 m6A-CP1_S12.read_starts_log2fc.bedgraph 
sort -nrk4 m6A-CP1_S12.read_starts_log2fc.bedgraph | head
grep 'chr15' m6A-CP1_S12.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph | grep '43792329' | grep '43792330'
grep 'chr15' RNAseq-CP1_S21.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph | grep '43792329' | grep '43792330'
head RNAseq-CP1_S21.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph 
grep 'chr15' RNAseq-CP1_S21.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph | grep '4379232'
clear
python scripts/calculate_foldchange_read_starts.py m6A-CP1_S12.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph RNAseq-CP1_S21.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph m6A-CP1_S12.manifest.txt.mapped_read_num m6A-CP1_S12.read_starts_log2fc.bedgraph
clear
aws s3 cp s3://eclipsebio-turkey-m6a040319/RNAseq-CP1_S21.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam .
aws s3 cp s3://eclipsebio-turkey-m6a040319/m6A-CP1_S12.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam .
bedtools genomecov -ibam m6A-CP1_S12.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam -dz -5 -strand + | head
bedtools genomecov -ibam m6A-CP1_S12.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam -dz -5 -strand + -bg | head
cler
clear
python scripts/calculate_foldchange_read_starts.py m6A-CP1_S12.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph RNAseq-CP1_S21.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph m6A-CP1_S12.manifest.txt.mapped_read_num m6A-CP1_S12.read_starts_log2fc.bedgraph
clear
python scripts/calculate_foldchange_read_starts.py m6A-CP1_S12.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph RNAseq-CP1_S21.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph m6A-CP1_S12.manifest.txt.mapped_read_num m6A-CP1_S12.read_starts_log2fc.bedgraph
clear
python scripts/calculate_foldchange_read_starts.py m6A-CP1_S12.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph RNAseq-CP1_S21.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph m6A-CP1_S12.manifest.txt.mapped_read_num m6A-CP1_S12.read_starts_log2fc.bedgraph
head m6A-CP1_S12.read_starts_log2fc.bedgraph 
clear
sort -nrk4 m6A-CP1_S12.read_starts_log2fc.bedgraph | head
grep 'chr15' m6A-CP1_S12.read_starts_log2fc.bedgraph | '3959561'
grep 'chr15' m6A-CP1_S12.read_starts_log2fc.bedgraph | grep '3959561'
sort -nrk4 m6A-CP1_S12.read_starts_log2fc.bedgraph | less
vim scripts/perl_scripts/motifs_by_read_density_heatmap_for_pureclip_sites.pl 
aws s3 cp s3://eclipsebio-turkey-m6a040319/single_nucleotide_analysis/pureclip_w_read_start_cutoffs_analysis/m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.bed .
head m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.bed 
awk '{print $1,$2+1,$3+1}' m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.bed | head
awk '{print $1,$2+1,$3+1,$4,$5,$6}' m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.bed | head
awk '{print $1,$2+1,$3+1,$4,$5,$6}'OFS='\t' m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.bed | head
awk '{print $1,$2+1,$3+1,$4,$5,$6}' OFS='\t' m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.bed | head
awk '{print $1,$2+1,$3+1,$4,$5,$6}' OFS='\t' m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.bed > edited_m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.bed 
python scripts/calculate_foldchange_read_starts.py m6A-CP1_S12.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph RNAseq-CP1_S21.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.5.pos.bedgraph m6A-CP1_S12.manifest.txt.mapped_read_num m6A-CP1_S12.read_starts_log2fc.bedgraph
head m6A-CP1_S12.read_starts_log2fc.bedgraph 
awk '{print $1,$2,$3,".",$4,"+"}' OFS='\t' m6A-CP1_S12.read_starts_log2fc.bedgraph | head
awk '{print $1,$2,$3,".",$4,"+"}' OFS='\t' m6A-CP1_S12.read_starts_log2fc.bedgraph > edited_m6A-CP1_S12.read_starts_log2fc.bedgraph 
awk '{print $1,$2,$3,".",$4,"+"}' OFS='\t' m6A-CP1_S12.read_starts_log2fc.bedgraph > edited_m6A-CP1_S12.read_starts_log2fc.bed 
rm edited_m6A-CP1_S12.read_starts_log2fc.bedgraph 
bedtools intersect -a edited_m6A-CP1_S12.read_starts_log2fc.bed edited_m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.bed -v > m6A-CP1_S12.enriched_read_starts_no_pureclip.bed
bedtools intersect -a edited_m6A-CP1_S12.read_starts_log2fc.bed -b edited_m6A-CP1_S12.crosslink_sites_w_read_starts_cutoff10.bed -v > m6A-CP1_S12.enriched_read_starts_no_pureclip.bed
head m6A-CP1_S12.enriched_read_starts_no_pureclip.bed 
wc -l m6A-CP1_S12.enriched_read_starts_no_pureclip.bed 
clear
sort -nrk5 m6A-CP1_S12.enriched_read_starts_no_pureclip.bed | head
logout
ls scripts/
cp scripts/create_summary.py scripts/calculate_log2foldchange_read_starts.py
vim scripts/calculate_
vim scripts/calculate_log2foldchange_read_starts.py 
mv scripts/calculate_log2foldchange_read_starts.py scripts/calculate_foldchange_read_starts.py
vim scripts/calculate_foldchange_read_starts.py 
clear
ls
vim scripts/calculate_foldchange_read_starts.py 
rm RNAseq-*
rm m6A-*
clear
ls
rm edited_m6A-CP1_S12.*
clear
ls
ls single_nucleotide_analysis/
rm -r single_nucleotide_analysis/
clear
ls
logout
source activate eclipsebio
pureclip
conda env list
conda create --name pureclip_env
source activate pureclip_env
conda install pureclip
clear
aws s3 cp s3://eclipsebio-m6a051419/m12-IP1_S1.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam .
aws s3 cp s3://eclipsebio-m6a051419/m12-IP1_S1.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.bai .
aws s3 cp s3://eclipsebio-m6a051419/m12-inp1_S11.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam .
aws s3 cp s3://eclipsebio-m6a051419/m12-inp1_S11.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.bai .
clear
pureclip -h
clear
pureclip -i m12-IP1_S1.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam -bai m12-IP1_S1.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.bai -g genomes/fasta_files/mm10.fa -nt 36 -o m12-IP1_S1.crosslink_sites.bed -ibam m12-inp1_S11.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam -ibai m12-inp1_S11.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam.bai 
clear
ls
rm m12-*
clear
ls
logout
cp scripts/create_summary.py scripts/switch_strand.py
vim scripts/switch_strand.py 
source activate eclipsebio
samtools index EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam 
aws s3 cp EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam.bai s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/
logout
cd rep_element_pipeline/eclipsebio-m6a-neb-kit100919/
ls
vim RM_EpiMark_7IP_S51_L002_R1_001.sh 
aws s3 cp s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/rep_element_pipeline/EpiMark-7IP_S51_L002_R1_001.fastq.gz.mapped_vs_bt2_hg38.sam.combined_w_uniquemap.prermDup.sam.gz .
gunzip EpiMark-7IP_S51_L002_R1_001.fastq.gz.mapped_vs_bt2_hg38.sam.combined_w_uniquemap.prermDup.sam.gz 
head EpiMark-7IP_S51_L002_R1_001.fastq.gz.mapped_vs_bt2_hg38.sam.combined_w_uniquemap.prermDup.sam 
clear
cut -f 2 EpiMark-7IP_S51_L002_R1_001.fastq.gz.mapped_vs_bt2_hg38.sam.combined_w_uniquemap.prermDup.sam | sort | uniq
vim RM_EpiMark_7IP_S51_L002_R1_001.sh 
logout
source activate eclipsebio
aws s3 cp s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam .
samtools vie EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam | head
samtools view EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam | head
samtools view -h EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam > EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.sam 
python scripts/switch_strand.py EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.sam EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.sam
clear
python scripts/switch_strand.py EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.sam EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.sam
head EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.sam 
head EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.sam
python scripts/switch_strand.py EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.sam EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.sam
clear
samtools view -bS EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.sam > EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam
samtools index EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam 
aws s3 cp EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/
aws s3 cp EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam.bai s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/
aws s3 cp s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/EpiMark-7input_S49_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam .
samtools view -h EpiMark-7input_S49_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam > EpiMark-7input_S49_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.sam 
python scripts/switch_strand.py EpiMark-7input_S49_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.sam EpiMark-7input_S49_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.sam
samtools view -bS EpiMark-7input_S49_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.sam > EpiMark-7input_S49_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam 
samtools index EpiMark-7input_S49_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam
samtools index EpiMark-7input_S49_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam
clear
ls EpiMark-7in*
aws s3 cp EpiMark-7input_S49_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam.bai s3://eclipsebio-internal001/eclipsebio-slamseq101419/CIMS_analysis/
aws s3 cp EpiMark-7input_S49_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam s3://eclipsebio-internal001/eclipsebio-slamseq101419/CIMS_analysis/
aws s3 cp EpiMark-7input_S49_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam.bai s3://eclipsebio-internal001/eclipsebio-slamseq101419/CIMS_analysis/
aws s3 cp EpiMark-7input_S49_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam.bai s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/
aws s3 cp EpiMark-7input_S49_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/
aws s3 cp EpiMark-7input_S49_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam.bai s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/
aws s3 cp s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/EpiMark-8IP_S52_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam .
clear
samtools view -h EpiMark-8IP_S52_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam > EpiMark-8IP_S52_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.sam 
python scripts/switch_strand.py EpiMark-8IP_S52_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.sam EpiMark-8IP_S52_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.sam
samtools view -bS EpiMark-8IP_S52_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.sam > EpiMark-8IP_S52_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam 
samtools index EpiMark-8IP_S52_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam
samtools index EpiMark-8IP_S52_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam
aws s3 cp EpiMark-8IP_S52_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam.bai s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/
aws s3 cp EpiMark-8IP_S52_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/
aws s3 cp EpiMark-8IP_S52_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam.bai s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/
aws s3 cp s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/EpiMark-8input_S50_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam .
clear
samtools view -bS EpiMark-8input_S50_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam > EpiMark-8input_S50_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.sam 
python scripts/switch_strand.py EpiMark-8input_S50_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.sam EpiMark-8input_S50_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.sam
samtools view -h EpiMark-8input_S50_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam > EpiMark-8input_S50_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.sam 
clear
python scripts/switch_strand.py EpiMark-8input_S50_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.sam EpiMark-8input_S50_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.sam
aws s3 cp scripts/switch_strand.py s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/
samtools view -bS EpiMark-8input_S50_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.sam > EpiMark-8input_S50_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam
samtools index EpiMark-8input_S50_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam
samtools index EpiMark-8input_S50_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam
aws s3 cp EpiMark-8input_S50_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.bam.bai s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/
aws s3 cp EpiMark-8input_S50_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/
aws s3 cp EpiMark-8input_S50_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam.bai s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/
clear
cd rep_element_pipeline/
ls
mkdir eclipsebio-m6a-neb-kit100919
aws s3 cp s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/rep_element_pipeline/EpiMark_7IP_S51_L002_R1_001.rep_element_manifest.txt eclipsebio-m6a-neb-kit100919/
cd eclipsebio-m6a-neb-kit100919/
vim EpiMark_7IP_S51_L002_R1_001.rep_element_manifest.txt 
cd ../
perl scripts/wrapper_FULLversion_includemultifamily_submit.pl /home/ec2-user/rep_element_pipeline/eclipsebio-m6a-neb-kit100919/ /home/ec2-user/rep_element_pipeline/eclipsebio-m6a-neb-kit100919/EpiMark_7IP_S51_L002_R1_001.rep_element_manifest.txt /home/ec2-user/rep_element_pipeline/eclipsebio-m6a-neb-kit100919/ hg38
cd eclipsebio-m6a-neb-kit100919/
vim RM_EpiMark_7IP_S51_L002_R1_001.sh 
cd ../
clear
vim Snakefile 
clear
python scripts/make_bigwig_files.py --bam EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam --genome genomes/chrom_sizes/hg38.chrom.sizes --bw_pos EpiMark-7IP_S51_L002_R1_001.switch_strand.CombinedID.merged.r2.norm.pos.bw --bw_neg EpiMark-7IP_S51_L002_R1_001.switch_strand.CombinedID.merged.r2.norm.neg.bw
python scripts/make_bigwig_files.py EpiMark-7input_S49_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam --genome genomes/chrom_sizes/hg38.chrom.sizes --bw_pos EpiMark-7input_S49_L002_R1_001.switch_strand.CombinedID.merged.r2.norm.pos.bw --bw_neg EpiMark-7input_S49_L002_R1_001.switch_strand.CombinedID.merged.r2.norm.neg.bw
python scripts/make_bigwig_files.py --bam EpiMark-7input_S49_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam --genome genomes/chrom_sizes/hg38.chrom.sizes --bw_pos EpiMark-7input_S49_L002_R1_001.switch_strand.CombinedID.merged.r2.norm.pos.bw --bw_neg EpiMark-7input_S49_L002_R1_001.switch_strand.CombinedID.merged.r2.norm.neg.bw
python scripts/make_bigwig_files.py --bam EpiMark-8IP_S52_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam --genome genomes/chrom_sizes/hg38.chrom.sizes --bw_pos EpiMark-8IP_S52_L002_R1_001.switch_strand.CombinedID.merged.r2.norm.pos.bw --bw_neg EpiMark-8IP_S52_L002_R1_001.switch_strand.CombinedID.merged.r2.norm.neg.bw
python scripts/make_bigwig_files.py --bam EpiMark-8input_S50_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam --genome genomes/chrom_sizes/hg38.chrom.sizes --bw_pos EpiMark-8input_S50_L002_R1_001.switch_strand.CombinedID.merged.r2.norm.pos.bw --bw_neg EpiMark-8input_S50_L002_R1_001.switch_strand.CombinedID.merged.r2.norm.neg.bw
clear
ls *bw
aws s3 cp . s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/ --recursive --exclude "*" --include "EpiMark-*bw"
cd rep_element_pipeline/eclipsebio-m6a-neb-kit100919/
aws s3 cp s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.fastq.gz .
aws s3 cp s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam .
vim RM_EpiMark_7IP_S51_L002_R1_001.sh 
gunzip EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.fastq.gz 
cp ../scripts/parse_bowtie2_output_realtime_includemultifamily.pl ../scripts/parse_bowtie2_output_realtime_includemultifamily_switch_strand.pl 
vim ../scripts/parse_bowtie2_output_realtime_includemultifamily_switch_strand.pl 
clear
perl ../scripts/parse_bowtie2_output_realtime_includemultifamily_switch_strand.pl /home/ec2-user/rep_element_pipeline/eclipsebio-m6a-neb-kit100919/EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.fastq /home/ec2-user/rep_element_pipeline/bt2_hg38 /home/ec2-user/rep_element_pipeline/eclipsebio-m6a-neb-kit100919/EpiMark-7IP_S51_L002_R1_001.fastq.gz.mapped_vs_bt2_hg38.sam /home/ec2-user/rep_element_pipeline/annotation_files/homo_sapiens_repeats_dfam_enst2id.txt
vim ../scripts/parse_bowtie2_output_realtime_includemultifamily_switch_strand.pl 
ls
clear
ls
rm *sam*
clear
ls
perl ../scripts/parse_bowtie2_output_realtime_includemultifamily_switch_strand.pl /home/ec2-user/rep_element_pipeline/eclipsebio-m6a-neb-kit100919/EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.fastq /home/ec2-user/rep_element_pipeline/bt2_hg38 /home/ec2-user/rep_element_pipeline/eclipsebio-m6a-neb-kit100919/EpiMark-7IP_S51_L002_R1_001.fastq.gz.mapped_vs_bt2_hg38.sam /home/ec2-user/rep_element_pipeline/annotation_files/homo_sapiens_repeats_dfam_enst2id.txt
clear
perl ../scripts/duplicate_removal_inline_paired_count_region_other_reads.pl /home/ec2-user/rep_element_pipeline/eclipsebio-m6a-neb-kit100919/EpiMark-7IP_S51_L002_R1_001.fastq.gz.mapped_vs_bt2_hg38.sam /home/ec2-user/rep_element_pipeline/eclipsebio-m6a-neb-kit100919/EpiMark-7IP_S51_L002_R1_001.adapterTrim.round2.rmRep.sorted.mapped_only.switch_strand.bam /home/ec2-user/genomes/gtf_files/gencode.v31.chr_patch_hapl_scaff.annotation.gtf /home/ec2-user/genomes/gtf_files/gencode.v31.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat ../annotation_files/RepeatMask.bed /home/ec2-user/rep_element_pipeline/annotation_files/homo_sapiens_repeats_dfam_enst2id.txt 
clear
ls
aws s3 cp RM_EpiMark_7IP_S51_L002_R1_001.sh s3://eclipsebio-internal001/eclipsebio-m6a/eclipsebio-m6a-neb-kit100919/rep_element_pipeline/
ls ../scripts/
ls ../../scripts/
ckear
ls
ckear
ks
clear
ls
rm *
cd ../
clear
ls
rm EpiMark-*
clear
ls
logout
gunzip gencode.v32.chr_patch_hapl_scaff.annotation.gtf.gz 
mv gencode.v32.chr_patch_hapl_scaff.annotation.gtf genomes/gtf_files/
clear
ls
rm -r files_for_dropbox/
ls genomes/gtf_files/
gunzip gencode.v32lift37.annotation.gtf.gz 
mv gencode.v32lift37.annotation.gtf genomes/gtf_files/
head genomes/gtf_files/gencode.v32lift37.annotation.gtf 
:q
clear
ls genomes/gtf_files/
vim scripts/perl_scripts/annotate_peaks_bedformat_wproxdistal_lncRNA.pl 
vim scripts/perl_scripts/annotate_peaks_fullformat_wproxdistal_lncRNA.pl 
vim scripts/perl_scripts/calculate_saturation_by_transcript_eachRBPseparately_submit.pl 
vim scripts/perl_scripts/count_reads_broadfeatures_frombamfi_SRmap.pl 
vim scripts/perl_scripts/find_genomic_features.pl 
vim rep_element_pipeline/scripts/wrapper_FULLversion_includemultifamily_submit.pl 
vim scripts/perl_scripts/annotate_peaks_bedformat_wproxdistal_lncRNA.pl 
vim scripts/perl_scripts/annotate_peaks_fullformat_wproxdistal_lncRNA.pl 
vim scripts/perl_scripts/calculate_saturation_by_transcript_eachRBPseparately_submit.pl 
vim scripts/perl_scripts/count_reads_broadfeatures_frombamfi_SRmap.pl 
vim scripts/perl_scripts/find_genomic_features.pl 
vim rep_element_pipeline/scripts/wrapper_FULLversion_includemultifamily_submit.pl 
clear
ls
clear
ls
gunzip gencode.vM23.chr_patch_hapl_scaff.annotation.gtf.gz
clear
ls
mv gencode* genomes/gtf_files/
clear
ls
vim scripts/perl_scripts/annotate_peaks_bedformat_wproxdistal_lncRNA.pl 
vim scripts/perl_scripts/annotate_peaks_fullformat_wproxdistal_lncRNA.pl 
vim scripts/perl_scripts/calculate_saturation_by_transcript_eachRBPseparately_submit.pl 
vim scripts/perl_scripts/count_reads_broadfeatures_frombamfi_SRmap.pl 
vim scripts/perl_scripts/find_genomic_features.pl 
vim rep_element_pipeline/scripts/wrapper_FULLversion_includemultifamily_submit.pl 
ckear
ks
clear
ls
logout
sudo yum update
source activate eclipsebio
vim HTML_Snakefile 
vim print_html_snakemake.sh 
vim run_html_snakemake.sh 
./print_html_snakemake.sh 
clear
./run_html_snakemake.sh 
aws s3 cp s3://eclipsebio-customers2019/eclipsebio-skyhawk0072-111519/data_analysis/0072-IP6-Treatment3_vs_IP3_DMSO3.report.html .
mv 0072-IP6-Treatment3_vs_IP3_DMSO3.report.html 0072-IP6-Treatment3_vs_IP3-DMSO3.report.html
aws s3 cp 0072-IP6-Treatment3_vs_IP3-DMSO3.report.html s3://eclipsebio-customers2019/eclipsebio-skyhawk0072-111519/data_analysis/
rm 0072-IP6-Treatment3_vs_IP3-DMSO3.report.html 
vim Snakefile 
clear
ls
rm -r eclipsebio-customers2019/
cd rep_element_pipeline/
ls
rm -r eclipsebio-m6a-neb-kit100919/
cd ../
clear
ls
vim HTML_Snakefile 
vim scripts/create_encode_idr_html_report.py 
vim scripts/create_edited_rep_element_html_report.py 
clear
ls
logout
source activate eclipsebio
python scripts/perl_scripts/Peak_normalization_wrapper.py /home/ec2-user/eclipsebio-m6a051419/m12-IP6_w_all_IP4_peaks.manifest.txt /home/ec2-user/normalization_outputs/ mm10
cd eclipsebio-m6a051419/
ls
aws s3 cp . s3://eclipsebio-m6a051419/ --recursive --exclude "*" --include "m12-IP6_w_all_IP4_peaks.manifest.txt*"
cd ../normalization_outputs/
clear
ls
rm -r IDR/
rm *bam
rm m12-IP6_S6.adapterTrim.round2.rmRep.sorted.rmDup.sorted.peaks.bed 
aws s3 cp . s3://eclipsebio-m6a051419/normalization_outputs/ --recursive
cd ../
python scripts/get_full_read_info_IPvIP_peaks.py m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP4_S4_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA normalization_outputs/m12-IP6_w_all_IP4_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA
head m12-IP4_vs_IP6_peaks_w_full_read_nums.tsv 
grep 'chr18' m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA | grep '34863919'
grep 'chr18' m12-IP4_S4_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA | grep '34863919'
grep 'chr18' normalization_outputs/m12-IP6_w_all_IP4_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA | grep '34863919'
aws s3 cp m12-IP4_vs_IP6_peaks_w_full_read_nums.tsv s3://eclipsebio-m6a051419/data_analysis/
head m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
wc -l m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
ls
rm m12-IP*
clear
ls
rm -r eclipsebio-m6a051419/
clear
ls
logout
exit
source activate eclipsebio
aws s3 cp s3://eclipsebio-m6a051419/normalization_outputs/m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA .
aws s3 cp s3://eclipsebio-m6a051419/normalization_outputs/m12-IP4_S4_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA .
aws s3 cp s3://eclipsebio-m6a051419/normalization_outputs/m12-IP6_S6_w_IP4_significant_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA .
cp scripts/create_summary.py scripts/get_full_read_info_IPvIP_peaks.py
clear
python scripts/get_full_read_info_IPvIP_peaks.py m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP4_S4_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP6_S6_w_IP4_significant_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
head m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
python scripts/get_full_read_info_IPvIP_peaks.py m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP4_S4_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP6_S6_w_IP4_significant_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
head m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
python scripts/get_full_read_info_IPvIP_peaks.py m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP4_S4_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP6_S6_w_IP4_significant_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
python scripts/get_full_read_info_IPvIP_peaks.py m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP4_S4_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP6_S6_w_IP4_significant_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
head m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
python scripts/get_full_read_info_IPvIP_peaks.py m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP4_S4_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP6_S6_w_IP4_significant_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
ls
head m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
python scripts/get_full_read_info_IPvIP_peaks.py m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP4_S4_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP6_S6_w_IP4_significant_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
wc -l m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
python scripts/get_full_read_info_IPvIP_peaks.py m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP4_S4_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP6_S6_w_IP4_significant_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
python scripts/get_full_read_info_IPvIP_peaks.py m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP4_S4_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP6_S6_w_IP4_significant_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
python scripts/get_full_read_info_IPvIP_peaks.py m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP4_S4_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP6_S6_w_IP4_significant_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
python scripts/get_full_read_info_IPvIP_peaks.py m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP4_S4_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP6_S6_w_IP4_significant_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
python scripts/get_full_read_info_IPvIP_peaks.py m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP4_S4_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP6_S6_w_IP4_significant_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
python scripts/get_full_read_info_IPvIP_peaks.py m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP4_S4_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP6_S6_w_IP4_significant_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
clear
rm m12-IP6_S6_w_IP4_significant_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
aws s3 cp s3://eclipsebio-m6a051419/normalization_outputs/m12-IP6_S6_w_IP4_significant_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA /
aws s3 cp s3://eclipsebio-m6a051419/normalization_outputs/m12-IP6_S6_w_IP4_significant_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA .
clear
ls
head m12-IP4_vs_IP6_peaks_w_full_read_nums.tsv 
python scripts/get_full_read_info_IPvIP_peaks.py m12-IP4_vs_IP6_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP4_S4_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA m12-IP6_S6_w_IP4_significant_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA 
head m12-IP4_vs_IP6_peaks_w_full_read_nums.tsv 
grep 'chr18' m12-IP6_S6_w_IP4_significant_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compressed2.bed.full.annotated_proxdist_miRlncRNA | grep '34863919'
grep 'chr18' m12-IP6_S6_w_IP4_significant_peaks_01.basedon.peaks.l2inputnormnew.bed.full.compres | grep '34863472'
mkdir eclipsebio-m6a051419
aws s3 cp s3://eclipsebio-m6a051419/normalization_outputs/m12-IP4_S4_01.basedon.peaks.l2inputnormnew.bed.compressed.bed .
aws s3 cp s3://eclipsebio-m6a051419/m12-IP4_S4.adapterTrim.round2.rmRep.sorted.rmDup.sorted.peaks.bed .
aws s3 cp s3://eclipsebio-m6a051419/m12-IP6_S6.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam eclipsebio-m6a051419/
aws s3 cp s3://eclipsebio-m6a051419/m12-IP6_S6.manifest.txt .
aws s3 cp s3://eclipsebio-m6a051419/m12-inp6_S16.adapterTrim.round2.rmRep.sorted.rmDup.sorted.bam eclipsebio-m6a051419/
cp m12-IP4_S4.adapterTrim.round2.rmRep.sorted.rmDup.sorted.peaks.bed eclipsebio-m6a051419/m12-IP6_S6.adapterTrim.round2.rmRep.sorted.rmDup.sorted.peaks.bed
cp m12-IP6_S6.manifest.txt eclipsebio-m6a051419/m12-IP6_w_all_IP4_peaks.manifest.txt
vim eclipsebio-m6a051419/m12-IP6_w_all_IP4_peaks.manifest.txt 
clear
ls eclipsebio-m6a051419/
clear
ls eclipsebio-m6a051419/
screen
logout
vim scripts/52.34.123.135
clear
ls
ls scripts/
vim scripts/get_full_read_info_IPvIP_peaks.py 
logout
l
cd
ls -a
cp .bourne-apparish Git/Templates/HomeSettings/bourne-apparish
cp .bashrc Git/Templates/HomeSettings/bashrc 
cp .bash_generic Git/Templates/HomeSettings/bash_generic
cd Git/
bm git
cd
to git
cd /
to git
echo $SHELL
via
cd
cd Git/
l
cd Templates/
l
cd HomeSettings/
l
which gcc
rm bourne-apparish 
vi bashrc
l
cp bash_apparix ~/.bash_apparix
cd
vi .bashrc 
l
cd
exit
l
cd
mv apparix-latest.tar.gz local/unpack/
cd local/unpack/
tar xzvf apparix-latest.tar.gz 
cd apparix-11-062/
l
ls
./configure --prefix=${HOME}/local/
make
make install
cd
to git
cd
l
cd Git/
l
bm git
via
to git
exit
exit
sudo yum install apparix
sudo yum install tmux
mkdir Git
cd Git/
cd ll
cd ..
mkdir -p local/unpack
cd local/unpack
git clone https://github.com/micans/bash-utils.git
cd bash-utils/
l
ls -a
cp .bourne-apparish ~
to git
cd
cd Git/
git clone https://github.com/getopt/Templates.git 
cd Templates/
l
ls
cd HomeSettings/
cp vimrc           ~/.vimrc
cp profile         ~/.profile
cp inputrc         ~/.inputrc
cp bashrc          ~/.bashrc
cp bash_profile    ~/.bash_profile
cp bash_path       ~/.bash_path
cp bash_generic    ~/.bash_generic
cp bash_apparix    ~/.bash_apparix
cd
source .bashrc
l
cd
source .bash_profile 
cd
ls -a
sotc
vi .bashrc 
source .bashrc
vi .bourne-apparish 
sotc
a
via
sotcd
sotc
vi .bourne-apparish 
sotc
vi .bash_generic 
l
sotc
tmux
exit
cd local/unpack/
exit
to git
cd /
to git
c
l
cd
cd Git/
l
cd Templates/
c
l
soruce activate eclipse
source activate eclipse
cd
source activate eclipse
vi .bashrc
vi .bash_path 
sotc
source activate eclipseio
source activate eclipsebio
l
exit
source activate pureclip_env
which pureclip
exit
l
vi .bash_path 
sotc
which homer
source activate eclipsebio
exit
