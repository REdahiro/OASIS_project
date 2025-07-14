for i in 1e-6 1e-4 1e-2 1; do
for j in {1..22}; do
qsub -q short@node2[1-6] -q short@node0[2-9] -q short@node1[0-9] -l m_mem_free=8G -l s_vmem=8G -pe OpenMP 1 -N PRScsx_${i}.${j} -cwd -V \
-o log_dir/${i}.${j}_rmHLA.o -e log_dir/${i}.${j}_rmHLA.e -b y \
"python /work22/home/redahiro/PRScsx/PRScsx.py \
--ref_dir=/work22/home/redahiro/PRScsx/1kg_reference \
--bim_prefix=/work22/home/redahiro/analysis/OASIS/PRS/OASIS_base/plink/WGS_auto_1kg.ref_hg38_rmHLA \
--sst_file=GWAS_sumstats/B2_COVID19_sumstats_for_prs.txt,GWAS_sumstats/B2_COVID19_HGI_EUR_sumstats_for_prs.txt \
--n_gwas=58860,2077779 \
--pop=EAS,EUR \
--phi=${i} \
--chrom=${j} \
--meta=True \
--out_dir=PRScsx_out --out_name=B2_COVID19_rmHLA";
done
done

