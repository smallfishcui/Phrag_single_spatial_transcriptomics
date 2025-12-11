cd /scratch/Cui3/Salt/VariantWGS_1/Fst
realSFS EU_invasive_monor.saf.idx -P 24 > EU_invasive_monor.sfs
realSFS EU_source_minor.saf.idx -P 24 > EU_source_minor.sfs
realSFS saf2theta EU_invasive_monor.saf.idx -outname EU_invasive_monor -sfs EU_invasive_monor.sfs
realSFS saf2theta EU_source_minor.saf.idx -outname EU_source_minor -sfs EU_source_minor.sfs
thetaStat do_stat EU_invasive_monor.thetas.idx -win 10000 -step 2000  -outnames EU_invasive_monor.theta.thetasWindow.gz
thetaStat do_stat EU_source_minor.thetas.idx -win 10000 -step 2000  -outnames EU_source_minor.theta.thetasWindow.gz

cut -d")" -f4 EU_source_minor.theta.thetasWindow.gz.pestPG | awk '{ print $1"\t"$2"\t"$8 }' | grep 'PaChr24' > PaChr24_TajimaD.txt
