#cyrus
cd /data2/Cui/Phragmites/singlecell/paml
./bin/codeml codeml_M0.ctl 
grep omega homologues_PaChr24.43_dnds #omega (dN/dS) =  0.16543
#NULL model
./bin/codeml codeml_branch_site.ctl
grep omega homologues_PaChr24.43_dnds_branch-site.txt #omega = 1.000 fixed
#alternative test
./bin/codeml codeml_branch_site_selection.ctl
#grep omega  homologues_PaChr24.43_dnds_branch-site_altenative.txt
#calculate 2ΔlnL=2×(lnLalt​−lnLnull​), 20.82
#p<0.01


