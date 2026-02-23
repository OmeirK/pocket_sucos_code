INP_PDB=$1
INP_SDF=$2
#PDB100=/lus/lfs1aip2/scratch/s5h/omeir.s5h/databases/foldseek_pdb100/pdb
PDB100=/lus/lfs1aip2/scratch/s5h/omeir.s5h/databases/plinder/receptor_db/plinder_db

pdb_n=$(basename ${INP_PDB})
aln_out=aln_${pdb_n}.tsv

foldseek easy-search $1 ${PDB100} ${aln_out} tmp/ --format-mode 4 --format-output query,target,qaln,taln,qstart,tstart,qend,tend,evalue,pident,qseq,tseq,u,t -v 1 > log_search_${pdb_n}.log

#python get_d3i_alignment.py -i ${aln_out} > emboss_${pdb_n}.txt
#python get_d3i_alignment.py -i ${aln_out} --tsv-out > d3i_aln_${pdb_n}.tsv
python get_d3i_alignment.py -i ${aln_out} --tsv-out --format-output query,target,qaln,taln,qstart,tstart,qend,tend,evalue,pident,qseq,tseq,u,t > d3i_aln_${pdb_n}.tsv

python3 Py_calc_pocket_qcov.py -i=d3i_aln_${pdb_n}.tsv -pdb=${INP_PDB} -sdf=${INP_SDF} -o=qcov_${pdb_n}.tsv
python3 Py_calc_sucos.py -q=qcov_${pdb_n}.tsv  -sdf=${INP_SDF} -o=sucos_qcov_${pdb_n}.tsv

rm qcov_${pdb_n}.tsv

