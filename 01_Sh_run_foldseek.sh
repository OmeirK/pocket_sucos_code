INP_PDB=$1
INP_SDF=$2
OUT_DIR=$3
#FOLDSEEK_DB=/lus/lfs1aip2/scratch/s5h/omeir.s5h/databases/foldseek_pdb100/pdb
FOLDSEEK_DB=/lus/lfs1aip2/scratch/s5h/omeir.s5h/databases/plinder/receptor_db/plinder_db

mkdir $OUT_DIR
pdb_n=$(basename ${INP_PDB})
aln_out=${OUT_DIR}/aln_${pdb_n}.tsv

foldseek easy-search $1 ${FOLDSEEK_DB} ${aln_out} tmp/ --format-mode 4 --format-output query,target,qaln,taln,qstart,tstart,qend,tend,evalue,pident,qseq,tseq,qcov,u,t -v 1 > ${OUT_DIR}/log_search_${pdb_n}.log

#python get_d3i_alignment.py -i ${aln_out} > emboss_${pdb_n}.txt
#python get_d3i_alignment.py -i ${aln_out} --tsv-out > d3i_aln_${pdb_n}.tsv
python get_d3i_alignment.py -i ${aln_out} --tsv-out --format-output query,target,qaln,taln,qstart,tstart,qend,tend,evalue,pident,qseq,tseq,qcov,u,t > ${OUT_DIR}/d3i_aln_${pdb_n}.tsv

python3 Py_calc_pocket_qcov.py -i=${OUT_DIR}/d3i_aln_${pdb_n}.tsv -pdb=${INP_PDB} -sdf=${INP_SDF} -o=${OUT_DIR}/qcov_${pdb_n}.tsv
python3 Py_calc_sucos.py -q=${OUT_DIR}/qcov_${pdb_n}.tsv  -sdf=${INP_SDF} -o=${OUT_DIR}/sucos_qcov_${pdb_n}.tsv

#rm ${OUT_DIR}/qcov_${pdb_n}.tsv

