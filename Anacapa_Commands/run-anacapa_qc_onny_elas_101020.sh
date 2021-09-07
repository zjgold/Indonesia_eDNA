# EDIT THESE
BASEDIR="/data/anacapa" # change to folder you want shared into container
CONTAINER="/data/anacapa/anacapa-1.5.0.img" # change to full container .img path
DB="/data/anacapa/Anacapa-New-Master_08152019/Anacapa_db" # change to full path to Anacapa_db
DATA="/data/home/zgold/UCSB_2018/PB052918B_onny_indo_SB" # change to input data folder (default 12S_test_data inside Anacapa_db)
OUT="/data/home/zgold/zgold/Onny_indo/onny_elas_11102020" # change to output data folder

# OPTIONAL
FORWARD="/data/home/zgold/zgold/PSRB/forward_primers_bio_12S_co1.txt"
REVERSE="/data/home/zgold/zgold/PSRB/reverse_primers_bio_12S_co1.txt"
MINLENGTH="/data/home/zgold/zgold/PSRB/metabarcode_loci_min_merge_length_bio_12S_co1.txt"

cd $BASEDIR

# If you need additional folders shared into the container, add additional -B arguments below

time singularity exec -B $BASEDIR $CONTAINER /bin/bash -c "$DB/anacapa_QC_dada2.sh -i $DATA -o $OUT -d $DB -f $FORWARD -r $REVERSE -e $MINLENGTH -a nextera -t MiSeq -l"