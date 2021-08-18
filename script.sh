#!/usr/bin/bash
#Purpose (1/2): install funannotate 
#Creation date: May 28, 2021
#Revision date: June 28, 2021
#Tutorial: https://github.com/nextgenusfs/funannotate
#Read the docs: https://funannotate.readthedocs.io/en/latest/install.html
#manuscript: https://academic.oup.com/bioinformatics/article/33/18/2936/3861332
# troubleshooting reference: https://github.com/nextgenusfs/funannotate/issues/423
#							 https://github.com/nextgenusfs/funannotate/issues/242 

# This script assumes that FUNGAP.sh and its dependencies have been installed

######################################################################################
#Setting up FunGAP (dependencies for Funannotate included)
######################################################################################
################################################################################################
# We first need to create private python environment & install miniconda installer (Just Once if Miniconda not installed)
# HPC info: https://hpc.nih.gov/apps/python.html#envs 
cd /data/$USER
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
mkdir -p /scratch/$USER/temp
TMPDIR=/scratch/$USER/temp bash Miniconda3-latest-Linux-x86_64.sh -p /data/$USER/conda -b
rm Miniconda3-latest-Linux-x86_64.sh
##################################################################################################

#get readyy to intall with conda
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate base
conda update conda
which conda  # It should be $HOME/anaconda3/condabin/conda
                # it's: /data/proctordm/conda/bin/conda


#### define home and make a fungap directory -- probably should make this direrctory somewhere in $PATH like bin 
HOME=/home/proctordm
mkdir $HOME/FunGAP
cd $HOME/FunGAP


############Use mamba to install dependencies 

# Install Mamba package manager (faster!)
conda install mamba -n base -c conda-forge

# Create FunGAP environment and install dependencies using Mamba
conda create -y -n fungap
conda activate fungap
mamba install \
  braker2=2.1.5 trinity=2.12.0 repeatmodeler=2.0.1 hisat2=2.2.1 pfam_scan=1.6 busco=5.1.2 \
  -c bioconda -c conda-forge

#(editSkyler) I had a bash error w/ mamba command, but it was resolved by running this entire
  #chunk at once (C&P all three comands simultaneously)


#The config/ directory from AUGUSTUS can be accessed with the variable AUGUSTUS_CONFIG_PATH.
#BRAKER2 requires this directory to be in a writable location, so if that is not the case, copy this directory to a writable location, e.g.:
cp -r /data/proctordm/conda/envs/fungap/config/ /home/proctordm/annotation
export AUGUSTUS_CONFIG_PATH=/home/proctordm/annotation

        #Due to license and distribution restrictions, GeneMark, GenomeThreader and ProtHint should be additionally installed for BRAKER2 to fully work.
        #These packages can be either installed as part of the BRAKER2 environment, or the PATH variable should be configured to point to them.
        #The GeneMark key should be located in /home/proctordm/.gm_key and GENEMARK_PATH should include the path to the GeneMark executables.
 
# Install Python and Perl modules (within fungap environment)
pip install biopython bcbio-gff markdown2 matplotlib
cpanm YAML Hash::Merge Logger::Simple Parallel::ForkManager MCE::Mutex Thread::Queue threads

# Install Maker using Mamba (Maker installation is conflict with Busco)
conda deactivate
conda create -y -n maker
conda activate maker
mamba install maker=3.01.03 -c bioconda -c conda-forge


#2. Download and install FunGAP
#2.1. Download FunGAP

#Download FunGAP using GitHub clone. Suppose we are installing FunGAP in your $HOME directory, but you are free to change the location. $FUNGAP_DIR is going to be your FunGAP installation directory.

cd $HOME  # or wherever you want
git clone https://github.com/CompSynBioLab-KoreaUniv/FunGAP.git
export FUNGAP_DIR=$(realpath FunGAP/)
# You can put this export command in the your .bashrc file
# so that you don't need to type every time you run the FunGAP ***************TO DO 


#3. Download Pfam

mkdir -p $FUNGAP_DIR/db/pfam
cd $FUNGAP_DIR/db/pfam
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
gunzip Pfam-A.hmm.gz Pfam-A.hmm.dat.gz
conda activate fungap
hmmpress Pfam-A.hmm  # HMMER package (would be automatically installed in the above Anaconda step)

#4. install gene-mark 
#go here and download it: http://exon.gatech.edu/GeneMark/license_download.cgi
  # Pick: GeneMark-ES/ET/EP ver 4.65_lic LINUX 64 
#transfer the software AND BE SURE TO GET THE KEY!!! to biowulf: $FUNGAP_DIR/external/

mkdir $FUNGAP_DIR/external/
mv gmes_linux_64.tar.gz gm_key_64.gz $FUNGAP_DIR/external/  # Move your downloaded files to this directory
cd $FUNGAP_DIR/external/
tar -zxvf gmes_linux_64.tar.gz
#gunzip gm_key_64.gz #not sure why this was commented out, but I needed to run gunzip to use the next cp
cp gm_key_64 ~/.gm_key
#perl gmes_petap.pl to check proper installation (should show help menu)
 

#4.2. Change the perl path

#GeneMark forces to use /usr/bin/perl instead of conda-installed perl. You can change this by running change_path_in_perl_scripts.pl script.

cd $FUNGAP_DIR/external/gmes_linux_64/
perl change_path_in_perl_scripts.pl "/usr/bin/env perl"
#test path change by assessing shebang (#!) line 
  ## less $FUNGAP_DIR/external/gmes_linux_64/gmes_petap.pl

#4.3 Check GeneMark and its dependencies are correctly installed.

cd $FUNGAP_DIR/external/gmes_linux_64/
perl ./gmes_petap.pl #note you have to include perrl before the script name


#5. Download RepeatMasker databases
conda activate fungap
cd $(dirname $(which RepeatMasker))/../share/RepeatMasker
# ./configure command will download required databases
echo -e "\n2\n$(dirname $(which rmblastn))\n\n5\n" > tmp && ./configure < tmp

# It should look like this
ls $(dirname $(which RepeatMasker))/../share/RepeatMasker/Libraries
# Artefacts.embl  Dfam.hmm       RepeatAnnotationData.pm  RepeatMasker.lib.nin  RepeatPeps.lib      RepeatPeps.lib.psq
# CONS-Dfam_3.0   README.meta    RepeatMasker.lib         RepeatMasker.lib.nsq  RepeatPeps.lib.phr  RepeatPeps.readme
# Dfam.embl       RMRBMeta.embl  RepeatMasker.lib.nhr     RepeatMaskerLib.embl  RepeatPeps.lib.pin  taxonomy.dat


#6. Configure FunGAP

#This script allows users to set and test (by --help command) all the dependencies. If this script runs without any issue, you are ready to run FunGAP!

cd $FUNGAP_DIR
conda activate maker
export MAKER_DIR=$(dirname $(which maker))
echo $MAKER_DIR  # /home/ubuntu/anaconda3/envs/maker/bin
conda activate fungap
python ./set_dependencies.py \
  --pfam_db_path db/pfam/ \
  --genemark_path external/gmes_linux_64/ \
  --maker_path ${MAKER_DIR}



######################################################################################
#Setting up funannotate
######################################################################################
screen
sinteractive --mem=96g --cpus-per-task 8 --time 36:00:00

#install mamba into base environment
conda activate base
conda install -n base mamba

#then use mamba as drop in replacmeent & install funannotate
mamba create -n funannotate funannotate
#start up conda ENV
conda activate funannotate
mamba install \
  funannotate \
  -c bioconda -c conda-forge

#set ENV variable for $FUNANNOTATE_DB
echo "export FUNANNOTATE_DB=/your/path" > /data/sungskh/conda/envs/funannotate/etc/conda/activate.d/funannotate.sh
echo "unset FUNANNOTATE_DB" > /data/sungskh/conda/envs/funannotate/etc/conda/deactivate.d/funannotate.sh

# find which directories are in path
echo $PATH 
#make a copy of the gemes_linux executable in a directory in your path: note may want to change funGap to access this new target_dierctory
cp -r /home/sungskh/FunGAP/external/gmes_linux_64 /home/sungskh/bin
export PATH=/home/sungskh/bin/gmes_linux_64:$PATH 
export GENEMARK_PATH=/home/sungskh/bin/gmes_linux_64 #'funannotate check' runs when both of these paths are set
#Troubleshooting GeneMark
# 1) all dependencies are downloaded: see /your/path/gmes_linux_64/README.GeneMark-ES-suite
# 2) perl path is set: perl change_path_in_perl_scripts.pl "/usr/bin/env perl"
# 3) run 'check_install.bash' and/or perl gmes_petap.pl to see if GM is operating


#check that all modules are installed
funannotate check --show-versions #signalp & eggmapper can be loaded from biowulf

#download/setup databases to a writable/readable location
funannotate setup -d $HOME/funannotate_db -f -w -l 
export FUNANNOTATE_DB=/home/sungskh/funannotate_db #location of databse: /home/sungskh/funannotate_db

#run a test annotation 
funannotate test -t all --cpus 8 #took me ~3 hrs



#####################################################################################
#Running funannotate
#####################################################################################
#!/usr/bin/bash
#Purpose (2/2): go through Funannotate pipeline w/ RNA-seq data Caur007, Caur008, Caur009
#Creation date: July 21, 2021
#Update date: August 2, 2021
#Tutorial: https://funannotate.readthedocs.io/en/latest/prepare.html

set -e 
sinteractive --mem=96g --cpus-per-task 8 --time 36:00:00

source /data/$USER/conda/etc/profile.d/conda.sh
conda activate base
conda activate funannotate
mkdir /data/sungskh/Files4RF
# directory I stored assembly data is in /data/sungskh/Files4RF

#get DNA seq
cp -r /data/sungskh/data/caur_genomes/Caur_3166.Tw_acuzr.Nano.fasta /data/sungskh/Files4RF
cd /data/sungskh/Files4RF 

#for some reason, must reset paths each time (error in: funannotate check --show-versions)
export PATH=/home/sungskh/bin/gmes_linux_64:$PATH
export GENEMARK_PATH=/home/sungskh/bin/gmes_linux_64 
export TRINITYHOME=/data/sungskh/conda/envs/funannotate/opt/trinity-2.8.5 

#re-run to get funannotate species 
funannotate setup -d $HOME/funannotate_db -f -w -l 
export FUNANNOTATE_DB=/home/sungskh/funannotate_db


echo "step 1: prepare assembly"
### Preparing Assembly ###
cd /data/sungskh/Files4RF 
# Clean (clean repetitive contigs, though this step may not be necessary for this .fa)
funannotate clean -i Caur_3166.Tw_acuzr.Nano.fasta -o Caur_3166.Tw_acuzr.Nano_clean.fasta

#sorting/rename FASTA headers
funannotate sort -i Caur_3166.Tw_acuzr.Nano_clean.fasta -o Caur_3166.Tw_acuzr.Nano_clean_sort.fasta

#RepeatMasking Assmebly  
funannotate mask -i Caur_3166.Tw_acuzr.Nano_clean_sort.fasta -o Caur_3166.Tw_acuzr.Nano_clean_sort_mask.fasta


### RNA Alignment (funannotate train) ###
echo "step 2: retrive & train RNA Seq data"
cd /data/Segrelab/data/novaseq-met/novaseq25
cp -r /data/Segrelab/data/novaseq-met/novaseq25/Caur007* /data/sungskh/Files4RF
cp -r /data/Segrelab/data/novaseq-met/novaseq25/Caur008* /data/sungskh/Files4RF
cp -r /data/Segrelab/data/novaseq-met/novaseq25/Caur009* /data/sungskh/Files4RF
cd /data/sungskh/Files4RF
gzip -d Caur007.nohuman.R1.fastq.gz
gzip -d Caur008.nohuman.R1.fastq.gz
gzip -d Caur009.nohuman.R1.fastq.gz
gzip -d Caur007.nohuman.R2.fastq.gz
gzip -d Caur008.nohuman.R2.fastq.gz
gzip -d Caur009.nohuman.R2.fastq.gz
gzip -d Caur007.nohuman.S.fastq.gz
gzip -d Caur008.nohuman.S.fastq.gz
gzip -d Caur009.nohuman.S.fastq.gz

module load signalp
module load eggNOG-mapper
funannotate train -i Caur_3166.Tw_acuzr.Nano_clean_sort_mask.fasta -o fun \
    --left Caur007.nohuman.R1.fastq Caur008.nohuman.R1.fastq Caur009.nohuman.R1.fastq \
    --right Caur007.nohuman.R2.fastq Caur008.nohuman.R2.fastq Caur009.nohuman.R2.fastq \
    --single Caur007.nohuman.S.fastq Caur008.nohuman.S.fastq Caur009.nohuman.S.fastq \
    --jaccard_clip --species "Candida auris" --stranded RF\  #check strandedness info if using RNASeq libraries
    --strain 3166_m3 --cpus 36 
module unload funannotate signalp eggNOG-mapper python/3.7
#keeping these modules loaded somehow conflcits with 'funannotate predict'

### Gene Prediction ### 
echo "step 3: Predict"
funannotate predict -i Caur_3166.Tw_acuzr.Nano_clean_sort_mask.fasta -o fun \
        --species "Candida auris" --strain 3166_m3 --cpus 36 --genemark_mode ES
#training parameters: fun/predict_results/candida_auris_3166_m3.parameters.json
#add species parameters to database
funannotate species -s candida_auris_3166_m3 -a fun/predict_results/candida_auris_3166_m3.parameters.json


### Add UTRs & Update Annotation w/ PASA (funannotate update) ###
echo "step 4: Update"
funannotate update -i fun --cpus 36


### Functional Annotation ###
echo "step 5: Protein Programs"
#Interproscanc
module load interproscan
cp -r /usr/local/apps/interproscan /data/sungskh/Files4RF
cd /data/sungskh/Files4RF
funannotate iprscan -i fun --cpus 36 -m local --iprscan_path /data/sungskh/Files4RF/interproscan/5.42-78.0/interproscan_app
# results are here fun/annotate_misc/iprscan.xml
# this will hang a long time at 0%, but should finish running after a few hours

#Eggnog
module load eggNOG-mapper
cp -r /usr/local/apps/eggNOG-mapper /data/sungskh/Files4RF
#Signalp
module load signalp
#antiSMASH
funannotate remote -i fun -m antismash -e sungskh@nih.gov 
#Phobius
funannotate remote -i fun -m phobius -e sungskh@nih.gov #usually doesn't run for some reason

#Annotation
echo "step 6: Annotate"
funannotate annotate -i fun --cpus 36
#GFF3 data to be visualized in R


rm -rf Caur007* Caur008* Caur009* Caur_3166.Tw_acuzr.Nano.fasta
rm -rf interproscan eggNOG-mapper