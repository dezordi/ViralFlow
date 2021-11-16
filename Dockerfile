FROM dezordi/viralflow:0.0.6

ENV INPUT_DIR='' \
    REFERENCE='' \
    THREADS='' \
    DEPTH='' \
    MIN_LEN='' \
    ADAPTERS='' \
    DP_INTRAHOST='' \
    TRIMM_LEN=''

#RUN /bin/bash -c "source /home/pango_update"
RUN "# activate conda
 __conda_setup="$('/root/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
 echo "Activating conda via bash"
 if [ $? -eq 0 ]; then
     eval "$__conda_setup"
 else
     if [ -f "/root/miniconda3/etc/profile.d/conda.sh" ]; then
         . "/root/miniconda3/etc/profile.d/conda.sh"
     else
         export PATH="/root/miniconda3/bin/:$PATH"
     fi
 fi
 unset __conda_setup
 # <<< conda initialize <<<

 echo "Activating pangolin"
 conda activate pangolin"
#RUN /bin/bash "echo DOCKER image creation"
#ENTRYPOINT /bin/bash -c "source /home/sars2_assembly_docker $REFERENCE $FASTQ1 $FASTQ2 $PREFIXOUT $THREADS $DEPTH $MIN_LEN $ADAPTERS $DP_INTRAHOST $TRIMM_LEN"
