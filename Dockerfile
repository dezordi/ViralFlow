FROM dezordi/iam_sarscov2:0.0.5

ENV INPUT_DIR='' \
    REFERENCE='' \
    THREADS='' \
    DEPTH='' \
    MIN_LEN='' \
    ADAPTERS='' \
    DP_INTRAHOST='' \
    TRIMM_LEN=''

#RUN /bin/bash -c "source /home/pango_update"
#RUN /bin/bash "echo DOCKER image creation"
#ENTRYPOINT /bin/bash -c "source /home/sars2_assembly_docker $REFERENCE $FASTQ1 $FASTQ2 $PREFIXOUT $THREADS $DEPTH $MIN_LEN $ADAPTERS $DP_INTRAHOST $TRIMM_LEN"
