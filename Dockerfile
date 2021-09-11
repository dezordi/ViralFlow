FROM dezordi/iam_sarscov2:0.0.5

ENV REFERENCE='' \
    FASTQ1='' \
    FASTQ2='' \
    PREFIXOUT='' \
    THREADS='' \
    DEPTH='' \
    MIN_LEN='' \
    ADAPTERS='' \
    DP_INTRAHOST='' \
    TRIMM_LEN=''

RUN /bin/bash -c "source /home/pango_update"

COPY $REFERENCE $FASTQ1 $FASTQ2 $ADAPTERS /home/

ENTRYPOINT /bin/bash -c "source /home/sars2_assembly_docker $REFERENCE $FASTQ1 $FASTQ2 $PREFIXOUT $THREADS $DEPTH $MIN_LEN $ADAPTERS $DP_INTRAHOST $TRIMM_LEN"