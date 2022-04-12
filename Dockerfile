FROM dezordi/viralflow:0.0.6.2

ENV INPUT_DIR='' \
    REFERENCE='' \
    THREADS='' \
    DEPTH='' \
    MIN_LEN='' \
    ADAPTERS='' \
    DP_INTRAHOST='' \
    TRIMM_LEN=''

RUN /bin/bash -c "source /home/pango_update" && \
    cd /app/ && \ 
    git clone https://github.com/dezordi/ViralFlow.git && \
    cd ViralFlow/ && \
    /bin/bash -c "source /app/ViralFlow/conda_activate" && \
    pip install -e ./

ENTRYPOINT /bin/bash -c "source /app/ViralFlow/viralflow_docker $FASTA $FASTQ1 $FASTQ2 $PREFIXOUT $THREADS $DEPTH $MIN_LEN $ADAPTERS $DP_INTRAHOST $TRIMM_LEN"
