FROM dezordi/iam_sarscov2:0.0.3

ENV REFERENCE='' \
    FASTQ1='' \
    FASTQ2='' \
    PREFIXOUT='' \
    THREADS='' \
    DEPTH='' \
    MIN_LEN='' \
    ADAPTERS='' 

COPY $REFERENCE $FASTQ1 $FASTQ2 $ADAPTERS /home/IAM_SARSCOV2/

RUN /bin/bash -c "source /home/IAM_SARSCOV2/pango_update"

ENTRYPOINT /bin/bash -c "source /home/IAM_SARSCOV2/sars2_assembly /home/IAM_SARSCOV2/$REFERENCE /home/IAM_SARSCOV2/$FASTQ1 /home/IAM_SARSCOV2/$FASTQ2 $PREFIXOUT $THREADS $DEPTH $MIN_LEN /home/IAM_SARSCOV2/$ADAPTERS && cp -r $PREFIXOUT.results/ /store_results"