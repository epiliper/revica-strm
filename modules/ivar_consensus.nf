process IVAR_CONSENSUS {
    tag "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    label 'process_high'
    // container 'quay.io/biocontainers/ivar:1.4--h6b7c446_1'

    input:
    tuple val(meta), val(ref_info), path(bam), path(bai)
    tuple val(meta), val(ref_info), path(ref)

    output:
    tuple val(meta), val(ref_info), path("*.fa"),       optional: true, emit: consensus
    tuple val(meta), val(ref_info), path("*.qual.txt"), optional: true, emit: qual
    tuple val(meta), val(ref_info), path("*.mpileup"),  optional: true, emit: mpileup

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: ''
    """

    if [[ \$(basename "$bam") = "FAILED.sorted.bam" ]]; then
        echo "Skipping $prefix consensus with $ref; failed depth/coverage previously"
        rm *.fa
        exit 0 # shouldn't cause fail if the outputs are optional
    fi

    samtools faidx $ref
    viggo $bam -a -A -Q 0 -f $ref \\
        | ivar \\
            consensus \\
            $args \\
            -p $prefix \\
            -i $prefix

    # get rid of linebreaks except the header line
    awk '/^>/ {printf "%s\\n", \$0; next} {printf "%s", \$0} END {print ""}' ${prefix}.fa > ${prefix}_temp.fa
    # removing leading Ns
    sed '/^>/!s/^N\\+//' ${prefix}_temp.fa > ${prefix}_temp_frontNtrimmed.fa
    # remove trailing Ns
    sed '/^>/!s/N\\+\$//' ${prefix}_temp_frontNtrimmed.fa > ${prefix}.fa
    
    rm ${prefix}_temp.fa
    rm ${prefix}_temp_frontNtrimmed.fa
    """
}
