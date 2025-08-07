process get_db_seqs {
    tag "${speci}"
    executor "local"
    cpus 1

    input:
    val(speci)
    val(seqdb_tar)

    output:
    tuple val(speci), path("${speci}.ffn.gz"), emit: sequences, optional: true
    tuple val(speci), path("${speci}.GETDBS_DONE"), emit: done_sentinel

    script:
    """
    tar xvf ${seqdb_tar} ${speci}.genes.ffn.gz || :
    if [[ -f ${speci}.genes.ffn.gz ]]; then
        mv -v ${speci}.genes.ffn.gz ${speci}.ffn.gz
    fi
    touch ${speci}.GETDBS_DONE
    """
    
}
