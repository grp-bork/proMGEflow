process get_db_seqs {
    tag "${speci}"
    executor "local"
    cpus 1

    input:
    val(speci)
    val(seqdb_tar)

    output:
    tuple val(speci), path("${speci}.ffn.gz"), emit: sequences

    script:
    """
    tar xvf ${seqdb_tar} ${speci}.genes.ffn.gz
    mv -v ${speci}.genes.ffn.gz ${speci}.ffn.gz
    """
    
}
