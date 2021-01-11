#make_personal_txdb("/DATA/projects/j.bhin/reference/GRCh38_gencode_v32_CTAT_lib_Dec062019.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf",
#                   "/home/j.bhin/FGFR/Daniel/R/txdb.GRCh38_genecode_v32_CTAT_lib.sqlite")
make_personal_txdb = function(gtf_file, output){
  # output : "out_path/out_file.rda"
  txdb = GenomicFeatures::makeTxDbFromGFF(file = gtf_file, format = "gtf")
  AnnotationDbi::saveDb(txdb, output)
}
