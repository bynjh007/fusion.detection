#' @export

# mapping chimeric reads that were excluded for STAR-Fusion to genes and add the information to the original table
breaks_to_gene_table = function(chimeric_before_filter, chimeric_after_filter, target_symbol, chimMultimapNmax_0 = F){

  if(length(readLines(chimeric_after_filter)) != 0){

    chimeric_after = read.table(chimeric_after_filter, header = F, stringsAsFactors = F, sep = "\t", comment.char = "")
    chimeric_before = read.table(chimeric_before_filter, header = F, stringsAsFactors = F, sep = "\t", comment.char = "")

    # selecting the candidate reads
    # reads not used by STAR-Fusion because one or both of the breakpoints (left or right breakpoint) are not matched to exon regions
    # possible cases : exon-exon (one is sense and the other is antisense), exon-intron, exon-intergenic, intron-intron, intron-intergenic
    if(chimMultimapNmax_0 == F){
      l_map = 22
      r_map = 23
    } else {
      l_map = 16
      r_map = 17
    }
    unused = chimeric_before[which(chimeric_before[,l_map] == "." | chimeric_before[,r_map] == "."), ]
    t_after = chimeric_after[c(which(grepl(":sense", chimeric_after[,l_map]) & grepl(":antisense", chimeric_after[,r_map])),
                               which(grepl(":sense", chimeric_after[,r_map]) & grepl(":antisense", chimeric_after[,l_map]))),]

    unused = rbind(unused, t_after)
    # filter the reads with multi-mapping
    unused = unused[unused[,10] %in% names(which(table(chimeric_before[,10]) == 1)),]

    # remove PCR duplicates (take only one of the PCR duplicates) # should be considered differently in SE?
    read_map_token = apply(unused[, c(1:3,12, 4:6,14)], 1, paste, collapse = ":")
    unused = unused[match(names(table(read_map_token)), read_map_token),]

    # filter the fusions with mitochondrial chromosome
    unused = unused %>%
      filter(!(grepl("chrM", V1) | grepl("chrM", V4))) %>%
      filter(V1 %in% paste("chr", c(1:22, "X"), sep = ""), V4 %in% paste("chr", c(1:22, "X"), sep = ""))


    if(nrow(unused)!=0){
      # reads mapping to gene
      unused_map = mapping_reads_to_genes(chimeric_table = unused, trx_info = trx_info, target_symbol = target_symbol)

      # filter the reads both brkpts are mapping to the same target and IGH
      a = unused_map %>% mutate(index = 1:nrow(unused_map)) %>% filter(up.symbol == target_symbol, down.symbol == target_symbol) %>% pull(index)
      b = unused_map %>% mutate(index = 1:nrow(unused_map)) %>% filter(grepl("IGH", up.symbol) | grepl("IGH", down.symbol)) %>% pull(index)
      unused_map = unused_map[setdiff(1:nrow(unused_map), c(a,b)),]

    } else {
      unused_map = NULL
    }

  } else {
    unused_map = NULL
  }

  return(unused_map)
}

# mapping chimeric reads to genes
mapping_reads_to_genes = function(chimeric_table, trx_info, target_symbol){

  # 'target_symbol' is used to give prioritization to the gene

  chimeric_table$V2 = ifelse(chimeric_table$V3 == "-", chimeric_table$V2+1, chimeric_table$V2-1)
  chimeric_table$V5 = ifelse(chimeric_table$V6 == "-", chimeric_table$V5-1, chimeric_table$V5+1)
  #####################################
  # find the genes with the break point
  #####################################
  gr = list(gr_left = GRanges(seqnames = chimeric_table$V1,
                              ranges = IRanges(ifelse(chimeric_table$V3 == "-", chimeric_table$V2, chimeric_table$V2),
                                               ifelse(chimeric_table$V3 == "-", chimeric_table$V2, chimeric_table$V2)),
                                               strand = "*"),

            gr_right = GRanges(seqnames = chimeric_table$V4,
                               ranges = IRanges(ifelse(chimeric_table$V6 == "-", chimeric_table$V5, chimeric_table$V5),
                                                ifelse(chimeric_table$V6 == "-", chimeric_table$V5, chimeric_table$V5)),
                                                strand = "*"))

  gene_target = matrix(NA, nrow(chimeric_table), 2)
  for(i in 1:2){
    # check whether there are two genes (different strand) with the breakpoints
    a = data.frame(trx_info$trx_gene_AA)[findOverlaps(gr[[i]], trx_info$trx, ignore.strand = T, select = "first"), c("GENEID", "aa_length")]
    b = data.frame(trx_info$trx_gene_AA)[findOverlaps(gr[[i]], trx_info$trx, ignore.strand = T, select = "last"), c("GENEID", "aa_length")]
    aa_tab = cbind(a$aa_length, b$aa_length); aa_tab[is.na(aa_tab)] = 0
    gene_tab = cbind(a$GENEID, b$GENEID)

    ind = apply(aa_tab, 1, which.max)
    t_gene_target = c()
    for(j in 1:length(ind)){
      t_gene_target = c(t_gene_target, gene_tab[j, ind[j]])
    }
    gene_target[,i] = t_gene_target
  }

  gene_target_str = cbind(trx_info$trx_gene_AA_major$STRAND[match(gene_target[,1], trx_info$trx_gene_AA_major$GENEID)],
                          trx_info$trx_gene_AA_major$STRAND[match(gene_target[,2], trx_info$trx_gene_AA_major$GENEID)])
  read_target_str = cbind(ifelse(chimeric_table$V3 == gene_target_str[,1], "sense", "antisense"),
                          ifelse(chimeric_table$V6 == gene_target_str[,2], "sense", "antisense"))
  read_target_str[is.na(read_target_str)] = "NA"

  # sense-sense or antisense-antisense is in-strand rearrangement
  ind_sense = which(read_target_str[,1] == "sense" & read_target_str[,2] == "sense")
  ind_antisense = which(read_target_str[,1] == "antisense" & read_target_str[,2] == "antisense")
  ind_incon = setdiff(1:nrow(chimeric_table), c(ind_sense, ind_antisense))

  t_table = chimeric_table[,c(1:7, 10)]
  colnames(t_table) = paste("V", 1:8, sep = "")

  # sense=sense
  t_table_sense = t_table[ind_sense,]

  # swap the alignment info in case both first and second alignments are antisense mapping
  t_table_antisense = t_table[ind_antisense,]
  t_table_antisense = t_table_antisense[, c(4:6,1:3,7,8)]
  if(length(ind_antisense)>0){
    t_table_antisense[,3] = unlist(lapply(t_table_antisense[,3], function(X){setdiff(c("+","-"), X)}))
    t_table_antisense[,6] = unlist(lapply(t_table_antisense[,6], function(X){setdiff(c("+","-"), X)}))
  }
  colnames(t_table_antisense) = paste("V", 1:8, sep = "")

  # inconsistent mapping (can be in-strand for upstream and out-of-strand for downstream or vice versa)
  t_table_incon_1 = t_table[ind_incon,]
  t_table_incon_2 = t_table[ind_incon,]
  t_table_incon_2 = t_table_incon_2[, c(4:6,1:3,7,8)]
  if(length(ind_incon)>0){
    t_table_incon_2[,3] = unlist(lapply(t_table_incon_2[,3], function(X){setdiff(c("+","-"), X)}))
    t_table_incon_2[,6] = unlist(lapply(t_table_incon_2[,6], function(X){setdiff(c("+","-"), X)}))
  }
  colnames(t_table_incon_2) = paste("V", 1:8, sep = "")

  t_table_comb = list(t_table_sense, t_table_antisense, t_table_incon_1, t_table_incon_2)

  gene_target_reorder = list(gene_target[ind_sense,], gene_target[ind_antisense, c(2,1)], gene_target[ind_incon,], gene_target[ind_incon, c(2,1)])
  gr_reorder = list(sense = list(gr_left = GRanges(seqnames = t_table_sense$V1, ranges = IRanges(t_table_sense$V2, t_table_sense$V2), strand = t_table_sense$V3),
                                 gr_right = GRanges(seqnames = t_table_sense$V4, ranges = IRanges(t_table_sense$V5, t_table_sense$V5), strand = t_table_sense$V6)),
                    antisense = list(gr_left = GRanges(seqnames = t_table_antisense$V1, ranges = IRanges(t_table_antisense$V2, t_table_antisense$V2), strand = t_table_antisense$V3),
                                     gr_right = GRanges(seqnames = t_table_antisense$V4, ranges = IRanges(t_table_antisense$V5, t_table_antisense$V5), strand = t_table_antisense$V6)),
                    incons_1 = list(gr_left = GRanges(seqnames = t_table_incon_1$V1, ranges = IRanges(t_table_incon_1$V2, t_table_incon_1$V2), strand = t_table_incon_1$V3),
                                    gr_right = GRanges(seqnames = t_table_incon_1$V4, ranges = IRanges(t_table_incon_1$V5, t_table_incon_1$V5), strand = t_table_incon_1$V6)),
                    incons_2 = list(gr_left = GRanges(seqnames = t_table_incon_2$V1, ranges = IRanges(t_table_incon_2$V2, t_table_incon_2$V2), strand = t_table_incon_2$V3),
                                    gr_right = GRanges(seqnames = t_table_incon_2$V4, ranges = IRanges(t_table_incon_2$V5, t_table_incon_2$V5), strand = t_table_incon_2$V6)))

  # transcripts mapping to the target gene (for forward, reverse of out-of-strand RE)
  loc_in_trx_comb = vector("list", length = 4)
  for(i in 1:4){

    if(nrow(t_table_comb[[i]])>0){
      # for each alignment
      loc_in_trx = c()
      for(j in 1:nrow(t_table_comb[[i]])){

        ##################################################
        # start of RE - transcripts mapping to start read
        ##################################################

        t_range_l = data.frame(trx_info$gene_exon[which(names(trx_info$gene_exon) == matrix(gene_target_reorder[[i]], ncol = 2)[j,1])])
        t_range_l = t_range_l[!grepl("@", t_range_l$group_name),]

        # presence of overlapping gene
        if(nrow(t_range_l)>0){

          # overlapping gene is in "-" strand --> check whether junction is at the end of known exon
          if(t_range_l$strand[1] == "-"){
            t_exon_l = t_range_l$exon_name[which((t_range_l$seqnames %in% t_table_comb[[i]]$V1[j]) & ((t_range_l$start) %in% t_table_comb[[i]]$V2[j]))]
            # overlapping gene is in "+" strand --> check whether junction is at the end of known exon
          } else {
            t_exon_l = t_range_l$exon_name[which((t_range_l$seqnames %in% t_table_comb[[i]]$V1[j]) & ((t_range_l$end) %in% t_table_comb[[i]]$V2[j]))]
          }

          # reference transcript to annotate breakpoint
          t_trx_l = trx_info$trx_gene_AA %>% filter(TXNAME %in% trx_info$exon_info$TXNAME[which(trx_info$exon_info$EXONNAME %in% t_exon_l)]) %>%
            arrange(desc(aa_length), desc(trx_length))
          # junctions from the end of known exon - select the most representative transcript (aa length, trx length)
          # this is for the case whether the junctions from the exon of non-major transcript
          # identify the location of the breakpoint based on the transcript
          if(nrow(t_trx_l)>0){
            t_trx_id_l = t_trx_l$TXID[1]
            splice_type_l = "ref:in-strand"
            # if breakpoint is somewhere in the exon or intron, major transcript will be used as reference
          } else{
            # sorting the transcripts of target gene based on protein and trx length
            t_trx_id_l = (trx_info$trx_gene_AA[queryHits(findOverlaps(trx_info$trx, gr_reorder[[i]][[1]][j])),] %>%
                            arrange(desc(aa_length), desc(trx_length)) %>% pull(TXID))[1]

            # in case of out-of-strand, no trx is mapped
            if(is.na(t_trx_id_l)){
              splice_type_l = "non_ref:out-of-strand"
              t_trx_id_l = (trx_info$trx_gene_AA[queryHits(findOverlaps(trx_info$trx, gr_reorder[[i]][[1]][j], ignore.strand = T)),] %>%
                              arrange(desc(aa_length), desc(trx_length)) %>% pull(TXID))[1]
            } else {
              splice_type_l = "non_ref:in-strand"
            }
          }
          G_ref_l = create_Gref(genome_info = trx_info, TXID = t_trx_id_l)
          loc_in_trx_l = G_ref_l[queryHits(findOverlaps(G_ref_l, gr_reorder[[i]][[1]][j], ignore.strand = T))]

          # if fusion of the downstream occurs at 5' of the gene, brkpt is not overlapped with reference (brkpt + 1/-1)
          # target range is labeled as just first exon of the transcript as we don't know the exact breakpoint
          if(length(loc_in_trx_l)==0){
            loc_in_trx_l = tail(G_ref_l,1)
            if(t_range_l$strand[1] == "+" & (t_table_comb[[i]]$V2[j]) == tail(end(G_ref_l),1)){
              loc_in_trx_r$id = "3P"
            } else if (t_range_l$strand[1] == "-" & (t_table_comb[[i]]$V2[j]) == tail(start(G_ref_l),1)){
              loc_in_trx_l$id = "3P"
            } else {
              stop("undefined breakpoint")
            }
          }

          # no overlapping gene (intergenic)
        } else {
          t_trx_id_l = "NA"
          splice_type_l = "non_ref:intergenic"
          loc_in_trx_l = gr_reorder[[i]][[1]][j]
          loc_in_trx_l$id = "intergenic"
        }

        ##################################################
        # end of RE - transcripts mapping to end read
        ##################################################
        t_range_r = data.frame(trx_info$gene_exon[which(names(trx_info$gene_exon) == matrix(gene_target_reorder[[i]], ncol = 2)[j,2])])
        t_range_r = t_range_r[!grepl("@", t_range_r$group_name),]

        # presence of overlapping gene
        if(nrow(t_range_r)>0){

          # overlapping gene is in "-" strand --> check whether junction is at the end of known exon
          if(t_range_r$strand[1] == "-"){
            t_exon_r = t_range_r$exon_name[which((t_range_r$seqnames %in% t_table_comb[[i]]$V4[j]) & ((t_range_r$end) %in% t_table_comb[[i]]$V5[j]))]
            # overlapping gene is in "+" strand --> check whether junction is at the end of known exon
          } else {
            t_exon_r = t_range_r$exon_name[which((t_range_r$seqnames %in% t_table_comb[[i]]$V4[j]) & ((t_range_r$start) %in% t_table_comb[[i]]$V5[j]))]
          }

          # reference transcript to annotate breakpoint
          t_trx_r = trx_info$trx_gene_AA %>% filter(TXNAME %in% trx_info$exon_info$TXNAME[which(exon_info$EXONNAME %in% t_exon_r)]) %>%
            arrange(desc(aa_length), desc(trx_length))
          # junctions from the end of known exon - select the most representative transcript (aa length, trx length)
          # this is for the case whether the junctions from the exon of non-major transcript
          # identify the location of the breakpoint based on the transcript
          if(nrow(t_trx_r)>0){
            t_trx_id_r = t_trx_r$TXID[1]
            splice_type_r = "ref:in-strand"
            # if breakpoint is somewhere in the exon or intron, major transcript will be used as reference
          } else{
            # sorting the transcripts of target gene based on protein and trx length
            t_trx_id_r = (trx_info$trx_gene_AA[queryHits(findOverlaps(trx_info$trx, gr_reorder[[i]][[2]][j])),] %>%
              arrange(desc(aa_length), desc(trx_length)) %>% pull(TXID))[1]

            # in case of out-of-strand, no trx is mapped
            if(is.na(t_trx_id_r)){
              splice_type_r = "non_ref:out-of-strand"
              t_trx_id_r = (trx_info$trx_gene_AA[queryHits(findOverlaps(trx_info$trx, gr_reorder[[i]][[2]][j], ignore.strand = T)),] %>%
                              arrange(desc(aa_length), desc(trx_length)) %>% pull(TXID))[1]
            } else {
              splice_type_r = "non_ref:in-strand"
            }
          }
          G_ref_r = create_Gref(genome_info = trx_info, TXID = t_trx_id_r)
          loc_in_trx_r = G_ref_r[queryHits(findOverlaps(G_ref_r, gr_reorder[[i]][[2]][j], ignore.strand = T))]

          # if fusion of the downstream occurs at 5' of the gene, brkpt is not overlapped with reference (brkpt + 1/-1)
          # target range is labeled as just first exon of the transcript as we don't know the exact breakpoint
          if(length(loc_in_trx_r)==0){
            loc_in_trx_r = G_ref_r[1]
            if(t_range_r$strand[1] == "+" & (t_table_comb[[i]]$V5[j]) == head(start(G_ref_r),1)){
              loc_in_trx_r$id = "5P"
            } else if (t_range_r$strand[1] == "-" & (t_table_comb[[i]]$V5[j]) == head(end(G_ref_r),1)){
              loc_in_trx_r$id = "5P"
            } else {
              stop("undefined breakpoint")
            }
          }

          # no overlapping gene (intergenic)
        } else {
          t_trx_id_r = "NA"
          splice_type_r = "non_ref:intergenic"
          loc_in_trx_r = gr_reorder[[i]][[2]][j]
          loc_in_trx_r$id = "intergenic"
        }

        loc_in_trx_l$brkpt = paste(t_table_comb[[i]][j, c("V1", "V2", "V3")], collapse = ":")
        loc_in_trx_r$brkpt = paste(t_table_comb[[i]][j, c("V4", "V5", "V6")], collapse = ":")
        loc_in_trx_l$splice_type = splice_type_l
        loc_in_trx_r$splice_type = splice_type_r

        loc_in_trx_l$geneid = matrix(gene_target_reorder[[i]], ncol = 2)[j,1]
        loc_in_trx_r$geneid = matrix(gene_target_reorder[[i]], ncol = 2)[j,2]

        loc_in_trx_l$symbol = ifelse(!is.na(loc_in_trx_l$geneid), trx_info$trx_gene_AA_major %>% filter(GENEID %in% loc_in_trx_l$geneid) %>% pull(SYMBOL), NA)
        loc_in_trx_r$symbol = ifelse(!is.na(loc_in_trx_r$geneid), trx_info$trx_gene_AA_major %>% filter(GENEID == loc_in_trx_r$geneid) %>% pull(SYMBOL), NA)

        loc_in_trx_l$txname = ifelse(!is.na(loc_in_trx_l$geneid), trx_info$trx_gene_AA %>% filter(TXID == t_trx_id_l) %>% pull(TXNAME), NA)
        loc_in_trx_r$txname = ifelse(!is.na(loc_in_trx_r$geneid), trx_info$trx_gene_AA %>% filter(TXID == t_trx_id_r) %>% pull(TXNAME), NA)

        # check whether target transcript is major isoform or not
        loc_in_trx_l$tx_major = loc_in_trx_l$txname %in% trx_info$trx_gene_AA_major$TXNAME
        loc_in_trx_r$tx_major = loc_in_trx_r$txname %in% trx_info$trx_gene_AA_major$TXNAME

        loc_in_trx_l = data.frame(loc_in_trx_l) %>% mutate(target_region = paste(seqnames, start, end, strand, sep = ":")) %>%
          select(brkpt, geneid, symbol, txname, tx_major, target_region, id, splice_type)
        loc_in_trx_r = data.frame(loc_in_trx_r) %>% mutate(target_region = paste(seqnames, start, end, strand, sep = ":")) %>%
          select(brkpt, geneid, symbol, txname, tx_major, target_region, id, splice_type)

        loc_in_trx = rbind(loc_in_trx, data.frame(up = loc_in_trx_l, down = loc_in_trx_r, junc_type = t_table_comb[[i]]$V7[j], read_name = t_table_comb[[i]]$V8[j]))

        print(paste(j, " alignment at ", i, sep = ""))
      }

      loc_in_trx_comb[[i]] = loc_in_trx
    }
  }

  # selecting the alignments for the inconsistent alignments
  ind_a = which(loc_in_trx_comb[[3]]$up.splice_type == "ref:in-strand" | loc_in_trx_comb[[3]]$down.splice_type == "ref:in-strand")
  ind_b = setdiff(which(loc_in_trx_comb[[4]]$up.splice_type == "ref:in-strand" | loc_in_trx_comb[[4]]$down.splice_type == "ref:in-strand"), ind_a)
  ind_up_fgfr2_1 = setdiff(which(loc_in_trx_comb[[3]]$up.symbol == target_symbol), c(ind_a, ind_b))
  ind_up_fgfr2_2 = setdiff(which(loc_in_trx_comb[[4]]$up.symbol == target_symbol), c(ind_a, ind_b, ind_up_fgfr2_1))

  loc_in_trx_incons = rbind(loc_in_trx_comb[[3]][ind_a,], loc_in_trx_comb[[4]][ind_b,], loc_in_trx_comb[[3]][ind_up_fgfr2_1,], loc_in_trx_comb[[4]][ind_up_fgfr2_2,])
  if(!is.null(loc_in_trx_incons)){
    loc_in_trx_incons = loc_in_trx_incons %>% mutate(junc_type = ifelse(junc_type>0, 1, junc_type))
  }


  loc_in_trx_final = rbind(loc_in_trx_comb[[1]], loc_in_trx_comb[[2]], loc_in_trx_incons)

  return(loc_in_trx_final)
}

# summarizing FGFR2 breakpoint
summarizing_fusions = function(chimeric_annot, num_total_reads){

  chimeric_annot_junc = chimeric_annot %>%
    mutate(fus_id = do.call("paste", c(chimeric_annot[, c("up.brkpt", "down.brkpt")], sep = "_"))) %>%
    filter(junc_type != -1) %>% select(-c(read_name))

  fus_tab = table(chimeric_annot_junc$fus_id)
  chimeric_annot_junc_summary = unique(chimeric_annot_junc)
  chimeric_annot_junc_summary = chimeric_annot_junc_summary %>%
    mutate(count = fus_tab[match(chimeric_annot_junc_summary$fus_id, names(fus_tab))], FFPM = (count * 10^6) / num_total_reads) %>%
    arrange(desc(count))

  return(chimeric_annot_junc_summary)
}

# Creating GRange object for the given transcript
create_Gref = function(genome_info, TXID){
  ind = as.numeric(TXID)
  t_exon = genome_info$trx_exon[[ind]]
  t_intron = genome_info$trx_intron[[ind]]
  t_str = unique(as.character(strand(t_exon)))

  if(length(t_intron)!=0){
    if(t_str == "-"){
      t_intron = rev(t_intron)
    }
    df1 = data.frame(data.frame(t_exon)[, c(1,2,3,5)], id = paste("Exon", 1:length(t_exon)))
    df2 = rbind(data.frame(data.frame(t_intron)[, c(1,2,3,5)], id = paste("Intron", 1:length(t_intron))), NA)
    df_comb = (gdata::interleave(df1, df2))
    df_comb = df_comb[-nrow(df_comb),]
    G_ref = makeGRangesFromDataFrame(df_comb)
    G_ref$id = df_comb$id
  } else {
    df1 = data.frame(data.frame(t_exon)[, c(1,2,3,5)], id = paste("Exon", 1:length(t_exon)))
    G_ref = makeGRangesFromDataFrame(df1)
    G_ref$id = df1$id
  }
  return(G_ref)
}
