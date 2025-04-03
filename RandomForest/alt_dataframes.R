## trying to make files from brocc_randomforest that can be piped into lefse or microviz

justgrps <- MB_wgroups %>%
  dplyr::select(group) %>%
  rownames_as_column("subject_id") %>%
  dplyr::filter(group %in% c("A", "B"))

newthang <- MB_asvtab %>%
  inner_join(justgrps, by = "subject_id")
newerthang <- as.data.frame(t(newthang))

write.table(newthang, "input_lefse.txt", sep = "\t", 
            row.names = F, quote = F)

write.table(justgrps, "justgrps_metadata.txt", sep = "\t", 
            row.names = F, quote = F)
