network_dir <- "networks/"
check_minimal <- FALSE

###### DO NOT MODIFY BELOW THIS LINE ######
library(dplyr) 

md5_sums <- c("gr_network_human_21122021.rds", "7b881b2d3f41595409c4b99ca515b171",
              "gr_network_mouse_21122021.rds", "a189561341dc5514774abdb31b2ddcc3",
              "ligand_target_matrix_nsga2r_final.rds", "b09606b04b2d4490418d9028c0e58b9f",
              "ligand_target_matrix_nsga2r_final_mouse.rds","ac80d846fe0bfc4879a5b52ca85ffeb9",
              "ligand_tf_matrix_nsga2r_final.rds","d8951b123e4f7ec6bd5fb50069050972",
              "ligand_tf_matrix_nsga2r_final_mouse.rds","444056f809a6cd9c49f22dc34a517b36",
              "lr_network_human_21122021.rds", "2a155f81e9ffffd5d5e709fe66bcc465",
              "lr_network_mouse_21122021.rds", "cf33ee8b6bf84bdf2d11cab9c8f94b9e",
              "signaling_network_human_21122021.rds", "b294c74b6ff6b959f5b75e50b39ae4c8",
              "signaling_network_mouse_21122021.rds", "bbfa9deb6529d79311f0ceab33c5118c",
              "weighted_networks_nsga2r_final.rds", "80014dee22df42e98d1c608731c685b5",
              "weighted_networks_nsga2r_final_mouse.rds", "6620b465a566a9e556ee0f49f5f19cfa")
md5_sums <- md5_sums[c(F, T)] %>% setNames(md5_sums[c(T, F)])

grep_pattern <- ifelse(check_minimal, "^[lw].*mouse.*rds", ".*rds")
files_to_check <- list.files(network_dir, grep_pattern)
expected_files <- grep(grep_pattern, names(md5_sums), value=TRUE)

if (length(files_to_check) != length(expected_files)){
  missing_networks <- setdiff(expected_files, files_to_check)
  stop(paste0("Missing network(s): ", paste0(missing_networks, collapse=", ")), call. = FALSE)
}

message(paste0("Checking MD5 sum of ", length(files_to_check), " files:\n",
        paste0(files_to_check, collapse="\n")))

for (file in files_to_check){
  if (tools::md5sum(paste0(network_dir, file)) != md5_sums[file]){
    stop(paste0("MD5 sum of ", file, " does not match the expected value."), call. = FALSE)
  }
}

message("All MD5 sums are correct.")


