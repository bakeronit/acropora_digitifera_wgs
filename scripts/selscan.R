
read_ehh_stat <- function(path){
  parse_path_info <- function(path){
    if ( grepl("ihs",path)){
      stat_name="ihs"
      scan_name <- path %>% basename() %>% str_match("[a-z]+")
    } else {
      scan_name <- path %>% basename() %>% str_match("[a-z]+_vs_[a-z]+")
      stat_name <- path %>% basename() %>% str_extract("xp[a-z]{3}")
    }
    # For iHS scan_name will be the population, but for xp* stats it will be a contrasting pair
    c(scan_name,stat_name)
  }
  
  info <- parse_path_info(path)
  read_tsv(path,col_names = c("chr","chr_pos","norm_value"),show_col_types=FALSE) %>% 
    mutate(pval = pnorm(abs(norm_value))) %>% 
    add_column(stat=info[2]) %>% 
    add_column(pop=info[1])
}


read_windowed_xpstats <- function(path){
  scan_name <- path %>% basename() %>% str_match("[a-z]+_vs_[a-z]+") %>% str_split("_vs_") %>% unlist %>% dplyr::first()
  stat_name <- path %>% basename() %>% str_extract("xp[a-z]{3}")
  info <- c(scan_name,stat_name)

  read_tsv(path,col_names = c("chr","start","end","nsnp","fracA","fracB","percentileA","percentileB","max","min")) %>% 
    mutate(pval = pnorm(abs(max))) %>% 
    add_column(stat=info[2]) %>% 
    add_column(pop=info[1])
}

read_windowed_ihs <- function(path) {
  pop <- path %>% basename() %>% str_match("[a-z]+shore") %>% as.character()
  read_tsv(path, col_names = c("chr","start","end","nsnp","frac","percentile","sd")) %>% add_column(pop=pop)
}
