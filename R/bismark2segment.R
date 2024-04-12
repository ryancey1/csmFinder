bismark2segment <- function(files, file_type = "regular", CpG_file, tmp_folder = "temp", split_by_chrom = FALSE) {
  options(stringsAsFactors = FALSE)
  reticulate::source_python(
    file.path(system.file(package = "csmFinder"), "bismark2segment.py")
  )

  if (!file.exists(tmp_folder)) {
    dir.create(tmp_folder)
  }

  if (file_type == "regular") {
    CpgReport.gz <- files
    print("reading bismark report file...")
    print(CpgReport.gz)
    file <- basename(CpgReport.gz)
    CpgReport <- file.path(tmp_folder, gsub(".gz", "", file))
    unfold <- R.utils::gunzip(
      filename = CpgReport.gz,
      remove = FALSE,
      destname = CpgReport
    )
    split_bismark_file(CpgReport)
    file.remove(CpgReport)
    print("generating 4-CpG segments:")
    print(" - handle_bismark")
    handle_bismark(tmp_folder)
    file.remove(list.files(tmp_folder, pattern = "bis_$", full.names = TRUE))
    print(" - process_segment")
    process_segment(CpG_file, tmp_folder)
    file.remove(list.files(tmp_folder, pattern = "segment_$", full.names = TRUE))
    print(" - merge_segment")
    merge_segment(tmp_folder)
    file.remove(list.files(tmp_folder, pattern = "filter_$", full.names = TRUE))
    files <- list.files(tmp_folder, pattern = "final_$", full.names = TRUE)
    segment <- lapply(files, function(seg_file) {
      as.data.frame(data.table::fread(seg_file, sep = "\t", header = FALSE, col.names = c("segment", "patterns")))
    })
    segment <- do.call(rbind, segment)
    file.remove(files)
    if (length(dir(tmp_folder) == 0)) {
      unlink(tmp_folder, recursive = TRUE)
    }
    return(segment)
  } else {
    print("reading bismark report file...")
    for (CpgReport.gz in files)
    {
      print(CpgReport.gz)
      file <- basename(CpgReport.gz)
      CpgReport <- file.path(tmp_folder, gsub(".gz", "", file))
      unfold <- R.utils::gunzip(
        filename = CpgReport.gz,
        remove = FALSE,
        destname = CpgReport
      )
      handle_bismark_single(CpgReport)
      file.remove(unfold)
    }
    print("processing 4-CpG segments...")
    process_segment_single(CpG_file, tmp_folder)
    file.remove(list.files(tmp_folder, pattern = "segment_$", full.names = TRUE))
    merge_segment(getwd())
    file.remove(list.files(tmp_folder, pattern = "filter_$", full.names = TRUE))
    segment <- make.bete.mixture.input(path2segments = "./", outdir = "./", split_by_chrom = split_by_chrom)
    file.remove(list.files(tmp_folder, pattern = "beta$", full.names = TRUE))
    file.remove(list.files(tmp_folder, pattern = "final_$", full.names = TRUE))
    if (length(dir(tmp_folder) == 0)) {
      unlink(tmp_folder, recursive = TRUE)
    }
    return(segment)
  }
}

make.bete.mixture.input <- function(path2segments, outdir, split_by_chrom = FALSE) {
  options(stringsAsFactors = F)
  require(stringr, quietly = TRUE, warn.conflicts = FALSE)
  require(data.table, quietly = TRUE, warn.conflicts = FALSE)
  require(reticulate, quietly = TRUE, warn.conflicts = FALSE)
  reticulate::source_python(paste(system.file(package = "csmFinder"), "read_segment_for_beta_mixture.py", sep = "/"))
  read_segment_for_beta_mixture(path2segments, outdir, as.numeric(split_by_chrom))
  if (split_by_chrom == TRUE) {
    segments_info <- list()
    files <- list.files(outdir, pattern = "beta")
    for (file in files) {
      chr <- str_extract(file, "chr[0-9a-zA-Z]{1,2}")
      # print(paste0(outdir,"/",file))
      segments_info[[chr]] <- as.data.frame(fread(paste0(outdir, "/", file)), sep = "\t", header = F)
    }
    return(segments_info)
  } else {
    file <- list.files(outdir, pattern = "beta")
    # print(paste0(outdir,"/",file))
    segments_info <- as.data.frame(fread(paste0(outdir, "/", file)), sep = "\t", header = F)
    colnames(segments_info) <- c("CHR_ID", "START", "END", "LEN", "CELL_COUNT", "CELL_LIST", "CELL_INFO", "CELL_POSI_INFO")
    return(segments_info)
  }
}
