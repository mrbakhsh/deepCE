  #' getCPX
  #' @title Fetch CORUM Complexes from the CORUM Database
  #' @param org Mammalian model organisms including: "Human", "Mouse",
  #' "Pig", "Bovine", "Rat", "Mammalia", "Rabbit", "Dog", "Hamster",
  #' and "MINK".
  #' @param tpath A character string indicating the path to the project
  #' directory. If the directory is
  #' missing, it will be stored in the Temp directory.
  #' @return A list containing protein complexes for mammalian organisms.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @importFrom utils read.delim
  #' @importFrom utils download.file
  #' @importFrom utils unzip
  #' @importFrom dplyr filter
  #' @description This function retrieves CORUM complexes directly from the
  #' CORUM database.
  #' @export

  getCPX <- function(org = "Human",
                     tpath = tempdir()) {

    Organism <- NULL

    url <-
      "http://mips.helmholtz-muenchen.de/corum/download"

    query <-
      paste0(url, "/","coreComplexes.txt.zip")



    temp <- tempfile()

    # Download the file
    download.file(query,temp)


    # Extract the target file from temp
    RawFile <- read.delim(unzip(file.path(temp)))
    file.remove("coreComplexes.txt")

    # Extract the cpx for specific organism
    RawFile <-
      filter(RawFile, Organism %in% org)

    # Extract cpx with >= 2 members
    rawCpx <-
      strsplit(RawFile[, 6], ";")
    names(rawCpx) <-
      RawFile[, 2]
    CORUMcpx <-
      rawCpx[unlist(lapply(rawCpx, length)) >= 3]


    fname <- file.path(tpath, "CORUMcpx.rda")
    save(CORUMcpx, file = fname)

    return(CORUMcpx)
  }


