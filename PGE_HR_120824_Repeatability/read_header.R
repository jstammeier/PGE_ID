# library(anytime)

## read dat time from dat file header

# fn = "Z:\\Laufende_Projekte\\Projekt_Standard\\databases\\Test_wt_Error\\GL 925 926 MULTI 080719 2\\10 ppb.dat"

extract_timestamp = function(fn) {
  tryCatch(
    {
    ## read first 200 bytes of file
    bytes <- readBin(fn, "raw", 200)
    ##
    # Byte Offset if time stamp in dat is 176; R start counting at 1, thus 177
    ##
    d0 = as.integer(bytes[177]) #nolint
    d1 = as.integer(bytes[178]) #nolint
    d2 = as.integer(bytes[179]) #nolint
    d3 = as.integer(bytes[180]) #nolint

    # dat file is encoded in little endian numbers,
    # smallest digit begins from left,
    # i.e. least significant byte (big endian = begin from right)
    epoch = d0 * 2^0 + d1 * 2^8 + d2 * 2^16 + d3 * 2^24
    #
    epoch_last <<- epoch
    #
    return(epoch)
    #epoch time stamp in Datum + Zeit umwandeln
    message("epoch imported")
          
    },
  error = function(cond) {
    #
    message(cond)
    return(epoch_last)
    message("Last epoch was used")
    }
  )
}
