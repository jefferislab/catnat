#' Read Amira surface (aka HxSurface or HyperSurface) files into hxsurf object
#'
#' @description Read Amira surface (aka HxSurface or HyperSurface) files into hxsurf object. Modified version of similar nat function, nat::read.hxsurf()
#'
#' @param filename file path from which to read
#' @param RegionNames Character vector specifying which regions should be read from file. Default value of NULL => all regions
#' @param RegionChoice Whether the Inner or Outer material, or both (default), should define the material of the patch. See details
#' @param FallbackRegionCol Colour to set regions when no colour is defined
#' @param Verbose Print status messages during parsing when TRUE
#' @param ... additional arguments passed to methods
#'
#' @details Note that when RegionChoice="both" or RegionChoice=c("Inner", "Outer") both polygons in inner and outer regions will be added to named regions. To understand the significance of this, consider two adjacent regions, A and B, with a shared surface. For the polygons in both A and B, Amira will have a patch with (say) InnerRegion A and OuterRegion B. This avoids duplication in the file. However, it might be convenient to add these polygons to both regions when we read them into R, so that regions A and B in our R object are both closed surfaces. To achieve this when RegionChoice="both", read.hxsurf adds these polygons to region B (as well as region A) but swaps the order of the vertices defining the polygon to ensure that the surface directionality is correct. As a rule of thumb, stick with RegionChoice="both". If you get more regions than you wanted, then try switching to RegionChoice="Inner" or RegionChoice="Outer"
#' @return Someneuronlist with cell sidedness in the metadata
#' @export
#' @importFrom utils read.table
get.hxsurf <- function (filename, RegionNames = NULL, RegionChoice = "both",
                          FallbackRegionCol = "grey", Verbose = FALSE, ...)
{
  firstLine = readLines(filename, n = 1)
  if (!any(grep("#\\s+hypersurface\\s+[0-9.]+\\s+ascii", firstLine,
                ignore.case = T, perl = T))) {
    stop(filename, " does not appear to be an Amira HyperSurface ASCII file!")
  }
  initialcaps <- function(x) {
    substr(x, 1, 1) = toupper(substr(x, 1, 1))
    x
  }
  RegionChoice = match.arg(initialcaps(RegionChoice), c("Inner",
                                                        "Outer", "Both"), several.ok = TRUE)
  if (RegionChoice[1] == "Both")
    RegionChoice = c("Inner", "Outer")
  t = readLines(filename)
  nLines = length(t)
  if (Verbose)
    cat(nLines, "lines of text to parse\n")
  dataStart = grep("^\\s*Vertices\\s*", t)[1]
  if (Verbose)
    cat("Data start line =", dataStart, "\n")
  headerLines = t[seq(dataStart - 1)]
  trim = function(x) sub("^\\s+", "", sub("\\s+$", "", x, perl = TRUE),
                         perl = TRUE)
  getfield = function(fName, textLines = headerLines, pos = 2) unlist(strsplit(trim(textLines[grep(fName,
                                                                                                   textLines)]), "\\s+", perl = TRUE))[pos]
  nVertices = as.numeric(getfield("Vertices", t[dataStart],
                                  2))
  if (Verbose)
    cat("nVertices =", nVertices, "\n")
  d = list()
  d$Vertices = read.table(filename, skip = dataStart, nrows = nVertices,
                          col.names = c("X", "Y", "Z"), colClasses = rep("numeric",
                                                                         3))
  d$Regions <- list()
  d$Vertices$PointNo = seq(nrow(d$Vertices))
  if (Verbose)
    cat("Finished processing Vertices\n")
  linesSkipped = dataStart + nVertices - 1
  remainingLines = t[(dataStart + nVertices):nLines]
  PatchDefLine = grep("^\\s*Patches\\s*", remainingLines, perl = TRUE)
  if (Verbose)
    cat("PatchDefLine =", PatchDefLine, "\n")
  nPatches = as.numeric(getfield("Patches", remainingLines[PatchDefLine],
                                 2))
  if (Verbose)
    cat("nPatches =", nPatches, "\n")
  PatchStarts = grep("^\\s*{", remainingLines[PatchDefLine:length(remainingLines)],
                     perl = TRUE) + PatchDefLine - 1
  if (length(PatchStarts) > nPatches)
    PatchStarts = PatchStarts[1:nPatches]
  PatchEnds = grep("^\\s*}", remainingLines[PatchDefLine:length(remainingLines)],
                   perl = TRUE) + PatchDefLine - 1
  if (length(PatchEnds) > nPatches)
    PatchEnds = PatchEnds[1:nPatches]
  TriangleDeflines <- grep("Triangles", remainingLines)
  if (length(TriangleDeflines) != nPatches)
    stop("Incorrect number of Triangle definition lines in",
         filename, "\n")
  for (i in 1:nPatches) {
    if (Verbose)
      cat("TriangleDefline =", TriangleDeflines[i], "\n")
    PatchHeader <- remainingLines[PatchStarts[i]:TriangleDeflines[i]]
    if (Verbose)
      cat("PatchHeader is", length(PatchHeader), "lines long\n")
    for (RegChoice in RegionChoice) {
      RegionName = getfield(paste(RegChoice, "Region",
                                  sep = ""), PatchHeader, 2)
      nTriangles = as.numeric(getfield("Triangles", PatchHeader,
                                       2))
      if (nTriangles < 0 || nTriangles > 1e+05)
        warning("Bad triangle number: ", nTriangles)
      if (Verbose)
        cat("nTriangles =", nTriangles, "for patch =",
            i, "\n")
      if (is.null(RegionNames) || RegionName %in% RegionNames) {
        if (nTriangles == 0 || RegionName == "Exterior")
          next
        thispatch = read.table(filename, skip = linesSkipped +
                                 TriangleDeflines[i], nrows = nTriangles, quote = "",
                               colClasses = "integer", blank.lines.skip = FALSE,
                               fill = FALSE, comment.char = "", col.names = c("V1",
                                                                              "V2", "V3"))
        if (getfield(paste(RegChoice, "Region", sep = ""),
                     PatchHeader, 1) == "OuterRegion") {
          thispatch <- thispatch[, c(1, 3, 2)]
          if (Verbose)
            message("Permuting vertices for ", RegionName,
                    "...")
          colnames(thispatch) <- c("V1", "V2", "V3")
        }
        if (RegionName %in% names(d$Regions)) {
          if (Verbose)
            cat("Adding to patch name", RegionName, "\n")
          d[["Regions"]][[RegionName]] = rbind(d[["Regions"]][[RegionName]],
                                               thispatch)
        }
        else {
          if (Verbose)
            cat("Making new patch name", RegionName,
                "\n")
          d[["Regions"]][[RegionName]] = thispatch
        }
      }
    }
  }
  d$RegionList = names(d$Regions)
  d$RegionColourList <- vector(length = length(d$RegionList))
  closeBraces <- grep("}", headerLines)
  for (regionName in d$RegionList) {
    headerSecStart <- grep(paste0(" ", regionName, " \\{"),
                           headerLines)[1]
    headerSecEnd <- closeBraces[closeBraces > headerSecStart][1]
    colorLine <- grep("Color", headerLines[headerSecStart:headerSecEnd],
                      value = T)
    #if (length(colorLine) > 0) {
#      rgbValues <- strsplit(regmatches(colorLine, gregexpr("[0-9]$|[0-9][^\\.]|[0-9]\\.[0-9]+",
 #                                                          colorLine, perl = T))[[1]], " ")
#      rgbValues = as.character(unlist(strsplit(gsub("[^0-9]", "", unlist(rgbValues)), "")))
 #     color <- rgb(rgbValues[[1]], rgbValues[[2]], rgbValues[[3]])
  #  }
   # else {
      color <- FallbackRegionCol
 #   }
    d$RegionColourList[which(d$RegionList == regionName)] <- color
  }
  class(d) <- c("hxsurf", class(d))
  return(d)
}
