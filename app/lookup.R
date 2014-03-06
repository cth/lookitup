# Handling identifiers for genomic regions such as range expressions, gene names and rs numbers.
library(NCBI2R)
library(memoise)
library(GenomicRanges)
library(shiny)

if(class(try({are.we.inside.shiny.context <- isolate})) == "try-error") {
	isolate <- do.call
}

gene.info <- function(gene) {
	if (is.null(gene)) {
		return(try(a.little.bit.harder))
	}	
	cat(paste("gene.info: ", gene, "\n"))
	lookup.str <- paste0(gene, "[sym]")
	cat(paste("lookup str: ", lookup.str, "\n"))
	id <- try(GetIDs(lookup.str))
	if (class(id) == "try-error") {
		return(id)
	} else {
		cat(paste("got id:", id[1], "\n"))
		return(try(GetGeneInfo(id[1])))
	}
}

snp.info <- function(snp) {
	if (is.null(snp)) {
		return(try(a.little.bit.harder))
	}	
	return(try( GetSNPInfo(snp) ))
}

gene.info.persist <- memoise(gene.info)
snp.info.persist <- memoise(snp.info)

gene.strand <- function(input) {
		if (input$gene %in% labels(test.genes)) {
			test.genes[[input$gene]]$strand
		} else if (!is.null(input$gene.lookup) && input$gene.lookup != 0) {
			gi <- gene.info.persist(input$gene)
			if (class(gi)=="try-error") {
				"NA"
			} else {
				gi$ori
			}
		}
}

gene.description <- function(input) {
		if (input$gene %in% labels(test.genes)) {
			test.genes[[input$gene]]$description
		} else if (!is.null(input$gene.lookup) && input$gene.lookup != 0) {
			gi <- gene.info.persist(input$gene)
			if (class(gi)=="try-error") {
				"NA"
			} else {
				gi$genesummary
			}
		}
}

#Lookup logical interpreter
lookupInterpreter <- function(inputLookup){
    #Function that interpret the user input into three types and always output it as range, summary, name, strand, chr and type
    returnOutput <- data.frame(NA,NA,NA,NA,NA,NA)
    colnames(returnOutput) <- c('range','summary','name','strand','chr','type') 
       #range as chr1-22,X,Y:123-456
    if(grepl(paste0(Reduce(function(...) {paste(...,sep="|") },paste0('chr',1:22)),'chrX|chrY'),inputLookup)){       
        returnOutput$range <- inputLookup
        chrInputted <- strsplit(inputLookup,":")[[1]][1]
        returnOutput$chr <- substring(chrInputted,4,nchar(chrInputted))
        returnOutput$type <- 'Range'
    } else if(grepl('rs[0-9]',inputLookup)){ #snp with rs name rs[0-9]
        si <- isolate({ snp.info.persist(inputLookup) })
        returnOutput$range <- paste0('chr',si$chr,':',si$chrpos,"-",si$chrpos)
        returnOutput$summary <- list(tags$ul(tags$li(paste0('Gene: ',si$genesymbol)),tags$li(paste0('Class: ',si$fxn_class)),tags$li(paste0('\nSpecies: ',si$species))))
        returnOutput$name <- inputLookup
        returnOutput$chr <- si$chr
        returnOutput$type <- 'SNP'
    } else {
            gi <- isolate({gene.info.persist(inputLookup)})
            if (class(gi)=="try-error") {
                print("Input not understood, gene is not found or servers are down")
                returnOutput$name <- inputLookup
                returnOutput$summary <- "Input not understood, not found or servers are down"
                returnOutput$type <- 'try-error'
            }
            else{
                returnOutput$range <- paste0("chr",gi$chr,":",gi$GeneLowPoint,'-',gi$GeneHighPoint)
                returnOutput$name <- inputLookup
                returnOutput$strand <- gi$ori
                returnOutput$summary <- gi$genesummary
                returnOutput$chr <- gi$chr
                returnOutput$type <- 'Gene'
            }
        }
    return(returnOutput)
}

# Parse a range expression in to a Genomic Range
genomicRangeParser <- function(rangeExpr,name="") {
	rangeObj <- NULL
	# Range starts with chr
	if (grepl("chr([0-9]+):([1-9]+[0-9]*)(-([1-9]+[0-9]*))?",rangeExpr)) {
		r <- regexec("chr([0-9]+):([1-9]+[0-9]*)(-([1-9]+[0-9]*))?", rangeExpr)
		matches <- regmatches(rangeExpr,r)

		print(matches)
				
		if (matches[[1]][[5]] == "") { # Singular range
			rangeObj <- GRanges(seqnames = matches[[1]][2],
								ranges = IRanges(start=as.numeric(matches[[1]][3]), end=as.numeric(matches[[1]][3])))
		} else { # Assume full range
			rangeObj <- GRanges(seqnames = matches[[1]][2],
								ranges = IRanges(start=as.numeric(matches[[1]][3]), end=as.numeric(matches[[1]][5])))
		}
		values(rangeObj) <- list(name=name)
	}
	return(rangeObj)
}

rangesFromInputList <- function(inputs,reduce.ranges=F) {
	lookupObjs <- lapply(inputs,lookupInterpreter)
	ranges <- sapply(lookupObjs, function(x) { genomicRangeParser(x$range,x$name) }) 
	oneGR <- Reduce(c,sapply(ranges, function(x) { x }))
	if (reduce.ranges) 
		reduce(oneGR)
	else
		oneGR
}

test.rangesFromInputList <- function() {
	rangesFromInputList(c("ID3","FTO","chr1:2300000-23557960", "chr1:2300001"),F)
}
