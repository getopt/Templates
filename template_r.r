library(R.utils)
args = R.utils::commandArgs(asValues=TRUE)

example.file  =    NULL
outfi.file    =    "outfi.txt"

if(!is.null(args[["example"]]))
    example.file = args$example
if(!is.null(args[["outfi"]]))
    outfi.file   = args$outfi

if (is.null(example.file))
    stop("must have --example=<example.file> option")
