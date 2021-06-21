library(docopt)

doc <- "
Usage:
  predictXGB.R --data <FILE> --models <FILE> [--out_csv <FILE>] [--out_rdata <FILE>]

Options:
  --data FILE file with expression data: csv, tsv, txt, rds or RData (containing data.frame named 'expr')
  --models FILE file with list of XGBoost models (.RData)
  --out_csv FILE output csv file
  --out_rdata FILE output RData file
  --help               show this help text"
opt <- docopt(doc)


error.state<-FALSE
if(!file.exists(opt$data)){
  message("File not found: --data ", opt$data)
  error.state<-TRUE
}
if(!file.exists(opt$models)){
  message("File not found: --models ", opt$models)
  error.state<-TRUE
}
if(is.null(opt$out_csv) & is.null(opt$out_rdata)){
  message("Output file not specified: --out_csv or --out_rdata (or both)")
  error.state<-TRUE
}

if(error.state){
  message("Program aborted")
  quit()
}

library(dplyr)
library(readr)  #read_delim
library(reader) #get.delim
library(xgboost)
library(EnvStats) #geoMean


ingestData<-function(fname){
  expr<-NULL
  ext <- tools::file_ext(fname)
  if(ext %in% c("csv", "tsv", "txt")){
    fdelim<-get.delim(fname)
    expr <- read_delim(file = fname, delim = fdelim)
    # removes genes where gene_symbol==NA
    expr<-expr[!is.na(expr[,1]),]
    
    #removes columns with no data (contains NA values only)
    is.empty.col<-apply(expr,2,function(x){sum(is.na(x))==nrow(expr)})
    expr<-expr[,!is.empty.col]
    
    # removes duplicate genenames, keeps the one with highest mean value (i.e. less zero values)
    genenames.dupe<-unique(expr[duplicated(expr[,1]),1])
    expr.keep<-NULL
    for(genename in genenames.dupe){
      df<-expr[expr[,1]==genename,]
      row.mean<-apply(df[,-1],1, mean)
      expr.keep<-rbind(expr.keep, df[row.mean==max(row.mean),])
    }
    expr<-rbind(expr[!expr[,1] %in% genenames.dupe,], expr.keep)
  }
  if(ext=="RData"){
    tmp_env <- new.env()
    load(fname, tmp_env)
    expr<-get(ls(tmp_env), envir=tmp_env) # returns the first object found in tmp_env
    rm(tmp_env)
  }
  if(ext=="rds"){
    expr=data.frame(readRDS(file=fname))
  }
  if(ext %in% c("RData","rds")){ # moves gene names in row.names to first column
    if(sum(is.finite(unlist(expr[,1])))==nrow(expr)){ # unlist required for tibbles
      expr<-data.frame(GENE_SYMBOL=row.names(expr),expr)
      row.names(expr)<-NULL
    }
  }
  if(is.null(expr)){
    message("ingestData: No data has been loaded")
  }else{
    message("ingestData: loaded ",ncol(expr)," samples, ",nrow(expr)," genes")
  }
  
  # until here, data.frame or tibble, with first column "GENE_SYMBOL"
  # final format: data.frame with gene symbols in row.names
  expr<-data.frame(expr)
  rownames(expr)<-expr[,1]
  expr<-expr[,-1]
  
  return(expr)
}

scaleData<-function(data, data.type=NULL){ # valid data.type values: prot, rna
  geomean.center = NULL
  if(data.type=="prot"){
    geomean.center = 1e7
  }
  if(data.type=="rna"){
    geomean.center = 2.5
  }
  data.scaled<-apply(data,2,function(x){
    x*geomean.center/geoMean(x[x!=0])
  })
  return(data.scaled)
}

loadXGB<-function(fname){
  load(file=fname) #list.xgb
  return(list.xgb)
}

inferXGB<-function(expr, list.xgb){
  list.pred<-list()
  no.genes.model<-NULL
  ratio.unmapped<-NULL
  geneset.names<-names(list.xgb)
  for(geneset.name in geneset.names){
    xgb<-list.xgb[[geneset.name]]
    expr.found<-expr[rownames(expr) %in% xgb$feature_names,]
    genes.missing<-xgb$feature_names[!xgb$feature_names %in% rownames(expr)]
    expr.missing<-data.frame(matrix(rep(0,length(genes.missing)*ncol(expr.found)),ncol=ncol(expr.found)))
    rownames(expr.missing)<-genes.missing
    colnames(expr.missing)<-colnames(expr.found)
    expr.new<-rbind(expr.found,expr.missing)
    expr.new<-expr.new[order(rownames(expr.new)),]
    identical(rownames(expr.new),xgb$feature_names)
    pred<-predict(xgb, t(expr.new)) - 10 # models were trained with +10 offset to avoid negative values
    pred<-format(round(pred, 4), nsmall = 4)
    names(pred)<-colnames(expr.new)
    list.pred[[geneset.name]]<-pred
    no.genes.model<-c(no.genes.model,length(xgb$feature_names))
    ratio.unmapped.curr<-length(genes.missing)/length(xgb$feature_names)
    ratio.unmapped.curr<-format(round(ratio.unmapped.curr, 3), nsmall = 3)
    ratio.unmapped<-c(ratio.unmapped,ratio.unmapped.curr)
  }
  pred<-data.frame(bind_rows(list.pred))
  rownames(pred)<-geneset.names
  
  pred<-data.frame(GENESET=row.names(pred),
                   NUMBER_GENES_MODEL=no.genes.model,
                   RATIO_UNMAPPED=ratio.unmapped, 
                   pred)
  row.names(pred)=NULL
  return(pred)
}


expr<-ingestData(fname=opt$data)
expr.norm<-scaleData(data=expr, data.type="prot")
list.xgb<-loadXGB(fname=opt$models)
res.inf<-inferXGB(expr=expr.norm, list.xgb=list.xgb)


if(!is.null(opt$out_csv)){
  message("saving... ", opt$out_csv)
  line.info<-paste0("# PARAMS: predictXGB.R --data ",opt$data, " --models ", opt$models, " --out_csv ", opt$out_csv, " --out_rdata ", opt$out_rdata)
  write.table(line.info, file=opt$out_csv, row.names=FALSE, col.names=FALSE)
  write_csv(res.inf, file=opt$out_csv, col_names=TRUE, append=TRUE)
}

if(!is.null(opt$out_rdata)){
  message("saving... ", opt$out_rdata)
  info=list(data=opt$data, 
            models=opt$models, 
            out_csv=opt$out_csv, 
            out_rdata=opt$out_rdata)
  mapping=res.inf[,c('GENESET', "NUMBER_GENES_MODEL", "RATIO_UNMAPPED")]
  pred=res.inf[,!colnames(res.inf) %in% c("GENESET","NUMBER_GENES_MODEL", "RATIO_UNMAPPED")]
  rownames(pred)<-res.inf$GENESET
  predictXGB<-list(info=info, mapping=mapping, pred=pred)
  save(predictXGB, file=opt$out_rdata, compress=TRUE)
}



message("End of program.")

