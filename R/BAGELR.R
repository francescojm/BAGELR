bagelR.computeAllGuidesBFs_v2<-function(cellLine,
                                        NUM_BOOTSTRAPS = 1000,
                                        ESSENTIAL_GENES,
                                        NON_ESSENTIAL_GENES,
                                        inputFolder='../../DATA/R/normalisedCounts_and_FCs/',
                                        outputFolder='../../RESULTS/BAGEL-R_output/',
                                        refGuidesLibrary,diagnosticPlots=TRUE,
                                        whatToTest='newFC'){

  fn<-paste(inputFolder,'/',cellLine,'/_CCRoutput.RData',sep='')
  load(fn)

  print('Computing sgRNA boostrapped Bayesian Factors and classification performances for corrected FCs...')

  correctedFCs<-correctedFCs$corrected_logFCs

  correctedFCs<-cbind(rownames(correctedFCs),correctedFCs)
  colnames(correctedFCs)[1]<-'sgRNA'

  rep1<-bagelR.sgRNA_bootStrapped_BFactors(FOLD_CHANGES = correctedFCs,
                                           NUM_BOOTSTRAPS = NUM_BOOTSTRAPS,
                                           ESSENTIAL_GENES = ESSENTIAL_GENES,
                                           NON_ESSENTIAL_GENES = NON_ESSENTIAL_GENES,
                                           SAMPLE_TO_TEST = whatToTest,
                                           refGuidesLibrary = refGuidesLibrary,
                                           compPerformances = TRUE)
  sgRNA_BFs<-rep1
  print('DONE')

  save(sgRNA_BFs,file=paste(outputFolder,cellLine,'_sgRNAs_BFs.Rdata',sep=''))
  if (diagnosticPlots){
    bagelR.diagnosticPlots(sgRNA_BFs,cellLine = cellLine,outDir = outputFolder)
  }

  return(sgRNA_BFs)
}

bagelR.sgRNA_bootStrapped_BFactors<-function(FOLD_CHANGES,
                                             ESSENTIAL_GENES,
                                             NON_ESSENTIAL_GENES,
                                             NUM_BOOTSTRAPS=1000,
                                             SAMPLE_TO_TEST,
                                             FC_DENS_THRESHOLD=2^-7,
                                             percTrainingSample=80,
                                             refGuidesLibrary,
                                             compPerformances=TRUE){

  rownames(FOLD_CHANGES)<-FOLD_CHANGES$sgRNA

  nguides<-nrow(FOLD_CHANGES)

  genes<-unique(FOLD_CHANGES$gene)
  guides<-FOLD_CHANGES$sgRNA

  ngenes<-length(genes)

  print(paste('n. of sgRNA =',nguides))
  print(paste('n. of unique genes =',ngenes))

  print(paste('n. of considered essential genes = ',length(intersect(ESSENTIAL_GENES,FOLD_CHANGES$gene)),
              ' (out of ',length(ESSENTIAL_GENES),' in reference set)',sep=''))
  print(paste('          targeted by a total number of ',length(which(is.element(FOLD_CHANGES$gene,ESSENTIAL_GENES))),
              ' sgRNAs',sep=''))


  print(paste('n. of considered non essential genes = ',length(intersect(NON_ESSENTIAL_GENES,FOLD_CHANGES$gene)),
              ' (out of ',length(NON_ESSENTIAL_GENES),' in reference set)',sep=''))
  print(paste('          targeted by a total number of ',length(which(is.element(FOLD_CHANGES$gene,NON_ESSENTIAL_GENES))),
              ' sgRNAs',sep=''))

  BFs_across_loops<-matrix(-Inf,nguides,NUM_BOOTSTRAPS,dimnames = list(guides,1:NUM_BOOTSTRAPS))
  BFs_across_loops_including_training<-matrix(-Inf,nguides,NUM_BOOTSTRAPS,dimnames = list(guides,1:NUM_BOOTSTRAPS))

  EssentialGuides<-as.character(FOLD_CHANGES$sgRNA[is.element(FOLD_CHANGES$gene,ESSENTIAL_GENES)])
  nonEssentialGuides<-as.character(FOLD_CHANGES$sgRNA[is.element(FOLD_CHANGES$gene,NON_ESSENTIAL_GENES)])

  nEssGuides<-length(EssentialGuides)
  nNonEssGuides<-length(nonEssentialGuides)

  print('Bootstrap iterations in progress...')

  pb <- txtProgressBar(min=1,max=NUM_BOOTSTRAPS,style=3)

  for (loop in 1:NUM_BOOTSTRAPS){
    setTxtProgressBar(pb, loop)
    guide_train_ess<-EssentialGuides[sample(nEssGuides)[1:ceiling(nEssGuides*percTrainingSample/100)]]
    guide_train_non_ess<-nonEssentialGuides[sample(nNonEssGuides)[1:ceiling(nNonEssGuides*percTrainingSample/100)]]

    guide_training<- union(guide_train_ess,guide_train_non_ess)
    guide_test <- setdiff(guides,guide_training)

    ess_train_fc <- as.numeric(FOLD_CHANGES[guide_train_ess,SAMPLE_TO_TEST])
    non_ess_train_fc <- as.numeric(FOLD_CHANGES[guide_train_non_ess,SAMPLE_TO_TEST])

    kess<-density(ess_train_fc, kernel = "gaussian")
    knon<-density(non_ess_train_fc, kernel = "gaussian")

    x <- seq(-10,2,0.01)
    nonfitx <- approx(knon$x,knon$y,x)$y

    f <- which(nonfitx > FC_DENS_THRESHOLD)
    xmin <- bagelR.round_to_hundredth( min(x[f]) )

    subx <- seq(xmin,max(x[f]),0.01)
    logratio_sample <- log2( approx(kess$x,kess$y,subx)$y / approx(knon$x,knon$y,subx)$y )

    f <- which(logratio_sample == min(logratio_sample,na.rm = TRUE))
    xmax <- bagelR.round_to_hundredth(subx[f])

    RANGE<-seq(xmin,xmax+0.01,0.01)
    RANGE<-round(RANGE,digits = 2)

    logratio_lookup <- log2(approx(kess$x,kess$y,RANGE)$y / approx(knon$x,knon$y,RANGE)$y)
    names(logratio_lookup)<-RANGE

    foldchanges <- FOLD_CHANGES[guide_test,SAMPLE_TO_TEST]

    foldchanges[foldchanges>xmax]<-xmax
    foldchanges[foldchanges<xmin]<-xmin
    foldchanges<-round(foldchanges,digits=2)

    currentLoop_bf<-rep(NA,length(guide_test))

    currentLoop_bf<-logratio_lookup[as.character(foldchanges)]

    names(currentLoop_bf)<-guide_test
    BFs_across_loops[names(currentLoop_bf),loop]<-currentLoop_bf
    BFs_across_loops_including_training[names(currentLoop_bf),loop]<-currentLoop_bf

    foldchanges <- FOLD_CHANGES[guide_training,SAMPLE_TO_TEST]

    foldchanges[foldchanges>xmax]<-xmax
    foldchanges[foldchanges<xmin]<-xmin
    foldchanges<-round(foldchanges,digits=2)

    currentLoop_bf<-rep(NA,length(guide_test))

    currentLoop_bf<-logratio_lookup[as.character(foldchanges)]

    names(currentLoop_bf)<-guide_training
    BFs_across_loops_including_training[names(currentLoop_bf),loop]<-currentLoop_bf

  }

  print('')
  print('DONE')

  Sys.sleep(1)
  close(pb)

  if(compPerformances){
    PERFORMANCES<-
      bagelR.sgRNA_bootStrapped_performances(BS_BF = BFs_across_loops,
                                             refGuidesLibrary = refGuidesLibrary,
                                             ESSENTIAL_GENES = ESSENTIAL_GENES,
                                             NON_ESSENTIAL_GENES = NON_ESSENTIAL_GENES)
  }else{
    PERFORMANCES<-NULL
  }

  BFs_across_loops[BFs_across_loops==-Inf]<-NA

  sgRNA_BFs<-apply(BFs_across_loops,MARGIN = 1,'mean',na.rm=TRUE)
  sgRNA_BFs_sd<-apply(BFs_across_loops,MARGIN = 1,'sd',na.rm=TRUE)

  sgRNA_BFs<-cbind(sgRNA_BFs,sgRNA_BFs_sd)
  colnames(sgRNA_BFs)<-c('avg_bootstr_BFs','sd_bootstr_BFs')


  sgRNA_BFs_inclTr<-apply(BFs_across_loops_including_training,MARGIN = 1,'mean',na.rm=TRUE)
  sgRNA_BFs_sd_inclTr<-apply(BFs_across_loops_including_training,MARGIN = 1,'sd',na.rm=TRUE)

  sgRNA_BFs_inclTr<-cbind(sgRNA_BFs_inclTr,sgRNA_BFs_sd_inclTr)
  colnames(sgRNA_BFs_inclTr)<-c('avg_bootstr_BFs','sd_bootstr_BFs')

  return(list(sgRNA_BFs=sgRNA_BFs,sgRNA_BFs_inclTr=sgRNA_BFs_inclTr,boostPERF=PERFORMANCES))
}

bagelR.sgRNA_bootStrapped_performances<-function(BS_BF,ESSENTIAL_GENES,
                                                 NON_ESSENTIAL_GENES,refGuidesLibrary){


  print('Computing Bootstrapped performances across iterations...')

  essentialGuides<-
    refGuidesLibrary[match(intersect(ESSENTIAL_GENES,refGuidesLibrary[,2]),refGuidesLibrary[,2]),1]

  nonessentialGuides<-
    refGuidesLibrary[match(intersect(NON_ESSENTIAL_GENES,refGuidesLibrary[,2]),refGuidesLibrary[,2]),1]

  tmp<-BS_BF[intersect(rownames(BS_BF),union(essentialGuides,nonessentialGuides)),]

  minimum<-min(tmp[tmp>-Inf])
  maximum<-max(tmp[tmp>-Inf])

  nboostLoops<-ncol(BS_BF)

  range<-seq(minimum,maximum,(maximum-minimum)/999)

  SENSITIVITY<-matrix(NA,1000,nboostLoops)
  SPECIFICITY<-matrix(NA,1000,nboostLoops)
  THRESHOLD<-matrix(NA,1000,nboostLoops)
  PPV<-matrix(NA,1000,nboostLoops)
  AUC<-rep(NA,nboostLoops)

  pb <- txtProgressBar(min=1,max=nboostLoops,style=3)

  for (i in 1:nboostLoops){
    setTxtProgressBar(pb, i)
    currentTestSet<-names(which(BS_BF[,i]!='-Inf'))

    testEssential<-intersect(essentialGuides,currentTestSet)
    testNonessential<-intersect(nonessentialGuides,currentTestSet)

    currentTestSet<-union(testEssential,testNonessential)

    predictions<-BS_BF[currentTestSet,i]

    essentiality<-is.element(currentTestSet,testEssential)+0

    ROC<-roc(essentiality,predictions)

    AUC[i]<-ROC$auc

    COORDS<-coords(ROC,x = range,ret = c('sensitivity','specificity','ppv'))

    SENSITIVITY[,i]<-COORDS['sensitivity',]
    SPECIFICITY[,i]<-COORDS['specificity',]
    PPV[,i]<-COORDS['ppv',]
  }

  THRESHOLD<-range

  avgSens<-apply(SENSITIVITY,MARGIN = 1,'mean',na.rm=TRUE)
  avgSpec<-apply(SPECIFICITY,MARGIN = 1,'mean',na.rm=TRUE)
  avgPpv<-apply(PPV,MARGIN=1,'mean',na.rm=TRUE)

  sdSens<-apply(SENSITIVITY,MARGIN = 1,'sd',na.rm=TRUE)
  sdSpec<-apply(SPECIFICITY,MARGIN = 1,'sd',na.rm=TRUE)
  sdPpv<-apply(PPV,MARGIN=1,'sd',na.rm=TRUE)

  print('')
  print('DONE')
  Sys.sleep(1)
  close(pb)


  return(list(th=THRESHOLD,
              ppv=avgPpv,sens=avgSens,spec=avgSpec,
              sd_ppv=sdPpv,sd_sens=sdSens,sd_spec=sdSpec,AUROC=mean(AUC)))

}

bagelR.diagnosticPlots<-function(sgRNA_BFs,cellLine,outDir){

  fn<-paste(outDir,cellLine,'.pdf',sep='')
  pdf(fn)

  plot(sgRNA_BFs$boostPERF$th,frame.plot = FALSE,
       sgRNA_BFs$boostPERF$ppv,type='l',lwd=2,
       xlab='sgRNA Bayesian Factor',col='red',
       main=paste(cellLine,': avg precision/sensitivity across boostrap iterations',sep=''),
       ylab='')

  par(new=TRUE)
  par(mar=c(4,4,4,5))
  plot(sgRNA_BFs$boostPERF$th,frame.plot = FALSE,
       sgRNA_BFs$boostPERF$sens,type='l',lwd=2,
       xaxt='n',yaxt='n',ylab='',xlab='',col='blue')
  axis(4)

  plot(1-sgRNA_BFs$boostPERF$spec,
       sgRNA_BFs$boostPERF$sens,
       type='l',lwd=2,
       xlab='FPR',col='blue',
       main=paste(cellLine,': avg FPR/TPR across boostrap iterations (avg AUC ',
                  format(sgRNA_BFs$boostPERF$AUROC,digits=4),')',sep=''),
       ylab='TPR')
  lines(x = c(0,1),y=c(0,1))
  dev.off()

}

bagelR.round_to_hundredth<-function(x){
  return (round(x*100) / 100.0)
}

bagelR.geneLevel_BFs<-function(allguidesBFs,refGuidesLibrary){
  nreplicates<-length(allguidesBFs)
  allBFsmatrix<-matrix(NA,nrow(allguidesBFs$sgRNA_BFs),nreplicates,
                       dimnames = list(rownames(allguidesBFs$sgRNA_BFs),names(allguidesBFs)))

  allBFsmatrix[,1]<-allguidesBFs$sgRNA_BFs_inclTr[,1]

  #red_refGuidesLibrary<-
  #    refGuidesLibrary[which(is.element(refGuidesLibrary$sgRNA,rownames(allBFsmatrix))),]

  genes<-unique(refGuidesLibrary[,2])
  ngenes<-length(genes)

  GENElevel_BF<-rep(NA,length(genes))
  GENElevel_BF_sd<-rep(NA,length(genes))
  GENElevel_BF<-cbind(GENElevel_BF,GENElevel_BF_sd)
  rownames(GENElevel_BF)<-genes
  colnames(GENElevel_BF)<-c('Avg','SD')

  print('Computing gene-level Bayes Factors...')
  pb <- txtProgressBar(min=1,max=ngenes,style=3)

  for (i in 1:ngenes){
    setTxtProgressBar(pb, i)
    #guidesetids<-intersect(rownames(allBFsmatrix),as.character(KY_Library_v1.0_list[[genes[i]]]))

    guidesetids<-intersect(rownames(allBFsmatrix),
                           refGuidesLibrary[which(refGuidesLibrary[2]==genes[i]),1])


    GENElevel_BF[i,1]<-mean(c(allBFsmatrix[guidesetids,1]))
    GENElevel_BF[i,2]<-sd(c(allBFsmatrix[guidesetids,1]))
  }

  print('')
  print('DONE')

  Sys.sleep(1)
  close(pb)

  return(GENElevel_BF)
}
