    #######################################################################################
    #                                                                                     #
    #                  LSx v1.1: An Rscript for reducing lineage evolutionary             #
    #            rate heterogeneity in a gene-by-gene, criterion guided manner            #
    #                                                                                     #
    #######################################################################################



#    Copyright, Carlos J. RIVERA-RIVERA and Juan I. MONTOYA-BURGOS, December 2017 onwards
#
#            LSx is free software: you can redistribute it and/or modify
#            it under the terms of the GNU General Public License as published by
#            the Free Software Foundation, either version 3 of the License, or
#            (at your option) any later version.
#
#            This program is distributed in the hope that it will be useful,
#            but WITHOUT ANY WARRANTY; without even the implied warranty of
#            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#            GNU General Public License for more details.
#
#            You should have received a copy of the GNU General Public License
#            along with this program.  If not, see <http://www.gnu.org/licenses/>.


#              If you use this code for your publications, please cite:

#              Rivera-Rivera and Montoya-Burgos (2017). LSX: Automated
#                  reduction of gene-specific lineage evolutionary rate
#                  heterogeneity for muti-gene phylogeny inference.
#                  bioRxiv. DOI: https://doi.org/10.1101/220053

#              Rivera-Rivera and Montoya-Burgos (2016). LS³: A Method
#                  for Improving Phylogenomic Inferences When Evolutionary
#                  Rates Are Heterogeneous among Taxa. Mol. Biol. Evol.
#                  33(6):1625-34. DOI: https://doi.org/10.1093/molbev/msw043


#                                      Thanks!


checkInputs<-function() {  #Checks that all the inputs are correct
  if (listFile %in% fileList) {
    print(paste("listFile ","FOUND: ",listFile,sep=""),quote=FALSE)
    geneList<<-read.csv(file=listFile,header=FALSE,colClasses="character") #Import the list of genes to be analyzed. below a command to give column headers
  } else {
    print(paste("Error: ",listFile," was not found. Aborting.",sep=""),quote=FALSE)
    quit()
  }

  if (guideTree %in% fileList) {
    print(paste("guideTree ","FOUND: ",guideTree,sep=""),quote=FALSE)
    guideTree<<-read.tree(file=guideTree) #Import guide tree
  } else {
    print(paste("Error: guideTree file ",guideTree," was not found. Aborting.",sep=""),quote=FALSE)
    quit()
  }

  if (CladeSpeciesFile %in% fileList) {
    print(paste("CladeSpeciesFile ","FOUND: ",CladeSpeciesFile,sep=""),quote=FALSE)
    generalSppCl<<-read.csv(file=CladeSpeciesFile,header=FALSE) #Import table with species and clade classification
  } else {
    print(paste("Error: clades/species file ",CladeSpeciesFile," was not found. Aborting.",sep=""),quote=FALSE)
    quit()
  }

  if (GenericPamlCtrlFile %in% fileList) {
    print(paste("GenericPamlCtrlFile ","FOUND: ",GenericPamlCtrlFile,sep=""),quote=FALSE)
    pamlCtrlFile<<-readLines(GenericPamlCtrlFile,-1) # read in generic baseml control file
  } else {
    print(paste("Error: the PAML control file ",GenericPamlCtrlFile," was not found. Aborting",sep=""),quote=FALSE)
    quit()
  }

  if (dataType == "NUC" || dataType == "AA") {
    print(paste("Datatype: ",dataType,sep=""),quote=FALSE)
     if (dataType == "NUC" ){
      basemlCtrlFile<<-pamlCtrlFile
     } else if (dataType == "AA" ) {      #For amino acid data
      if (length(grep("defaultAAMatrix",ls(.GlobalEnv)))) {  #Check if a default aa matrix was given
        if (defaultAAMatrix %in% fileList) {        #If given, check if the file is in the current folder
          print(paste("defaultAAMatrix ","FOUND: ",defaultAAMatrix,sep=""),quote=FALSE)
        } else {                                   #If given and not found, abort
          print(paste("Error: amino acid rate matrix ",defaultAAMatrix," was not found. Aborting.",sep=""),quote=FALSE)
          quit()
        }
      } else {                                     #If not given (i.e. variable was commented out), announce it.
        print("No default AA matrix file given, expecting to see it specified in the input gene list.",quote=FALSE)
      }
      aamlCtrlFile<<-pamlCtrlFile                   #In any case, if it's AA data, reads in the aaml control file
    }
  } else {
    print(paste("Datatype ",dataType," incorrect. Aborting.",sep=""),quote=FALSE)
    quit()
  }

  PamlExecutablePath<<-Sys.which(PamlExecutablePath)

  if (nchar(PamlExecutablePath)) {
    print(paste("PAML executable FOUND: ",PamlExecutablePath,sep=""),quote=FALSE)
    if (dataType == "NUC" && length(grep("baseml",PamlExecutablePath))) {
      print("Baseml executable matches with NUC dataType",quote=FALSE)
    } else if (dataType == "AA" && length(grep("codeml",PamlExecutablePath))) {
      print("Codeml executable matches with AA dataType",quote=FALSE)
    } else {
      print(paste("Error: dataType ",dataType," and PAML executable ",PamlExecutablePath," not compatible. Are you using baseml for NUC and codeml for AA? Aborting.",sep=""),quote=FALSE)
      quit()
    }
  } else {
    print("Error: PAML executable not found in address provided. Aborting.",quote=FALSE)
    if (PlatformOS == "windows") {
      print("Remember that the PAML executable path must be full, and each backslash must be double '\\'",quote=FALSE)
    }
    quit()
  }

  if (Flavor == "LS3" || Flavor == "LS4") {
    print(paste("Algorithm flavor: ",Flavor,sep=""),quote=FALSE)
  } else {
    print(paste("Error: Flavor ",Flavor," incorrect. Aborting.",sep=""),quote=FALSE)
    quit()
  }

  oldPThresh<-pThresh
  pThresh<<-as.numeric(pThresh)
  if(!is.na(pThresh)) {
    if (pThresh > 0 && pThresh <= 1) {
      print(paste("p-Value Threshold: ",pThresh,sep=""),quote=FALSE)
    } else {
      print(paste("Error: p-Value ",pThresh," out-of-bounds or incorrect. Aborting.",sep=""),quote=FALSE)
      quit()
    }
  } else {
    print(paste("Error: pThresh value ",oldPThresh," is not a number. Aborting.",sep=""),quote=FALSE)
    quit()
  }

  oldMinTaxa<-minTaxa
  minTaxa<<-as.integer(minTaxa)
  if (!is.na(minTaxa)) {
    if (minTaxa > 0) {
      print(paste("minTaxa Threshold: ",minTaxa," species/lineage.",sep=""),quote=FALSE)
    } else {
      print(paste("Error: minTaxa value ",minTaxa," out-of-bounds. Aborting.",sep=""),quote=FALSE)
      quit()
    }
  } else {
    print(paste("Error: minTaxa value ",oldMinTaxa," is not a number. Aborting.",sep=""),quote=FALSE)
    quit()
  }

  if (length(grep("numCores",ls(.GlobalEnv)))) {  #Check if multi-threading was requested
    numCores<<-as.numeric(numCores)
    if (numCores > 1) {
      parallelization<<-TRUE
      suppressPackageStartupMessages(library(parallel))
      print(paste("Multi-threading requested on ",numCores," cores.",sep=""),quote=FALSE)
      if (numCores > detectCores()) {
        numCores<<-detectCores()
        print(paste("Multi-threading on the ",numCores," detecteble cores.",sep=""),quote=FALSE)
      }
      cluster<<-makeCluster(numCores)
    } else if (numCores == 1) {
      parallelization<<-FALSE
      print("Number of requested cores is 1, running a sequential analysis.",quote=FALSE)
    } else {
      print("numCores has an invalid entry. Aborting.",quote=FALSE)
      quit()
    }
  } else {
    parallelization<<-FALSE
    print("No multi-threading was requested, running a sequential analysis.",quote=FALSE)
  }

} #End of checkInputs() function



readInData<-function() {  #Reads in the data from the input files
  colnames(generalSppCl)<<-c("Lineage","Species")
  lineagesOfInt<<-levels(generalSppCl$Lineage)[-length(levels(generalSppCl$Lineage))]
  degreesOfFreedom<<-length(lineagesOfInt)-1    #Determines the amount of degrees of freedom for the LRT (the amount of declared lineages - 1)

  if (length(lineagesOfInt) > 4) {
    print("I must warn you that running this algorithm for more than 4 lineages of interest is not recommended. I'll give you 5 seconds to think about this.", quote=FALSE)
    Sys.sleep(5)
  } else {
    print(paste(length(lineagesOfInt)," lineages (df=",degreesOfFreedom," for LRT)",sep=""),quote=FALSE)
  }

  if (is.rooted(guideTree)) {
    print("Input tree is rooted. Proceeding",quote=FALSE)
  } else {
    print("Input tree is not rooted. Aborting",quote=FALSE)
    quit(save="no")
  }

  if (dataType == "NUC") {    #Extracting data from listFile if nucleotide sequences

    if (ncol(geneList) == 1) {
      colnames(geneList)<<-"AlignmentFile"
    } else if (ncol(geneList) == 4) {
      colnames(geneList)<<-c("AlignmentFile","ModelNum","Gamma","GammaCats")
    } else {
      print("Gene list doesn't seem to be properly formatted. Make sure it's either (a) a 1-column CSV file (default sequence evol. model params), or (b) a 4-column CSV file (giving sequence evol. model params)",quote=FALSE)
      quit(save="no")
    }

  } else if (dataType == "AA") {    #Extracting data from listFile if amino acid sequences

    if (ncol(geneList) == 1) {
      colnames(geneList)<<-"AlignmentFile"
    } else if (ncol(geneList) == 4) {
      colnames(geneList)<<-c("AlignmentFile","SubstMatrix","Gamma","GammaCats")
      AAmatricesList<<-unique(geneList$SubstMatrix)

      if (length(grep("FALSE",AAmatricesList %in% fileList))) {
        print("Some subst. matrix file(s) given were not found in the current directory. Aborting.",quote=FALSE)
        quit()
      } else {
        print("Subst. matrix file(s) found in current directory. Proceeding.",quote=FALSE)
      }
    } else {
      print("Gene list doesn't seem to be properly formatted. Make sure it's either (a) a 1-column CSV file (default subst. matrix and model params), or (b) a 4-column CSV file (giving subst. matrix and model params.)",quote=FALSE)
      quit(save="no")
    }
  }

aliPresence<-geneList$AlignmentFile %in% fileList      #check all input alignments are there
  if (length(grep("FALSE",aliPresence))) {
    print("The following alignments were not found:",quote=FALSE)
    print(geneList$AlignmentFile[grep("FALSE",aliPresence)],quote=FALSE)
    print("Aborting.",quote=FALSE)
    quit()
  } else {
    print("All alignments were found in the current directory. Proceeding.",quote=FALSE)
  }

} #End of the readInData function



trimAlignment<-function(inpAli,living.sp) {   #Reduces alignment to only the species given in character vector living.sp
    inpAli[match(living.sp,labels(inpAli)[[1]]),]
} #End of the trimAlignment function



trimTree<-function(inpTree,living.sp) {       #Reduces the tree to only the species given in character vector living.sp
  droppedTips<-inpTree$tip.label[grep("FALSE",inpTree$tip.label %in% living.sp)]
  ladderize(drop.tip(inpTree,droppedTips,trim.internal=TRUE),right=TRUE)
} #End of the trimTree function



getRidOfEmptySeqs<-function(inpAli) {         #Removes the sequences that are composed of only missing data
  geneSpeciesList<-labels(inpAli)[[1]]
  seqEmptinessTable<-sapply(geneSpeciesList,calcPropMissingAndGapData,inpAli)
  emptySeqs<-labels(seqEmptinessTable[ seqEmptinessTable == 1 ])
  if (length(emptySeqs)) {
    nonEmptySeqs<-geneSpeciesList[ !geneSpeciesList %in% emptySeqs ]
    return(trimAlignment(inpAli,nonEmptySeqs))
  } else {
    return(inpAli)
  }
} #End of the getRidOfEmptySeqs function



calcPropMissingAndGapData<-function(sppName,inpAli) {    #Calculates for each species the proportion of their data coded as missing or gap
  if (dataType == "NUC") {
    totMissingSites<-length(grep("\\?",inpAli[sppName,])) + length(grep("-",inpAli[sppName,])) + length(grep("n",inpAli[sppName,])) + length(grep("N",inpAli[sppName,]))
  } else if (dataType == "AA") {
    totMissingSites<-length(grep("\\?",inpAli[sppName,])) + length(grep("-",inpAli[sppName,])) + length(grep("x",inpAli[sppName,])) + length(grep("X",inpAli[sppName,]))
  }
  totSites<-length(inpAli[sppName,])
  propMissing<-totMissingSites/totSites
  return(propMissing)
} #End of the calcPropMissingAndGapData function



makePamlAli<-function(inpAli,outFilePrefix) {          #Makes alignments "PAML friendly" (adds an I at the end of the first line)
  outFileName<-paste(outFilePrefix,".phy",sep="")
  sppNumber<-nrow(inpAli)
  sitesNumber<-ncol(inpAli)
  header<-paste(sppNumber," ",sitesNumber," I",sep="")

  write.dna(toupper(inpAli),file=outFileName,format="interleaved")
  speciesList<-labels(inpAli)[[1]]

  nonPamlFriendly<-readLines(outFileName,-1)
  pamlFriendly<-sub(head(nonPamlFriendly,1),header,nonPamlFriendly)
  for (s in speciesList) {
    pamlFriendly<-sub(s,paste(s," ",sep=""),pamlFriendly)
  }

  writeLines(pamlFriendly,outFileName)
} #End of the makePamlAli function



makePamlTrees<-function(inpTree,outFilePrefix,currSppCl) {      #Labels the branches of the trees for the LRT analysis
  write.tree(inpTree,file=paste(outFilePrefix,"_backbone.nwk",sep=""))
  inpTree<-makeNodeLabel(inpTree,method="number",prefix="Node")
  nodeParts<-treePart(inpTree)
  write.tree(inpTree,"protoTree.nwk")
  protoTree<-readLines("protoTree.nwk",-1)
  inpTreeSR<-protoTree
  inpTreeMR<-protoTree
  lineageNum<-1
  for (h in levels(currSppCl$Lineage)[grep("Clade",levels(currSppCl$Lineage))]) {
    cladeSppList<-currSppCl$Species[grep(h,currSppCl$Lineage)]
    cladeSppNum<-length(cladeSppList)

    if (cladeSppNum == 1) {
      inpTreeSR<-sub(as.character(cladeSppList[1]),paste(as.character(cladeSppList[1]),"$1",sep=""),inpTreeSR)
      inpTreeMR<-sub(as.character(cladeSppList[1]),paste(as.character(cladeSppList[1]),"$",lineageNum,sep=""),inpTreeMR)
    } else {
      cladeNode<-colnames(nodeParts)[max(grep(cladeSppNum,colSums(nodeParts[labels(nodeParts)[[1]] %in% cladeSppList,])))]

      inpTreeSR<-sub(paste(cladeNode,",",sep=""),"$1,",inpTreeSR)       #Each Node name can be before a comma, a closing parenthesis ")", or a semicolon ";" (although technically if it's the semicolon, the tree is not properly rooted).
      inpTreeSR<-sub(paste(cladeNode,")",sep=""),"$1)",inpTreeSR)
      inpTreeSR<-sub(paste(cladeNode,";",sep=""),"$1;",inpTreeSR)
      inpTreeMR<-sub(paste(cladeNode,",",sep=""),paste("$",lineageNum,",",sep=""),inpTreeMR)
      inpTreeMR<-sub(paste(cladeNode,")",sep=""),paste("$",lineageNum,")",sep=""),inpTreeMR)
      inpTreeMR<-sub(paste(cladeNode,";",sep=""),paste("$",lineageNum,";",sep=""),inpTreeMR)
    }
    lineageNum<-lineageNum+1
  }

  inpTreeSR<-gsub("Node[0-9]*","",inpTreeSR)
  inpTreeMR<-gsub("Node[0-9]*","",inpTreeMR)
  writeLines(inpTreeSR,paste(outFilePrefix,"_SR.nwk",sep=""))
  writeLines(inpTreeMR,paste(outFilePrefix,"_MR.nwk",sep=""))
  unlink("protoTree.nwk") 
} #End of the makePamlTrees function



makeBasemlCtrlFiles<-function(baseFile,outFilePrefix,modelNum,gamma,gammaRateCategories) {   #Prepares the files for the Baseml analysis (for nucleotide data)
  aliFile<-paste(outFilePrefix,".phy",sep="")
  BLOtreeFile<-paste(outFilePrefix,"_backbone.nwk",sep="")
  sRLRTtreeFile<-paste(outFilePrefix,"_SR.nwk",sep="")
  mRLRTtreeFile<-paste(outFilePrefix,"_MR.nwk",sep="")

  baseFile[grep("seqfile =",baseFile)]<-paste("      seqfile = ",aliFile,sep="")

  baseFile<-sub("model = [0-9]",paste("model = ",modelNum,sep=""),baseFile)
  baseFile<-sub("RateAncestor = [0-9]","RateAncestor = 0",baseFile)

  if (gamma) {
    baseFile<-sub("fix_alpha = [0-9]","fix_alpha = 0",baseFile)
    baseFile[grep("ncatG = ",baseFile)]<-paste("\tncatG = ",gammaRateCategories,sep="")
  } else {
    baseFile<-sub("fix_alpha = [0-9]","fix_alpha = 1",baseFile)
    baseFile<-sub(" alpha = [0-9].?[0-9]*"," alpha = 0",baseFile)
    baseFile<-sub("        ncatG = [0-9]{1,3}","*       ncatG = 8",baseFile)
  }

  baseFileBLO<-sub("clock = [0-9]","clock = 0",baseFile)
  baseFileBLO[grep("outfile = ",baseFileBLO)]<-paste("      outfile = ",outFilePrefix,"_branchopt",sep="")
  baseFileBLO[grep("treefile =",baseFileBLO)]<-paste("     treefile = ",BLOtreeFile,sep="")

  baseFileLRTsR<-sub("clock = [0-9]","clock = 2",baseFile)
  baseFileLRTsR[grep("outfile = ",baseFileLRTsR)]<-paste("      outfile = ",outFilePrefix,"_SR",sep="") #outfile name for the "single rate" likelihood estimation
  baseFileLRTsR[grep("treefile =",baseFileLRTsR)]<-paste("     treefile = ",sRLRTtreeFile,sep="")

  baseFileLRTmR<-baseFileLRTsR
  baseFileLRTmR[grep("outfile =",baseFileLRTmR)]<-paste("      outfile = ",outFilePrefix,"_MR",sep="") #outfile name for the "mutiple rate" likelihood estimation
  baseFileLRTmR[grep("treefile =",baseFileLRTmR)]<-paste("     treefile = ",mRLRTtreeFile,sep="")

  BLObaseFileName<<-tempfile("baseml-branchopts_", tmpdir=".", fileext=".ctl")
  LRTsRbaseFileName<<-tempfile("baseml-SR_", tmpdir=".", fileext=".ctl")
  LRTmRbaseFileName<<-tempfile("baseml-MR_", tmpdir=".", fileext=".ctl")

  writeLines(baseFileBLO,BLObaseFileName)
  writeLines(baseFileLRTsR,LRTsRbaseFileName)
  writeLines(baseFileLRTmR,LRTmRbaseFileName)
} #End of the makeBasemlCtrlFiles function



makeAamlCtrlFiles<-function(aamlBaseFile,outFilePrefix,substMatrixFile,gamma,gammaRateCategories) {   #Prepares the files for the Codeml analysis (for amino acid data)
  aliFile<-paste(outFilePrefix,".phy",sep="")
  BLOtreeFile<-paste(outFilePrefix,"_backbone.nwk",sep="")
  sRLRTtreeFile<-paste(outFilePrefix,"_SR.nwk",sep="")
  mRLRTtreeFile<-paste(outFilePrefix,"_MR.nwk",sep="")

  aamlBaseFile[grep("seqfile =",aamlBaseFile)]<-paste("      seqfile = ",aliFile,sep="")
  aamlBaseFile[grep("aaRatefile =",aamlBaseFile)]<-paste("      aaRatefile = ",substMatrixFile,sep="")

  aamlBaseFile<-sub("model = [0-9]",paste("model = 2"," ",sep=""),aamlBaseFile)  #CONSIDER ALLOWING THE USER TO CHOOSE
  aamlBaseFile<-sub("RateAncestor = [0-9]","RateAncestor = 0",aamlBaseFile)      #CONSIDER ALLOWING THE USER TO CHOOSE

  if (gamma) {
    aamlBaseFile<-sub("fix_alpha = [0-9]","fix_alpha = 0",aamlBaseFile)
    aamlBaseFile[grep("ncatG = ",aamlBaseFile)]<-paste("\tncatG = ",gammaRateCategories,sep="")
  } else {
    aamlBaseFile<-sub("fix_alpha = [0-9]","fix_alpha = 1",aamlBaseFile)
    aamlBaseFile<-sub(" alpha = [0-9].?[0-9]*"," alpha = 0",aamlBaseFile)
    aamlBaseFile<-sub("        ncatG = [0-9]{1,3}","*       ncatG = 8",aamlBaseFile)
  }

  baseFileBLO<-sub("clock = [0-9]","clock = 0",aamlBaseFile)
  baseFileBLO[grep("outfile = ",baseFileBLO)]<-paste("      outfile = ",outFilePrefix,"_branchopt",sep="")
  baseFileBLO[grep("treefile =",baseFileBLO)]<-paste("     treefile = ",BLOtreeFile,sep="")

  baseFileLRTsR<-sub("clock = [0-9]","clock = 2",aamlBaseFile)
  baseFileLRTsR[grep("outfile = ",baseFileLRTsR)]<-paste("      outfile = ",outFilePrefix,"_SR",sep="") #outfile name for the "single rate" likelihood estimation
  baseFileLRTsR[grep("treefile =",baseFileLRTsR)]<-paste("     treefile = ",sRLRTtreeFile,sep="")

  baseFileLRTmR<-baseFileLRTsR
  baseFileLRTmR[grep("outfile =",baseFileLRTmR)]<-paste("      outfile = ",outFilePrefix,"_MR",sep="") #outfile name for the "mutiple rate" likelihood estimation
  baseFileLRTmR[grep("treefile =",baseFileLRTmR)]<-paste("     treefile = ",mRLRTtreeFile,sep="")

  BLObaseFileName<<-tempfile("aaml-branchopts_", tmpdir=".", fileext=".ctl")
  LRTsRbaseFileName<<-tempfile("aaml-SR_", tmpdir=".", fileext=".ctl")
  LRTmRbaseFileName<<-tempfile("aaml-MR_", tmpdir=".", fileext=".ctl")

  writeLines(baseFileBLO,BLObaseFileName)
  writeLines(baseFileLRTsR,LRTsRbaseFileName)
  writeLines(baseFileLRTmR,LRTmRbaseFileName)
} #End of the makeAamlCtrlFiles function



LRT<-function(currPrefix) {  #Runs the LRT between the single- and multiple-rates models
  pValue<-1

  while (pValue == 1) {
    LRT_SROutFileName<-paste(currPrefix,"_SR",sep="")
    LRT_MROutFileName<-paste(currPrefix,"_MR",sep="")

    system(paste(PamlExecutablePath," ",LRTsRbaseFileName,sep=""),wait=TRUE)
    system(paste(PamlExecutablePath," ",LRTmRbaseFileName,sep=""),wait=TRUE)

    LRT_SROutFile<-readLines(LRT_SROutFileName,-1)
    SR_lnL<-as.numeric(strsplit(gsub(" +"," ",LRT_SROutFile[grep("^lnL",LRT_SROutFile)]),split=" ",fixed=TRUE)[[1]][5])

    LRT_MROutFile<-readLines(LRT_MROutFileName,-1)
    MR_lnL<-as.numeric(strsplit(gsub(" +"," ",LRT_MROutFile[grep("^lnL",LRT_MROutFile)]),split=" ",fixed=TRUE)[[1]][5])

    pValue<-pchisq(abs(SR_lnL-MR_lnL)*2,df=degreesOfFreedom,lower.tail=FALSE)
  }

  ratesString<-strsplit(gsub(" +"," ",LRT_MROutFile[grep("rates for branches",LRT_MROutFile)]),split=" ")[[1]]
  lineageRates<-as.numeric(ratesString[5:length(ratesString)])
  lineageRatesProportionalToFirst<-lineageRates/head(lineageRates,n=1)
  ratesVect<-c(lineageRates,lineageRatesProportionalToFirst)

  #Organizing output into a list, and making a vector for the labels of the list
  titlesVectRate<-vector()
  titlesVectPropRate<-vector()

  for (lin in 1:(length(ratesVect)/2)) {
    titlesVectRate[length(titlesVectRate)+1]<-paste("Clade",lin,".Rate",sep="")
    titlesVectPropRate[length(titlesVectPropRate)+1]<-paste("Clade",lin,".PropRate",sep="")
  }

  titlesVect<-c(titlesVectRate,titlesVectPropRate)
  lineageRatesList<-as.list(setNames(ratesVect,titlesVect))
  lineageRatesList$P.Value<-pValue
  lineageRatesList$Sub.Sample<-currPrefix

  return(lineageRatesList)
} #End of the LRT function



chooseToFlag<-function(aliPrefix,currSppCl) { #Chooses which will be the taxon to flag. Note: it now has two flavors ("LS3" and "LS4"), but it is easily extensible.
  outFile<-paste(aliPrefix,"_branchopt",sep="")
  ingroupList<-currSppCl$Species[grep("Clade",currSppCl$Lineage)]
  treeWBLengths<-NULL

  while (is.null(treeWBLengths)) {
    system(paste(PamlExecutablePath," ",BLObaseFileName,sep=""),wait=TRUE)
    BLengthOutput<-readLines(outFile,-1)
    writeLines(tail(BLengthOutput[grep("^\\(",BLengthOutput)],n=1),paste(aliPrefix,"_branchopt.nwk",sep=""))
    treeWBLengths<-read.tree(file=paste(aliPrefix,"_branchopt.nwk",sep=""))
  }

  SOBTable<-as.data.frame(distRoot(treeWBLengths,method="patristic"))
  SOBTable$Species<-factor(labels(SOBTable)[[1]])
  row.names(SOBTable)<-NULL
  colnames(SOBTable)[1]<-"SOB"

  currSOBTable<-merge(currSppCl,SOBTable)

  cladesAboveMinTaxa<-labels(which(head(summary(currSOBTable$Lineage),n=-1) > minTaxa))
  flaggableSpp<-data.frame(Species=factor(),Lineage=factor(),SOB=numeric())

  minTable<-numeric()

  for (lin in labels(head(summary(currSOBTable$Lineage),n=-1))) {
    cladeSOB<-currSOBTable[currSOBTable$Lineage == lin,]
    cladeMin<-min(cladeSOB$SOB)
    minTable<-rbind(cladeMin,minTable)
    if (lin %in% cladesAboveMinTaxa) {
      flaggableSpp<-rbind(currSOBTable[currSOBTable$Lineage == lin,],flaggableSpp)
    }
  }
  if (Flavor == "LS3") { #the taxon to be flgged is determined by the flavor. other flavors can be added here.
    toFlag<-flaggableSpp$Species[which.max(flaggableSpp$SOB)]
  } else if (Flavor == "LS4") {
    rateAim<-max(minTable)
    flaggableSpp$diffToAim<-abs(flaggableSpp$SOB-rateAim)
    toFlag<-flaggableSpp$Species[which.max(flaggableSpp$diffToAim)]
  }

  return(toFlag)
} #End of the choseToFlag function



makeHomRateAlis<-function(bestSsRow,geneOutDataFrame,geneAli,homRatesGeneOutFile,generalSppCl) { #If a subsample of taxa was found in which all taxa evolve homogeneously, it produces an alignment with these species only.
  totSpp<-generalSppCl$Species
  labelsGeneAli<-labels(geneAli)[[1]]
  if (bestSsRow > 1) {
    flaggedTaxa<-geneOutDataFrame$Flagged.Excluded[2:bestSsRow]
    unFlaggedFinal<-labelsGeneAli[!labelsGeneAli %in% flaggedTaxa]
    homRateAli<-trimAlignment(geneAli,unFlaggedFinal)
    write.dna(toupper(homRateAli),file=homRatesGeneOutFile,format="interleaved")
    nameFixer<-readLines(homRatesGeneOutFile,-1)
    for (s in labels(homRateAli)[[1]]) {
      nameFixer<-sub(s,paste(s," ",sep=""),nameFixer)
    }
    writeLines(nameFixer,homRatesGeneOutFile)

    presAbsVect<-sapply(totSpp,presAbs,flaggedTaxa,unFlaggedFinal)
    return(as.list(setNames(presAbsVect,totSpp)))

  } else {
    write.dna(toupper(geneAli),file=homRatesGeneOutFile,format="interleaved")
    nameFixer<-readLines(homRatesGeneOutFile,-1)
    for (s in labelsGeneAli) {
      nameFixer<-sub(s,paste(s," ",sep=""),nameFixer)
    }
    writeLines(nameFixer,homRatesGeneOutFile)

    flaggedTaxa<-NULL
    unFlaggedFinal<-labelsGeneAli
    presAbsVect<-sapply(totSpp,presAbs,flaggedTaxa,unFlaggedFinal)
    return(as.list(setNames(presAbsVect,totSpp)))
  }
}



presAbs<-function(someSpp,flaggedTaxa,unFlagged) {
  if (someSpp %in% flaggedTaxa) {
    return(0) #0 if species was flagged and removed
  } else if (someSpp %in% unFlagged) {
    return(1) #1 if the species is present in final dataset
  } else {
    return(NA) #NA if the species was not present in the original dataset
  }
} #End of makeHomRatesAli function



LSx_main<-function(geneOfInt,geneList) {  #This is the main LSx function that coordinates all the tests

  if (length(geneList[grep(paste("^",geneOfInt,sep=""),geneList[,1]),]) == 1) {
    currEntry<-geneList[grep(paste("^",geneOfInt,sep=""),geneList[,1]),]
    geneAli<-read.dna(file=currEntry,as.character=TRUE)
    geneFileName<-currEntry
  } else if (length(geneList[grep(paste("^",geneOfInt,sep=""),geneList[,1]),]) == 4) {
    currEntry<-geneList[grep(paste("^",geneOfInt,sep=""),geneList[,1]),]
    geneAli<-read.dna(file=currEntry$AlignmentFile,as.character=TRUE)
    geneFileName<-currEntry$AlignmentFile
  }

  geneAli<-getRidOfEmptySeqs(geneAli)
  genePrefix<-unlist(strsplit(geneFileName,split=".",fixed=TRUE))[1]
  homRatesGeneOutFile<-paste(genePrefix,"_homRates",Flavor,".phy",sep="")
  geneOutFile<-paste(genePrefix,"_",Flavor,"Out.csv",sep="")
  outDirName<-paste(genePrefix,"_",Flavor,sep="")

  dir.create(outDirName)

  if (dataType == "NUC") {
    modEvol<-if (length(currEntry) > 1) {
             currEntry$ModelNum
           } else {
             "7"
           } #Setting model of nucleotide sequence evolution (Note that the default is GTR)
  } else if (dataType == "AA") {
  substMatrix<-if (length(currEntry) > 1) {   #if there's something more than just an alignment name in the listFile...
             currEntry$SubstMatrix            #read subst matrix from listFile
           } else {
             if (length(grep("defaultAAMatrix",ls(.GlobalEnv)))) {  #if defaultAAMatrix was not commented out
               print(paste("Using default amino acid substitution matrix: ",defaultAAMatrix),quote=FALSE)
               defaultAAMatrix
             } else {          #if no subst matrix is given in the genelist, and defaultAAMatrix was commented out
               print("You have not specified a gene-specific amino acid substitution matrix NOR the defaultAAMatrix variable in the LSx control file. Aborting.",quote=FALSE)
               quit()
             }
           } #Setting matrix of amino acid substitution (Note that the default is defaultAAMatrix - given in the input text file)
  }

  gamma<-if (length(currEntry) > 1) {
           currEntry$Gamma == "G"
         } else {
           TRUE
         } #Checkng whether to use gamma or not (note that if only the filename is given, it defaults to having gamma)
  gammaRateCats<-if (gamma) {
                   if (length(currEntry) > 1) {
                     currEntry$GammaCats
                   } else {
                     "4"
                   }
                 } else {
                   FALSE
                 } #Sets the number of gamma rate categories. Note that if no number of gamma rate categories is given, but gamma is present, the default is 4. If gamma is absent, the variable defaults to "FALSE".
  currSppClTmp<-generalSppCl[match(labels(geneAli)[[1]],generalSppCl$Species),]

  if (length(levels(currSppClTmp$Lineage)) == length(lineagesOfInt)+1 ) {
    allLineagesPresent<-TRUE
  } else {
    allLineagesPresent<-FALSE
    print(paste(geneFileName,": This gene does not seem to have all inner clades represented",sep=""),quote=FALSE)
  }

  options(stringsAsFactors = FALSE)
  geneOutDataFrame<-data.frame()
  currSubsample<--1
  toFlag<-""
  currAli<-geneAli

  resubsample<-length(grep("TRUE",summary(currSppClTmp$Lineage)[-length(summary(currSppClTmp$Lineage))] > minTaxa))

#SUBSAMPLING LOOP
  while (resubsample && allLineagesPresent)
  {
    currSubsample<-currSubsample+1
    unFlagged<-labels(currAli)[[1]][!labels(currAli)[[1]] %in% toFlag]
    currAli<-trimAlignment(currAli,unFlagged)
    currAliPrefix<-paste(unlist(strsplit(geneFileName,split=".",fixed=TRUE))[1],"_ss",currSubsample,sep="") #Name for the run
    currSpp<-labels(currAli)[[1]] #List of species in current subsample
    currSppCl<-generalSppCl[match(currSpp,generalSppCl$Species),] #Gene's Sppcl file
    currGuideTree<-if (length(guideTree$tip.label) > length(currSpp) ) {
                     trimTree(guideTree,currSpp)
                   } else {
                     guideTree
                   } #Gene and subsample's guide tree for paml (works like this for branch length estimation, must be modified with "$" for LRT
    makePamlAli(currAli,currAliPrefix) #Produce phy alignment for paml's analyses (both branch length estimation AND likelihood ratio test)
    makePamlTrees(currGuideTree,currAliPrefix,currSppCl)

    if (dataType == "NUC") {
      makeBasemlCtrlFiles(basemlCtrlFile,currAliPrefix,modEvol,gamma,gammaRateCats)
    } else if (dataType == "AA") {
      makeAamlCtrlFiles(aamlCtrlFile,currAliPrefix,substMatrix,gamma,gammaRateCats)
    }

    LRTResult<-LRT(currAliPrefix)
    toFlagCol<-list(Flagged.Excluded=toFlag)
    geneOutDataFrame<-rbind(geneOutDataFrame,c(LRTResult,toFlagCol))
    toFlag<-chooseToFlag(currAliPrefix,currSppCl)

    resubsample<-length(grep("TRUE",summary(currSppCl$Lineage)[-length(summary(currSppCl$Lineage))] > minTaxa))

    filesToPutAway<-list.files(path=".",pattern=currAliPrefix,full.names=TRUE)
    file.copy(filesToPutAway,outDirName)
    file.remove(c(filesToPutAway,BLObaseFileName,LRTsRbaseFileName,LRTmRbaseFileName))
  } #####End of LSx_main subsampling loop

  geneOutDataFrameColNums<-ncol(geneOutDataFrame)
  rownames(geneOutDataFrame)<-1:nrow(geneOutDataFrame)
  geneOutDataFrame<-geneOutDataFrame[,c(geneOutDataFrameColNums-1,1:(geneOutDataFrameColNums-2),geneOutDataFrameColNums)]
  bestSsRow<-as.numeric(rownames(head(geneOutDataFrame[geneOutDataFrame$P.Value >= pThresh,],n=1)))

  if (length(bestSsRow)) {
    genePresAbsResults<-makeHomRateAlis(bestSsRow,geneOutDataFrame,geneAli,homRatesGeneOutFile,generalSppCl)
    homFound=TRUE
  } else {
    sppNamesVect<-generalSppCl$Species
    zeroVect<-rep.int(0,times=length(sppNamesVect))
    genePresAbsResults<-as.list(setNames(zeroVect,sppNamesVect))
    homFound=FALSE
  }

  write.csv(geneOutDataFrame,file=geneOutFile,quote=FALSE)
  totResultsList<-list(Final.DF=totGenesPresAbsDF,Final.Flagged=flaggedGenes)
  return(genePresAbsResults)
} ##End of LSx_main function



masterFunction<-function() {  #Runs the entire test and organizes and outputs the resulting files. Will probably work interactively if both checkInputs() and readInData() run without error.

  totGenesPresAbsDF<<-data.frame()
  flaggedGenes<<-data.frame(Flagged.Gene=character(0))

  if (parallelization) {
    clusterExport(cluster,ls(.GlobalEnv))
    clusterEvalQ(cluster,library(ape))
    clusterEvalQ(cluster,library(adephylo))
    totResultsDF<<-parSapply(cluster,geneList$AlignmentFile,LSx_main,geneList)
    stopCluster(cluster)
  } else {
    totResultsDF<<-sapply(geneList$AlignmentFile,LSx_main,geneList)
  }
#  totResultsDF<<-sapply(geneList$AlignmentFile,LSx_main,geneList)

  print(totResultsDF)
  fullRunSuffix<<-gsub(":","-",gsub(" ","",date()))
  write.csv(t(totResultsDF),file=paste(Flavor,"_presAbsTot_",fullRunSuffix,".csv",sep=""),quote=FALSE)

  #warnings()  #Uncomment if you're suspecting a bug!
} #End of masterFunction

#End of all Functions



#####-------------------------------------##CODE LAUNCHER##----------------------------------#####



if ( !interactive() ) {
  arg=commandArgs(trailingOnly=TRUE) #For handling arguments properly

  suppressPackageStartupMessages(library(adephylo))
  suppressPackageStartupMessages(library(ape))

  PlatformOS<-.Platform$OS.type
  if (PlatformOS == "unix" || PlatformOS == "windows") {
    print(paste("OS detected: ",PlatformOS,sep=""),quote=FALSE)
  } else {
    print(paste("Sorry, but I don't know how to run this analysis on a ",PlatformOS," system. Aborting",sep=""),quote=FALSE)
    quit()
  }

  fileList<-list.files()
  source(arg[1]) #import variables

  checkInputs()
  readInData()
  masterFunction()
}
