bounds=function (nums, dup.inds = FALSE, pre.sort = TRUE) 
{
    if (dup.inds) {
        s.ind <- which(!duplicated(nums))
        e.ind <- c(s.ind[-1] - 1, length(nums))
        return(cbind(1:length(s.ind), start = s.ind, end = e.ind, 
            length = (e.ind - s.ind + 1)))
    }
    else {
        if (!is.numeric(nums)) 
            stop("must supply a numeric vector")
        # CHANGES ARE DONE HERE, ORIGINAL ==0, WE CHANGED <=0
        # When DSSP secondary st analysis gives output such as HHHHEEETGGG
        # Look at T, only one is there, bounds fails bcoz no start, end and length.
        if (length(nums) == 0 ) 
            return(nums)
        if (length(nums) == 1 ) 
            return(which(nums=="TP"))
        # We should learn how to get integer(0), here to get we have tricked but u must learn. 
        if (length(nums) <= 2) {
            bounds <- sort(nums)
            bounds <- c(bounds, (bounds[2] - bounds[1] + 1))
            names(bounds) <- c("start", "end", "length")
            return(t(as.matrix(bounds)))
        }
        if (pre.sort) {
            nums <- sort(unique(nums))
        }
        bounds <- nums[1]
        nums.start <- nums[1]
        diff.i <- 1
        for (i in 2:length(nums)) {
            if ((nums[i] - diff.i) != nums.start) {
                bounds <- c(bounds, nums[i - 1], nums[i])
                nums.start <- nums[i]
                diff.i <- 1
            }
            else {
                diff.i <- diff.i + 1
            }
        }
        bounds <- c(bounds, nums[length(nums)])
        bounds <- matrix(bounds, ncol = 2, byrow = TRUE, dimnames = list(c(1:(length(bounds)/2)), 
            c("start", "end")))
        bounds <- cbind(bounds, length = (bounds[, 2] - bounds[, 
            1]) + 1)
        return(bounds)
    }
}
dssp_new=function (pdb, exepath = "",prefix, resno = TRUE) 
{
    print("Running DSSP");
    infile=paste(prefix,".pdb",sep="",collapse="");
    outfile=paste(prefix,".dssp",sep="",collapse="");
    write.pdb(pdb, file = infile)
    # If dssp command gives any error, then due to intern=TRUE we get
    # a attributed variable. So the variable has an attribute which is a list
    # If error, then is.list(attributes(error)) returns true and we return sorry saying dssp is not running
    # If error, this function returns "SORRY" else a list.
    # We test is.list(DSSP) and then call any function dependent on DSSP result (E.G. ionic_addDSSP)
    tryCatch({error=system(paste(exepath, " -i ", infile, " -o ", outfile, sep = ""),ignore.stderr = FALSE,intern=TRUE)},error=function(e){return("SORRY");});
    
    if(is.list(attributes(error))){return("SORRY");}
    if(!file.exists(outfile))
    {
      return("SORRY");
    }
    
    raw.lines <- readLines(outfile)
    unlink(infile)
    type <- substring(raw.lines, 1, 3)
    raw.lines <- raw.lines[-(1:which(type == "  #"))]
    ter <- substring(raw.lines, 14, 15);
    ter.index=which(ter=="!*");
    if(length(ter.index)>0){raw.lines <-    raw.lines[-ter.index]}
#    670        !*
# Sometimes when residues are missing ! comes as we can see below for number 17, we have to remove these
#   16   19 A L     <        0   0    4     -4,-3.2     7,-2.0    -5,-0.2    -2,-0.2   0.886 360.0 360.0 -54.9 360.0   31.7   19.5   21.5
#   17        !              0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0
#   18   21 A G     >        0   0   40     -4,-2.0     4,-2.7    -5,-0.2     5,-0.2   0.000 360.0 360.0 360.0 167.2   27.4   18.3   18.2

    cha <- substring(raw.lines, 12, 12)
    aa <- substring(raw.lines,14,14);
    sse <- substring(raw.lines, 17, 17)
    sse[which(sse==" ")]="C";
    phi <- as.numeric(substring(raw.lines, 104, 109))
    psi <- as.numeric(substring(raw.lines, 110, 115))
    acc <- as.numeric(substring(raw.lines, 35, 38))
    seqres <- substring(raw.lines, 6, 11) # includes insert 67P where P=insertion code
seqres=gsub(" ","",seqres)#remove space


    h.ind <- bounds(which(sse == "H"), pre.sort = FALSE)
    g.ind <- bounds(which(sse == "G"), pre.sort = FALSE)
    e.ind <- bounds(which(sse == "E"), pre.sort = FALSE)
    t.ind <- bounds(which(sse == "T"), pre.sort = FALSE)
    h.res <- h.ind
    g.res <- g.ind
    e.res <- e.ind
    t.res <- t.ind
    if (resno) {
        res.num <- substring(raw.lines, 6, 11) # includes insert
	res.num=gsub(" ","",res.num) #remove space
        if (length(h.ind) > 0) {
            h.res[, "start"] <- res.num[h.ind[, "start"]]
            h.res[, "end"] <- res.num[h.ind[, "end"]]
        }
        if (length(g.ind) > 0) {
            g.res[, "start"] <- res.num[g.ind[, "start"]]
            g.res[, "end"] <- res.num[g.ind[, "end"]]
        }
        if (length(e.ind) > 0) {
            e.res[, "start"] <- res.num[e.ind[, "start"]]
            e.res[, "end"] <- res.num[e.ind[, "end"]]
        }
        if (length(t.ind) > 0) {
            t.res[, "start"] <- res.num[t.ind[, "start"]]
            t.res[, "end"] <- res.num[t.ind[, "end"]]
        }
    }
    sheet = list(start = NULL, end = NULL, length = NULL, chain = NULL)
    helix = list(start = NULL, end = NULL, length = NULL, chain = NULL, 
        type = NULL)
    turn = sheet
    if (length(h.res) > 1) {
        if (is.null(nrow(h.res))) 
            h.s <- as.matrix(t(h.res))
        helix$start = c(helix$start, h.res[, "start"])
        helix$end = c(helix$end, h.res[, "end"])
        helix$length = c(helix$length, h.res[, "length"])
        helix$chain = c(helix$chain, cha[h.ind[, "start"]])
        helix$type = c(helix$type, sse[h.ind[, "start"]])
    }
    if (length(g.res) > 1) {
        if (is.null(nrow(g.res))) 
            g.s <- as.matrix(t(g.res))
        helix$start = c(helix$start, g.res[, "start"])
        helix$end = c(helix$end, g.res[, "end"])
        helix$length = c(helix$length, g.res[, "length"])
        helix$chain = c(helix$chain, cha[g.ind[, "start"]])
        helix$type = c(helix$type, sse[g.ind[, "start"]])
    }
    if (length(e.res) > 1) {
        if (is.null(nrow(e.res))) 
            e.s <- as.matrix(t(e.res))
        sheet$start = c(sheet$start, e.res[, "start"])
        sheet$end = c(sheet$end, e.res[, "end"])
        sheet$length = c(sheet$length, e.res[, "length"])
        sheet$chain = c(sheet$chain, cha[e.ind[, "start"]])
    }
    if (length(t.res) > 1) {
        if (is.null(nrow(t.res))) 
            t.s <- as.matrix(t(t.res))
        turn$start = c(turn$start, t.res[, "start"])
        turn$end = c(turn$end, t.res[, "end"])
        turn$length = c(turn$length, t.res[, "length"])
        turn$chain = c(turn$chain, cha[t.ind[, "start"]])
    }
 print("Finished DSSP");
    out <- list(helix = helix, sheet = sheet, turn = turn, phi = phi, 
        psi = psi, acc = acc, aa=aa, res=seqres, cha=cha,ss=sse)
}

