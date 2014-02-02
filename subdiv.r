function(matrix, lag, rand, factor){
    #o is useful to see which bin the program is on. This value is set to 1 when no bins exist (ie. comparing the bin results to one another).
    #factor is the size of each bin, generally set to 4, but I would like to improve this to 5 after code optimization, this value will vary
    #for remainder solutions and for solving amongst different bins
    #seq is a single row from one of the three seq (normal bin), mseq (mod or remainder bin), or mainSeq (all samples) matrices. 
    #This is the top of the distance matrix [1,]
    #lag is self evident, its the lag you want to look at
    #dist, is the 3D distance array returned by ccf() aka. the autocorrelation function
    #random is the number of random iterations we want to do of the entire code, this important so we can output every solution made
    distance=function(o, factor, seq, lag, dist, random){
        #This first function is what turns a distance matrix with just the top row completed -> complete distance matrix -> correlation matrix -> avg correlation
        ret=NULL
        tempor=matrix(data=NA, ncol=factor, nrow=factor)
        if(o==1){
            seq=seq+-1*seq[1]
        }
        #this is to adjust for the fact that [1,1] needs to be set to zero in the distance matrix
        #this only becomes a problem with the global solution looking between bins
        tempor[1,]=seq
        tempor[,1]=-1*seq
        for(i in 1:factor){
            tempor[i,i]=0
        }
        #flag is to see if a value is changed to another value
        flag=0
        #This fills in all the unknown values through subtraction
        for(i in 1:factor){
            for(j in 1:factor){
                for(k in 1:factor){
                    if(!is.na(tempor[i,j])&&!is.na(tempor[i,k])){  
                        if(is.na(tempor[j,k])&&!(k==i || j == i || k == j)){
                            tempor[k,j]=tempor[i,j]-tempor[i,k]
                            tempor[j,k]=tempor[i,k]-tempor[i,j]
                        }
                        else if(tempor[k,j]!=tempor[i,j]-tempor[i,k]&&!(k==i || j == i || k == j)){
                            flag=flag+1
                        }
                    }
                }  
            }
        }
        if(any(is.na(tempor))){
            ret$avg=-3
            return(ret)
        }
        if(flag>0){
            ret$avg=-1
            return(ret)
        }
        summary=matrix(data=NA, nrow=factor, ncol=factor)
        for(v in 1:factor){
            for(j in 1:factor){
                if((tempor[v,j]+2*lag+1)>(4*lag+1)||(tempor[v,j]+2*lag+1)<1){
                    #This is in case the answer is out of subscript bounds or NA, we need to not try and solve this, it will crash the program
                    ret$avg=-4
                    ret$temp=tempor
                    return(ret)
                }
                else if(is.na(dist[(tempor[v,j]+2*lag+1),random[((o-1)*factor+v)],random[((o-1)*factor+j)]])){
                    ret$avg=-5
                    ret$temp=tempor
                    return(ret)
                }
                summary[v,j]=dist[(tempor[v,j]+2*lag+1),random[((o-1)*factor+v)],random[((o-1)*factor+j)]]
            }
        }
        if(any(is.na(summary))){
            ret$avg=-2
            return(ret)
        }
        ret$avg=sum(summary)/factor^2
        ret$tempor=tempor
        return(ret)
    }
    
    ##############################################################
    
    matrix=as.matrix(matrix)
    if(length(matrix[,1])<length(matrix[1,])){
        matrix=t(matrix)
    }
    wid=length(matrix[1,])
    hei=length(matrix[,1])
    #don't forget {pracma} library in R
    #gets out weird csv formatting "V1 V2 V3..."
    trix=matrix(data=NA, ncol=length(matrix[1,]), nrow=length(matrix[,1]))
    #detrending here
    for(i in 1:wid){
        temp=detrend(matrix[,i][!is.na(matrix[,i])])
        trix[1:length(temp),i]=temp
    }
    #return(trix)
    #The lag k value returned by ccf(x, y) estimates the correlation between x[t+k] and y[t].
    
    #changed coefficient to 4 so that it can accommodate the extra lagging so its harder for us to go out of subscript bounds
    
    dist=array(NA, dim=c(4*lag+1, wid, wid))
    for(i in 1:wid){
        for(j in 1:wid){
            temp=ccf(trix[,i], trix[,j], na.action=na.pass, lag.max=lag*2, plot=FALSE)
            dist[,j,i]=temp$acf
        }
    }
    slag=lag
    #lag=as.integer(lag/2)
    #I'm reducing lag here, so the software has a lower chance of going out of range when solving in the distance() function
    tempor=matrix(data=NA, nrow=wid, ncol=wid)
    #tempor is the distance matrix used by the distance() function
    prime=(lag*2+1)
    #prime is used in calculating how many iterations for building sequence matrices and for going through all possibilities
    rema=wid%%factor
    #rema stands for remainder
    if(rema){
        mfac=as.integer(wid/factor)+1
        #mfac is the number of bins to be solved
    }
    if(!rema){
        mfac=as.integer(wid/factor)
    }
    seq=matrix(data=-1*lag, ncol=factor, nrow=(prime^(factor-1)))
    if(rema){
        seqm=matrix(data=-1*lag, ncol=rema, nrow=(prime^(rema-1)))
        #You will see one of these for loops with similar code between seq (normal bin), mseq (mod or remainder bin), or mainSeq (all samples) matrices
        #They essentially are a clock, and the base number for each place or element is 2*lag+1
        for(t in 1:(prime^(rema-1)-1))
        {
            power=0
            for(q in 1:rema){
                if(t%%prime^q!=0){
                    break
                }
                power=power+1
            }
            seqm[t+1,]=seqm[t,]
            seqm[t+1,power+1]=seqm[t,power+1]+1
            if(power){
                for(jay in 1:power){
                    seqm[t+1,jay]=-1*lag
                }
            }
        }
        seqm[,2:rema]=seqm[,1:(rema-1)]
        seqm[,1]=0
    }
    
    seq[1,]=-1*lag
    for(t in 1:((prime^(factor-1))-1))
    {
        power=0
        for(q in 1:factor){
            if(t%%prime^q!=0){
                break
            }
            power=power+1
        }
        seq[t+1,]=seq[t,]
        seq[t+1,power+1]=seq[t,power+1]+1
        if(power){
            for(jay in 1:power){
                seq[t+1,jay]=-1*lag
            }
        }
        
    }
    seq[,2:factor]=seq[,1:(factor-1)]
    seq[,1]=0
    mainSeq=matrix(data=0, ncol=wid, nrow=prime^mfac)
    mainSeq[1,]=-1*lag
    for(t in 1:(prime^mfac-1)){
        power=0
        for(q in 1:mfac){
            if(t%%prime^q!=0){
                break
            }
            power=power+1
        }
        mainSeq[t+1,]=mainSeq[t,]
        #to differentiate if this is is a normal sized bin or the remainder bin
        if((factor*power+factor)>wid&&factor*power+1<=wid){
            mainSeq[t+1,(factor*power+1):wid]=mainSeq[t+1,(factor*power+1):wid]+1
            mainSeq[t+1,1:(wid-rema)]=-1*lag
        }
        else if((factor*power+factor)<=wid){
            mainSeq[t+1, (factor*power+1):(factor*power+factor)]=mainSeq[t+1, (factor*power+1):(factor*power+factor)]+1
            if(power){
                for(jay in 1:power){
                    mainSeq[t+1,(factor*(jay-1)+1):(jay*factor)]=-1*lag
                }
            }
        }
    }
    locAvg=0
    #loc(al)Avg is used to find the highest solution within each bin
    locTemp=NULL
    locSum=NULL
    gloVec=matrix(data=NA, ncol=wid, nrow=rand)
    gloAvg=matrix(data=0,nrow=rand,ncol=1)
    gloTemp=matrix(data=NA, nrow=wid, ncol=wid)
    #glo(bal) is used to store succesful solutions to be used for the alignment of the bins amongst each other
    truVec=matrix(data=NA, nrow=rand, ncol=wid)
    #tru(e)Vec(tor) is used to output the final solutions
    for(r in 1:rand)
    {
        print("~ROTATION~")
        random=sample(1:wid, wid, replace=FALSE)
        for(w in 1:mfac){
            locAvg=0
            if(w*factor<=wid){
                for(f in 1:(prime^(factor-1))){
                    temp=distance(w,factor,seq[f,],slag,dist,random)
                    if(temp$avg>locAvg){
                        print(temp$avg)
                        for(q in 1:factor){
                            gloVec[r,(q+(w-1)*factor)]=temp$tempor[1,q]
                        }
                        locAvg=temp$avg
                        locTemp=temp$tempor
                    }
                }
                print(locAvg)
            }
            else if(w*factor>wid&&rema!=1){
                for(f in 1:(prime^(rema-1))){
                    temp=distance(w,rema,seqm[f,],slag,dist,random)
                    if(temp$avg>locAvg){   
                        for(q in 1:rema){
                            gloVec[r,(q+(w-1)*factor)]=temp$tempor[1,q]
                        }
                        locAvg=temp$avg
                        locTemp=temp$tempor
                    }
                }
                print(locAvg)
            }
        }
        print(gloVec)
        if(rema==1){
            gloVec[r,random[mfac]]=0
        }
        #I need to adjust the main sequence matrix to adjust for the local solutions from within each bin, we want to keep the constrained distances
        #that are optimal within each bin set when comparing bins among each other
        for(o in 1:length(mainSeq[,1])){
            mainSeq[o,]=mainSeq[o,]+gloVec[r,]
        }
        locAvg=0
        mainSeq[,1]=0
        for(w in 1:(prime^mfac-1)){
            temp=distance(1,wid,mainSeq[w,],slag,dist,random)
            if(temp$avg>locAvg){
                gloVec[r,]=temp$tempor[1,]
                gloAvg[r]=temp$avg
                gloTemp=temp$tempor
                locAvg=temp$avg
                print(temp$avg)
            }
        }
        for(j in 1:wid){
            truVec[r,random[j]]=gloVec[r,j]
            #This removes the fact [1,] and [,1] or [n,] and [,n] doesn't refer to that sample number, it sorts them out for printing so that
            #you can easily identify which sample should be at which alignment
        }
    }    
    random=1:wid
    for(i in 1:rand){
        print("******************************************")
        temp=distance(1,wid,truVec[i,],slag,dist,random)
        print(temp)
        print("******************************************")
    }
}
